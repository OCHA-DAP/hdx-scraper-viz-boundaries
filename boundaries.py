import logging
import re
from geojson import loads
from geopandas import GeoDataFrame, read_file
from glob import glob
from json import dump
from mapbox import Uploader
from os import remove
from os.path import join
from pandas import concat, merge
from pandas.api.types import is_numeric_dtype
from shapely.geometry import MultiPolygon, box
from shapely.validation import make_valid
from time import sleep
from topojson import Topology
from unicodedata import normalize
from zipfile import BadZipFile, ZipFile

from hdx.data.dataset import Dataset
from hdx.data.hdxobject import HDXError
from hdx.location.country import Country
from hdx.utilities.downloader import DownloadError
from hdx.utilities.uuid import get_uuid

logger = logging.getLogger()


def drop_fields(df, keep_fields):
    df = df.drop(
        [f for f in df.columns if f not in keep_fields and f.lower() != "geometry"],
        axis=1,
    )
    return df


def replace_mapbox_tileset(
    mapid, mapbox_auth, name, json_to_upload, temp_folder
):
    service = Uploader(access_token=mapbox_auth)
    saved_file = join(temp_folder, "file_to_upload.geojson")
    json_to_upload.to_file(saved_file, driver="GeoJSON")
    with open(saved_file, "rb") as src:
        upload_resp = service.upload(src, mapid, name=name)
    if upload_resp.status_code == 422:
        for i in range(5):
            sleep(5)
            with open(saved_file, "rb") as src:
                upload_resp = service.upload(src, mapid, name=name)
            if upload_resp.status_code != 422:
                break
    if upload_resp.status_code == 422:
        logger.error(f"Could not upload {name}")
        return None
    return mapid


class Boundaries:
    def __init__(
        self, configuration, downloader, all_boundaries, mapbox_auth, temp_folder
    ):
        self.downloader = downloader
        self.boundaries = all_boundaries
        self.temp_folder = temp_folder
        self.exceptions = {
            "dataset": configuration["hdx_inputs"].get("dataset_exceptions", {}),
            "resource": configuration["hdx_inputs"].get("resource_exceptions", {}),
        }
        self.headers = configuration["shapefile_attribute_mappings"]
        self.mapbox = configuration["mapbox"]
        self.mapbox_auth = mapbox_auth
        self.countries = configuration["countries"]
        self.regional_info = configuration["regional"]
        self.bboxes = {}
        self.lookups = {}

    def update_subnational_boundaries(self, countries, do_not_process):
        for level in countries:
            req_fields = ["alpha_3"]
            for i in range(int(level[-1]) + 1):
                req_fields.append(f"ADM{i}_PCODE")
                req_fields.append(f"ADM{i}_REF")

            for iso in countries[level]:
                if iso in do_not_process:
                    logger.warning(f"{iso}: Not processing for now")
                    continue

                logger.info(f"{iso}: Processing {level} boundaries")

                # select single country boundary (including disputed areas), cut out water, and dissolve
                country_adm0 = self.boundaries["polbnda_int_1m"].copy(deep=True)
                country_adm0 = country_adm0.loc[country_adm0["ISO_3"] == iso]
                country_adm0 = country_adm0.overlay(self.boundaries["lake"], how="difference")
                country_adm0 = country_adm0.dissolve()
                country_adm0 = drop_fields(country_adm0, ["ISO_3"])
                country_adm0["ISO_3"] = iso
                if not country_adm0.crs:
                    country_adm0 = country_adm0.set_crs(crs="EPSG:4326")
                if not country_adm0["geometry"][0].is_valid:
                    country_adm0.loc[0, "geometry"] = make_valid(country_adm0.loc[0, "geometry"])

                # find the correct admin boundary dataset
                dataset_name = self.exceptions["dataset"].get(iso, f"cod-em-{iso.lower()}")
                dataset = Dataset.read_from_hdx(dataset_name)
                if not dataset:
                    dataset = Dataset.read_from_hdx(f"cod-ab-{iso}")
                if not dataset:
                    logger.error(f"{iso}: Could not find boundary dataset")
                    continue

                resource_name = self.exceptions["resource"].get(iso, "adm")
                boundary_resource = [
                    r
                    for r in dataset.get_resources()
                    if r.get_file_type() == "shp"
                    and bool(re.match(f".*{resource_name}.*", r["name"], re.IGNORECASE))
                ]
                if len(boundary_resource) > 1:
                    name_match = [
                        bool(re.match(f".*adm(in)?(\s)?(0)?{level[-1]}.*", r["name"], re.IGNORECASE))
                        for r in boundary_resource
                    ]
                    boundary_resource = [
                        boundary_resource[i]
                        for i in range(len(boundary_resource))
                        if name_match[i]
                    ]

                    if len(boundary_resource) != 1:
                        logger.error(f"{iso}: Could not distinguish between resources")
                        continue

                # find the correct admin boundary shapefile in the downloaded zip
                try:
                    _, resource_file = boundary_resource[0].download(folder=self.temp_folder)
                except DownloadError:
                    logger.error(f"{iso}: Could not download resource")
                    return None

                temp_folder = join(self.temp_folder, get_uuid())
                try:
                    with ZipFile(resource_file, "r") as z:
                        z.extractall(temp_folder)
                except BadZipFile:
                    logger.error(f"{iso}: Could not unzip file - it might not be a zip!")
                    continue
                boundary_shp = glob(join(temp_folder, "**", "*.shp"), recursive=True)

                if len(boundary_shp) == 0:
                    logger.error(f"{iso}: Did not find the file!")
                    continue

                if len(boundary_shp) > 1:
                    name_match = [
                        bool(re.match(f".*admbnda.*adm(in)?(0)?{level[-1]}.*", b, re.IGNORECASE))
                        for b in boundary_shp
                    ]
                    if any(name_match):
                        boundary_shp = [boundary_shp[i] for i in range(len(boundary_shp)) if name_match[i]]

                if len(boundary_shp) > 1:
                    simp_match = [bool(re.match(".*simplified.*", b, re.IGNORECASE)) for b in boundary_shp]
                    if any(simp_match):
                        boundary_shp = [boundary_shp[i] for i in range(len(boundary_shp)) if not simp_match[i]]

                if len(boundary_shp) != 1:
                    logger.error(f"{iso}: Could not distinguish between downloaded shapefiles")
                    continue

                boundary_lyr = read_file(boundary_shp[0])
                if not boundary_lyr.crs:
                    boundary_lyr = boundary_lyr.set_crs(crs="EPSG:4326")
                if not boundary_lyr.crs.name == "WGS 84":
                    boundary_lyr = boundary_lyr.to_crs(crs="EPSG:4326")

                # calculate fields, finding admin1 name and pcode fields from config
                boundary_lyr["alpha_3"] = iso.upper()
                boundary_lyr["ADM0_REF"] = Country.get_country_name_from_iso3(iso)
                boundary_lyr["ADM0_PCODE"] = Country.get_iso2_from_iso3(iso)

                fields = boundary_lyr.columns
                for l in range(1, int(level[-1]) + 1):
                    pcode_field = None
                    name_field = None
                    if f"ADM{l}_PCODE" in fields:
                        pcode_field = f"ADM{l}_PCODE"
                    if f"ADM{l}_EN" in fields:
                        name_field = f"ADM{l}_EN"
                    for field in fields:
                        if not pcode_field:
                            if field.upper() in self.headers["pcode"][f"adm{l}"]:
                                pcode_field = field
                        if not name_field:
                            if field.upper() in self.headers["name"][f"adm{l}"]:
                                name_field = field

                    if not name_field:
                        logger.error(f"{iso}: Could not map name field")
                        continue

                    # calculate text pcodes
                    # if pcod is not in field name or there are no codes in the boundaries, create pcodes
                    if not pcode_field:
                        boundary_lyr[f"ADM{l}_PCODE"] = ""
                    if pcode_field:
                        if is_numeric_dtype(boundary_lyr[pcode_field]):
                            boundary_lyr[f"ADM{l}_PCODE"] = (
                                boundary_lyr[pcode_field].astype(int).astype(str)
                            )
                        else:
                            boundary_lyr[f"ADM{l}_PCODE"] = boundary_lyr[pcode_field]

                    boundary_lyr[f"ADM{l}_REF"] = boundary_lyr[name_field]

                boundary_lyr = drop_fields(boundary_lyr, req_fields)
                boundary_lyr = boundary_lyr.dissolve(by=req_fields, as_index=False)

                # assign pcodes to highest level if they're missing
                if not pcode_field:
                    logger.error(f"{iso}: Could not map pcodes at {level} - assigning randomly!")
                    numrows = len(str(len(boundary_lyr.index)))
                    for i, _ in boundary_lyr.iterrows():
                        boundary_lyr.loc[i, f"ADM{level[-1]}_PCODE"] = boundary_lyr.loc[
                            i, "ADM0_PCODE"
                        ] + str(i + 1).zfill(numrows)

                na_count = boundary_lyr[f"ADM{level[-1]}_REF"].isna().sum()
                if na_count > 0:
                    logger.warning(f"{iso}: Found {na_count} null values at {level}")

                # simplify geometry of boundaries
                boundary_topo = Topology(boundary_lyr)
                eps = 0.0075
                if int(level[-1]) > 0:
                    eps = eps / int(level[-1])
                boundary_topo = boundary_topo.toposimplify(
                    epsilon=eps,
                    simplify_algorithm="dp",
                    prevent_oversimplify=True,
                )
                boundary_lyr = boundary_topo.to_gdf(crs="EPSG:4326")

                # make sure geometry is valid
                for i, _ in boundary_lyr.iterrows():
                    if not boundary_lyr.geometry[i].is_valid:
                        boundary_lyr.geometry[i] = make_valid(boundary_lyr.geometry[i])
                    if boundary_lyr.geometry[i].geometryType() == "GeometryCollection":
                        new_geom = []
                        for part in boundary_lyr.geometry[i]:
                            if part.geometryType() in ["Polygon", "MultiPolygon"]:
                                new_geom.append(part)
                        if len(new_geom) == 0:
                            logger.error(f"{iso}: Boundary found with no geometry")
                        if len(new_geom) == 1:
                            new_geom = new_geom[0]
                        else:
                            new_geom = MultiPolygon(new_geom)
                        boundary_lyr.geometry[i] = new_geom

                # clip international boundary to UN admin0 country boundary
                boundary_lyr = boundary_lyr.clip(mask=country_adm0, keep_geom_type=True)

                self.boundaries[f"polbnda_{level}"] = self.boundaries[f"polbnda_{level}"][
                    self.boundaries[f"polbnda_{level}"]["alpha_3"] != iso
                ]
                self.boundaries[f"polbnda_{level}"] = self.boundaries[f"polbnda_{level}"].append(boundary_lyr)
                self.boundaries[f"polbnda_{level}"].sort_values(by=[f"ADM{level[-1]}_PCODE"], inplace=True)
                logger.info(f"{iso}: Finished processing {level} boundaries")

            # convert polygon boundaries to point
            centroid = GeoDataFrame(self.boundaries[f"polbndp_{level}"].representative_point())
            centroid.rename(columns={0: "geometry"}, inplace=True)
            centroid[req_fields] = self.boundaries[f"polbndp_{level}"][req_fields]
            centroid = centroid.set_geometry("geometry")
            self.boundaries[f"polbndp_{level}"] = centroid

    def update_subnational_resources(self, dataset, levels):
        for level in levels:
            logger.info(f"Updating HDX datasets at {level}")
            polbnda_file = join(self.temp_folder, f"polbnda_{level}_1m_ocha.geojson")
            self.boundaries[f"polbnda_{level}"].to_file(polbnda_file, driver="GeoJSON")
            centroid_file = join(self.temp_folder, f"polbndp_{level}_1m_ocha.geojson")
            self.boundaries[f"polbndp_{level}"].to_file(centroid_file, driver="GeoJSON")

            resource_a = [
                r
                for r in dataset.get_resources()
                if r.get_file_type() == "geojson"
                and bool(re.match(f".*polbnda_{level}.*", r["name"], re.IGNORECASE))
            ][0]
            resource_a.set_file_to_upload(polbnda_file)
            resource_p = [
                r
                for r in dataset.get_resources()
                if r.get_file_type() == "geojson"
                and bool(re.match(f".*polbndp_{level}.*", r["name"], re.IGNORECASE))
            ][0]
            resource_p.set_file_to_upload(centroid_file)

            try:
                resource_a.update_in_hdx()
            except HDXError:
                logger.exception("Could not update polygon resource")
            try:
                resource_p.update_in_hdx()
            except HDXError:
                logger.exception("Could not update point resource")

    def merge_subn_boundaries(self, visualization):
        levels = [key for key in self.countries[visualization] if not key == "adm0"]
        merged_boundaries = dict()
        merged_boundaries["a"] = {
            level: self.boundaries[f"polbnda_{level}"].copy(deep=True) for level in levels
        }
        merged_boundaries["p"] = {
            level: self.boundaries[f"polbndp_{level}"].copy(deep=True)
            for level in levels
        }
        for t in merged_boundaries:
            for level in levels:
                merged_boundaries[t][level] = merged_boundaries[t][level][
                    merged_boundaries[t][level]["alpha_3"].isin(self.countries[visualization][level])
                ]
                merged_boundaries[t][level]["ADM_LEVEL"] = int(level[-1])
                merged_boundaries[t][level]["ADM_PCODE"] = merged_boundaries[t][level][f"ADM{level[-1]}_PCODE"]
                merged_boundaries[t][level]["ADM_REF"] = merged_boundaries[t][level][f"ADM{level[-1]}_REF"]
        merged_polygons = [merged_boundaries["a"][level] for level in levels]
        merged_points = [merged_boundaries["p"][level] for level in levels]
        merged_polygons = concat(merged_polygons)
        merged_points = concat(merged_points)
        if visualization == "hornafrica":
            merged_polygons.loc[merged_polygons["ADM_LEVEL"] == 3, "ADM_LEVEL"] = 2
            merged_points.loc[merged_points["ADM_LEVEL"] == 3, "ADM_LEVEL"] = 2
            excep = merged_polygons.copy(deep=True)
            excep = excep[excep["alpha_3"] == "KEN"]
            excep["ADM_LEVEL"] = 2
            merged_polygons = concat([merged_polygons, excep])
            excep = merged_points.copy(deep=True)
            excep = excep[excep["alpha_3"] == "KEN"]
            excep["ADM_LEVEL"] = 2
            merged_points = concat([merged_points, excep])
        return merged_polygons, merged_points

    def update_mapbox_tilesets(self, visualizations):
        for visualization in visualizations:
            logger.info(f"Updating MapBox tilesets for {visualization}")
            polygon_to_upload, point_to_upload = self.merge_subn_boundaries(visualization)
            replace_mapbox_tileset(
                self.mapbox[visualization]["polbnda_subn"]["mapid"],
                self.mapbox_auth,
                self.mapbox[visualization]["polbnda_subn"]["name"],
                polygon_to_upload,
                self.temp_folder,
            )
            replace_mapbox_tileset(
                self.mapbox[visualization]["polbndp_subn"]["mapid"],
                self.mapbox_auth,
                self.mapbox[visualization]["polbndp_subn"]["name"],
                point_to_upload,
                self.temp_folder,
            )
            # admin0 boundaries should include disputed areas and be dissolved to single features
            res = "1m"
            if len(self.countries[visualization]["adm0"]) > 25:
                res = "15m"
            to_upload = self.boundaries[f"polbnda_int_{res}"].copy(deep=True)
            to_upload = to_upload[to_upload["ISO_3"].isin(self.countries[visualization]["adm0"])]
            to_upload = to_upload.dissolve(by="ISO_3")
            to_upload["ISO_3"] = to_upload.index
            to_upload.reset_index(drop=True, inplace=True)
            to_upload = drop_fields(to_upload, ["ISO_3"])
            to_upload = merge(
                to_upload,
                self.boundaries[f"polbnda_int_{res}"][["ISO_3", "STATUS", "Color_Code", "Terr_ID", "Terr_Name"]],
                on="ISO_3",
            )
            replace_mapbox_tileset(
                self.mapbox[visualization]["polbnda_int"]["mapid"],
                self.mapbox_auth,
                self.mapbox[visualization]["polbnda_int"]["name"],
                to_upload,
                self.temp_folder,
            )
            replace_mapbox_tileset(
                self.mapbox[visualization]["polbndl_int"]["mapid"],
                self.mapbox_auth,
                self.mapbox[visualization]["polbndl_int"]["name"],
                self.boundaries[f"polbndl_int_{res}"][
                    (self.boundaries[f"polbndl_int_{res}"]["BDY_CNT01"].isin(self.countries[visualization]["adm0"]))
                    | (self.boundaries[f"polbndl_int_{res}"]["BDY_CNT02"].isin(self.countries[visualization]["adm0"]))
                ],
                self.temp_folder,
            )
            replace_mapbox_tileset(
                self.mapbox[visualization]["polbndp_int"]["mapid"],
                self.mapbox_auth,
                self.mapbox[visualization]["polbndp_int"]["name"],
                self.boundaries["polbndp_int_1m"][
                    self.boundaries["polbndp_int_1m"]["ISO_3"].isin(self.countries[visualization]["adm0"])
                ],
                self.temp_folder,
            )

    def update_lookups(self, visualizations, save=True):
        for visualization in visualizations:
            logger.info(f"Updating subnational lookups for {visualization}")
            subn_boundaries, _ = self.merge_subn_boundaries(visualization)
            attributes = dict()
            for _, row in subn_boundaries.iterrows():
                level = row["ADM_LEVEL"]
                new_name = (row["ADM_REF"].replace("-", " ").replace("`", "").replace("'", ""))
                new_name = (
                    normalize("NFKD", new_name)
                    .encode("ascii", "ignore")
                    .decode("ascii")
                )
                attribute = {
                    "country": row["ADM0_REF"],
                    "iso3": row["alpha_3"],
                    "pcode": row["ADM_PCODE"],
                    "name": new_name,
                }
                if level not in attributes:
                    attributes[level] = list()
                attributes[level].append(attribute)
            for level in attributes:
                attributes[level] = sorted(attributes[level], key=lambda i: (i["country"], i["name"]))

            self.lookups[visualization] = attributes
            if save:
                with open(join("saved_outputs", f"subnational-attributes-{visualization}.txt"), "w") as f:
                    for level in attributes:
                        f.write(f"adm{level}\n")
                        for row in attributes[level]:
                            if "," in row["name"]:
                                row["name"] = '"' + row["name"] + '"'
                            f.write("- %s\n" % str(row).replace("'", ""))

        logger.info("Updated subnational lookups")

    def update_bboxes(self, visualizations, HRPs, save=True):
        for visualization in visualizations:
            logger.info(f"Updating regional bbox jsons for {visualization}")
            regional_file = join("saved_outputs", f"ocha-regions-bbox-{visualization}.geojson")

            # select only territories or disputed areas in visualization
            adm0_region = self.boundaries["polbnda_int_1m"].copy(deep=True)
            adm0_region = adm0_region[adm0_region["ISO_3"].isin(self.countries[visualization]["adm0"])]
            # assign region and HRP status
            adm0_region["region"] = ""
            adm0_region["HRPs"] = ""
            adm0_region["All"] = "All"
            adm0_region.loc[adm0_region["ISO_3"].isin(HRPs), "HRPs"] = "HRPs"
            dataset = Dataset.read_from_hdx(self.regional_info["dataset"])
            resource = [
                r
                for r in dataset.get_resources()
                if r.get_file_type() == self.regional_info["format"]
            ]
            _, iterator = self.downloader.get_tabular_rows(
                resource[0]["url"], dict_form=True
            )
            for row in iterator:
                adm0_region.loc[
                    adm0_region["ISO_3"] == row[self.regional_info["iso3_header"]],
                    "region",
                ] = row[self.regional_info["region_header"]]
            # dissolve boundaries by region and HRP
            adm0_dissolve = adm0_region.dissolve(by="region")
            adm0_dissolve_HRPs = adm0_region[adm0_region["HRPs"] == "HRPs"].dissolve(
                by="HRPs"
            )
            adm0_dissolve_all = adm0_region.dissolve(by="All")
            adm0_dissolve = adm0_dissolve.append(adm0_dissolve_HRPs)
            adm0_dissolve = adm0_dissolve.append(adm0_dissolve_all)
            # convert to bounding box
            adm0_dissolve = adm0_dissolve.bounds
            adm0_dissolve["geometry"] = [
                box(l, b, r, t)
                for l, b, r, t in zip(
                    adm0_dissolve["minx"],
                    adm0_dissolve["miny"],
                    adm0_dissolve["maxx"],
                    adm0_dissolve["maxy"],
                )
            ]
            adm0_dissolve = GeoDataFrame(adm0_dissolve["geometry"])
            adm0_dissolve["region"] = adm0_dissolve.index
            adm0_region = adm0_dissolve.to_json(show_bbox=True, drop_id=True)
            adm0_region = loads(adm0_region)
            # remove any "NO_COVERAGE" regions
            for i in reversed(range(len(adm0_region["features"]))):
                if (
                    adm0_region["features"][i]["properties"][
                        "region"
                    ]
                    == "NO COVERAGE"
                ):
                    adm0_region["features"].remove(adm0_region["features"][i])
            adm0_region["name"] = "ocha regions - bbox"
            adm0_region["crs"] = {
                "type": "name",
                "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"},
            }
            self.bboxes[visualization] = adm0_region
            if save:
                try:
                    remove(regional_file)
                except FileNotFoundError:
                    logger.info("File does not exist - creating!")
                with open(regional_file, "w") as f_open:
                    dump(adm0_region, f_open)

        logger.info("Updated regional bbox jsons")
