import logging
from geojson import loads
from geopandas import GeoDataFrame, read_file
from json import dump
from mapbox import Uploader
from os import remove
from os.path import join
from pandas import concat, merge
from shapely.geometry import box
from time import sleep
from unicodedata import normalize

from hdx.data.dataset import Dataset

logger = logging.getLogger()


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
        self, configuration, downloader, mapbox_auth, temp_folder
    ):
        self.downloader = downloader
        self.boundaries = dict()
        self.temp_folder = temp_folder
        self.mapbox = configuration["mapbox"]
        self.mapbox_auth = mapbox_auth
        self.countries = configuration["countries"]
        self.regional_info = configuration["regional"]
        self.bboxes = dict()
        self.lookups = dict()

    def download_boundary_inputs(self, dataset_name, levels):
        logger.info("Downloading boundaries")
        all_boundaries = dict()
        dataset = Dataset.read_from_hdx(dataset_name)
        for resource in dataset.get_resources():
            if "coastl" in resource["name"] or "lake" in resource["name"]:
                continue
            if "adm" in resource["name"] and resource["name"].split("_")[1] not in levels:
                continue
            _, resource_file = resource.download(folder=self.temp_folder)
            lyr = read_file(resource_file)
            if "polbnda_adm" in resource["name"] or "polbndp_adm" in resource["name"]:
                all_boundaries[resource["name"].replace("_1m_ocha.geojson", "")] = lyr
            if "_int_" in resource["name"]:
                all_boundaries[resource["name"].replace("wrl_", "").replace("_uncs.geojson", "")] = lyr

        self.boundaries = all_boundaries

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
            to_upload = to_upload.drop(
                [f for f in to_upload.columns if f != "ISO_3" and f.lower() != "geometry"],
                axis=1,
            )
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
            adm0_dissolve = concat([adm0_dissolve, adm0_dissolve_HRPs, adm0_dissolve_all])
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
