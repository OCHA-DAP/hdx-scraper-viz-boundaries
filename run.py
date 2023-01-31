import argparse
import logging
import warnings
from geopandas import read_file
from os import getenv
from os.path import expanduser, join
from shapely.errors import ShapelyDeprecationWarning

from hdx.api.configuration import Configuration
from hdx.data.dataset import Dataset
from hdx.facades.keyword_arguments import facade
from hdx.utilities.downloader import Download
from hdx.utilities.easy_logging import setup_logging
from hdx.utilities.path import temp_dir
from boundaries import Boundaries

setup_logging()
logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

lookup = "hdx-scraper-viz-boundaries"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-hk", "--hdx_key", default=None, help="HDX api key")
    parser.add_argument("-ua", "--user_agent", default=None, help="user agent")
    parser.add_argument("-pp", "--preprefix", default=None, help="preprefix")
    parser.add_argument("-hs", "--hdx_site", default=None, help="HDX site to use")
    parser.add_argument("-co", "--countries", default=None, help="Which countries to update")
    parser.add_argument("-vi", "--visualizations", default=None, help="Which visualizations to update")
    parser.add_argument("-ut", "--update_tilesets", default=None, help="Update mapbox tilesets")
    parser.add_argument("-ub", "--update_subn_bounds", default=None, help="Re-adjust country boundaries")
    parser.add_argument("-ma", "--mapbox_auth", default=None, help="Credentials for accessing MapBox data")
    args = parser.parse_args()
    return args


def main(
    countries,
    visualizations,
    mapbox_auth,
    update_tilesets=None,
    update_subn_bounds=None,
    **ignore,
):
    logger.info(f"##### hdx-viz-data-inputs ####")
    configuration = Configuration.read()
    with temp_dir(folder="TempDataExplorerInputs") as temp_folder:
        with Download(rate_limit={"calls": 1, "period": 0.1}) as downloader:

            if not visualizations:
                visualizations = [key for key in configuration["countries"]]

            levels = [l for viz in visualizations for l in configuration["countries"][viz] if l != "adm0"]
            countries_to_process = {level: set() for level in levels}
            for viz in visualizations:
                for level in configuration["countries"][viz]:
                    if level == "adm0":
                        continue
                    for country in configuration["countries"][viz][level]:
                        if not countries or country in countries:
                            countries_to_process[level].add(country)
            for level in levels:
                countries_to_process[level] = list(countries_to_process[level])

            # download boundaries
            logger.info("Downloading boundaries")
            all_boundaries = dict()
            dataset = Dataset.read_from_hdx(configuration["hdx_inputs"]["dataset"])
            for resource in dataset.get_resources():
                if "coastl" in resource["name"]:
                    continue
                if resource["name"][8:12] not in countries_to_process \
                        and ("polbnda_adm" in resource["name"] or "polbndp_adm" in resource["name"]):
                    continue
                _, resource_file = resource.download(folder=temp_folder)
                lyr = read_file(resource_file)
                if "polbnda_adm" in resource["name"] or "polbndp_adm" in resource["name"]:
                    all_boundaries[resource["name"].replace("_1m_ocha.geojson", "")] = lyr
                if "lake" in resource["name"]:
                    all_boundaries["lake"] = lyr
                if "_int_" in resource["name"]:
                    all_boundaries[resource["name"].replace("wrl_", "").replace("_uncs.geojson", "")] = lyr

            bounds = Boundaries(
                configuration,
                downloader,
                all_boundaries,
                mapbox_auth,
                temp_folder,
            )
            if update_subn_bounds:
                bounds.update_subnational_boundaries(
                    countries_to_process,
                    configuration["boundaries"].get("do_not_process", []),
                )
                bounds.update_subnational_resources(dataset, levels)
            if update_tilesets:
                bounds.update_mapbox_tilesets(visualizations)
            bounds.update_lookups(visualizations)
            bounds.update_bboxes(visualizations, configuration["HRPs"])

            logger.info("Finished processing!")


if __name__ == "__main__":
    args = parse_args()
    hdx_key = args.hdx_key
    if hdx_key is None:
        hdx_key = getenv("HDX_KEY")
    user_agent = args.user_agent
    if user_agent is None:
        user_agent = getenv("USER_AGENT")
    preprefix = args.preprefix
    if preprefix is None:
        preprefix = getenv("PREPREFIX")
    hdx_site = args.hdx_site
    if hdx_site is None:
        hdx_site = getenv("HDX_SITE", "prod")
    countries = args.countries
    if countries is None:
        countries = getenv("COUNTRIES", None)
    if countries:
        countries = countries.split(",")
    visualizations = args.visualizations
    if visualizations is None:
        visualizations = getenv("VISUALIZATIONS", "all")
    if visualizations:
        visualizations = visualizations.split(",")
    update_tilesets = args.update_tilesets
    if not update_tilesets:
        update_tilesets = getenv("UPDATE_TILESETS", None)
    update_subn_bounds = args.update_subn_bounds
    if not update_subn_bounds:
        update_subn_bounds = getenv("UPDATE_SUBN_BOUNDS", None)
    mapbox_auth = args.mapbox_auth
    if mapbox_auth is None:
        mapbox_auth = getenv("MAPBOX_AUTH", None)
    facade(
        main,
        hdx_key=hdx_key,
        user_agent=user_agent,
        user_agent_config_yaml=join(expanduser("~"), ".useragents.yml"),
        user_agent_lookup=lookup,
        preprefix=preprefix,
        hdx_site=hdx_site,
        project_config_yaml=join("config", "project_configuration.yml"),
        countries=countries,
        visualizations=visualizations,
        update_tilesets=update_tilesets,
        update_subn_bounds=update_subn_bounds,
        mapbox_auth=mapbox_auth,
    )
