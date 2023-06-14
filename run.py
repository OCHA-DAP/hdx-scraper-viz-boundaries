import argparse
import logging
import warnings
from os import getenv
from os.path import expanduser, join
from shapely.errors import ShapelyDeprecationWarning

from hdx.api.configuration import Configuration
from hdx.facades.keyword_arguments import facade
from hdx.utilities.downloader import Download
from hdx.utilities.easy_logging import setup_logging
from hdx.utilities.path import temp_dir
from boundaries import Boundaries

setup_logging()
logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

lookup = "hdx-scraper-viz-inputs"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-co", "--countries", default=None, help="Which countries to update")
    parser.add_argument("-vi", "--visualizations", default=None, help="Which visualizations to update")
    parser.add_argument("-ut", "--update_tilesets", default=False, help="Update mapbox tilesets")
    parser.add_argument("-ub", "--update_subn_bounds", default=False, help="Re-adjust country boundaries")
    parser.add_argument("-ma", "--mapbox_auth", default=None, help="Credentials for accessing MapBox data")
    args = parser.parse_args()
    return args


def main(
    mapbox_auth,
    countries=None,
    visualizations=None,
    update_tilesets=False,
    update_subn_bounds=False,
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

            boundaries = Boundaries(
                configuration,
                downloader,
                mapbox_auth,
                temp_folder,
            )

            boundaries.download_boundary_inputs(configuration["hdx_inputs"]["dataset"], levels)

            if update_subn_bounds:
                boundaries.update_subnational_boundaries(
                    countries_to_process,
                    configuration["boundaries"].get("do_not_process", []),
                )
                boundaries.update_subnational_resources(configuration["hdx_inputs"]["dataset"], levels)
            if update_tilesets:
                boundaries.update_mapbox_tilesets(visualizations)
            boundaries.update_lookups(visualizations)
            boundaries.update_bboxes(visualizations, configuration["HRPs"])

            logger.info("Finished processing!")


if __name__ == "__main__":
    args = parse_args()
    countries = args.countries
    if countries is None:
        countries = getenv("COUNTRIES")
    if countries:
        countries = countries.split(",")
    visualizations = args.visualizations
    if visualizations is None:
        visualizations = getenv("VISUALIZATIONS")
    if visualizations:
        visualizations = visualizations.split(",")
    update_tilesets = args.update_tilesets
    if not update_tilesets:
        update_tilesets = getenv("UPDATE_TILESETS")
    update_subn_bounds = args.update_subn_bounds
    if not update_subn_bounds:
        update_subn_bounds = getenv("UPDATE_SUBN_BOUNDS")
    mapbox_auth = args.mapbox_auth
    if mapbox_auth is None:
        mapbox_auth = getenv("MAPBOX_AUTH")
    facade(
        main,
        user_agent_config_yaml=join(expanduser("~"), ".useragents.yml"),
        user_agent_lookup=lookup,
        project_config_yaml=join("config", "project_configuration.yml"),
        countries=countries,
        visualizations=visualizations,
        update_tilesets=update_tilesets,
        update_subn_bounds=update_subn_bounds,
        mapbox_auth=mapbox_auth,
    )
