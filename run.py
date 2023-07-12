import argparse
import logging
from os import getenv
from os.path import expanduser, join

from hdx.api.configuration import Configuration
from hdx.facades.keyword_arguments import facade
from hdx.utilities.downloader import Download
from hdx.utilities.easy_logging import setup_logging
from hdx.utilities.path import temp_dir
from boundaries import Boundaries

setup_logging()
logger = logging.getLogger(__name__)

lookup = "hdx-scraper-viz-inputs"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-vi", "--visualizations", default=None, help="Which visualizations to update")
    parser.add_argument("-ma", "--mapbox_auth", default=None, help="Credentials for accessing MapBox data")
    args = parser.parse_args()
    return args


def main(
    mapbox_auth,
    visualizations=None,
    **ignore,
):
    logger.info(f"##### hdx-viz-data-inputs ####")
    configuration = Configuration.read()
    with temp_dir(folder="TempDataExplorerInputs") as temp_folder:
        with Download(rate_limit={"calls": 1, "period": 0.1}) as downloader:

            if not visualizations:
                visualizations = [key for key in configuration["countries"]]

            levels = [l for viz in visualizations for l in configuration["countries"][viz] if l != "adm0"]

            boundaries = Boundaries(
                configuration,
                downloader,
                mapbox_auth,
                temp_folder,
            )

            boundaries.download_boundary_inputs(configuration["hdx_dataset"], levels)
            boundaries.update_mapbox_tilesets(visualizations, configuration["mapbox"])
            boundaries.update_lookups(visualizations)
            boundaries.update_bboxes(visualizations, configuration["HRPs"])

            logger.info("Finished processing!")


if __name__ == "__main__":
    args = parse_args()
    visualizations = args.visualizations
    if visualizations is None:
        visualizations = getenv("VISUALIZATIONS")
    if visualizations:
        visualizations = visualizations.split(",")
    mapbox_auth = args.mapbox_auth
    if mapbox_auth is None:
        mapbox_auth = getenv("MAPBOX_AUTH")
    facade(
        main,
        user_agent_config_yaml=join(expanduser("~"), ".useragents.yml"),
        user_agent_lookup=lookup,
        project_config_yaml=join("config", "project_configuration.yml"),
        visualizations=visualizations,
        mapbox_auth=mapbox_auth,
    )
