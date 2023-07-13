### Updater for HDX data explorer boundaries

This repo updates the national and subnational boundaries used in HDX data explorers such as the [COVID viewer](https://data.humdata.org/visualization/covid19-humanitarian-operations/), [Ukraine Data Explorer](https://data.humdata.org/visualization/ukraine-humanitarian-operations/), and others.

### Usage

    python run.py

For the script to run, you will need to have a file called .hdx_configuration.yml in your home directory containing your HDX key eg.

    hdx_key: "XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX"
    hdx_read_only: false
    hdx_site: prod
    
You will also need to supply the universal .useragents.yml file in your home directory as specified in the parameter *user_agent_config_yaml* passed to facade in run.py. The collector reads the key **hdx-scraper-viz-inputs** as specified in the parameter *user_agent_lookup*.
 
Alternatively, you can set up environment variables: USER_AGENT, HDX_KEY, HDX_SITE.
 
Other needed environment variables are: PREPREFIX, VISUALIZATIONS, MAPBOX_AUTH.

### Process

National and subnational administrative boundaries are downloaded from a private dataset on HDX, subsetted to match the countries in each visualization, and Mapbox tilesets used in the data explorers are updated.

Then bounding box geojsons for OCHA regions are generated for each visualization along with subnational unit info text documents. The subnational info docs are used in the scrapers that generate the input data for the explorers.
