# Data
Sylvain Schmitt -
Jun 12, 2024

All data needed for the analyses and the scripts to retrieve them.

Scripts:

- get_data.R: this script takes advantage of Paracou Plot 6 in 2006
  (public data of the LoggingLab R package) to create dummy data as an
  example. To be run using `source("data/get_data.R").`

Data:

- dummy_data.tsv: dummy data as an example from Paracou Plot 6 in 2006
  (public data of the LoggingLab R package)
- JRC_TMF_AnnualChange_v1_2016_SAM_ID49_N10_W60.tif: forest cover
  evolution at 60W 10N in 2016 from the Tropical Moist Forest (TMF)
  product from the European Joint Research Centre (JRC) downloaded
  manually at <https://forobs.jrc.ec.europa.eu/TMF/data#downloads>
  (should be changes my automated script)
- French-Guiana_chelsa2_monthly-means_1980-2005.nc: was downloaded using
  the DownClim workflow: <https://github.com/sylvainschmitt/DownClim>
  (should be changes my automated script)

``` r
fs::dir_tree()
```

    .
    ├── French-Guiana_chelsa2_monthly-means_1980-2005.nc
    ├── JRC_TMF_AnnualChange_v1_2016_SAM_ID49_N10_W60.tif
    ├── README.md
    ├── README.qmd
    ├── README.rmarkdown
    └── harmonized_datasets_ss
        ├── output_MISIONES.csv
        ├── output_data-paracou-t2.csv
        ├── output_data-paracou-t3.csv
        ├── output_data-temoin.csv
        ├── output_data_paracou-t1.csv
        └── output_moju.csv
