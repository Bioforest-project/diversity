# Data
Sylvain Schmitt -
Jul 31, 2024

The workflow to get the data and the data needed for the analyses.

Workflow:

- Snakefile: main snakemake workflow file
- config/: all configuration files
- scripts/: all scripts files

Data:

- climate: climate data
- landscape: landscape data
- soil: soil data
- harmonized_dataset_ss: first harmonized shared inventories by Mithila

``` r
fs::dir_tree()
```

    .
    ├── README.md
    ├── README.qmd
    ├── README.rmarkdown
    ├── Snakefile
    ├── climate
    │   ├── missiones_climate.tsv
    │   └── paracou_climate.tsv
    ├── config
    │   ├── config.yml
    │   ├── logging_diversity.yml
    │   └── sites.tsv
    ├── harmonized_datasets_ss
    │   ├── output_MISIONES.csv
    │   ├── output_data-paracou-t2.csv
    │   ├── output_data-paracou-t3.csv
    │   ├── output_data-paracou-temoin.csv
    │   ├── output_data_paracou-t1.csv
    │   └── output_moju.csv
    ├── landscape
    │   ├── missiones_landscape.tsv
    │   └── paracou_landscape.tsv
    ├── scripts
    │   ├── get_climate.py
    │   ├── get_landscape.py
    │   └── get_soil.py
    └── soil
        ├── missiones_soil.tsv
        └── paracou_soil.tsv
