# Data
Sylvain Schmitt -
Apr 4, 2024

All data needed for the analyses and the scripts to retrieve them.

Scripts:

- get_data.R: this script takes advantage of Paracou Plot 6 in 2006
  (public data of the LoggingLab R package) to create dummy data as an
  example. To be run using `source("data/get_data.R").`

Data:

- dummy_data.tsv: dummy data as an example from Paracou Plot 6 in 2006
  (public data of the LoggingLab R package)

``` r
fs::dir_tree()
```

    .
    ├── README.md
    ├── README.qmd
    ├── README.rmarkdown
    ├── dummy_data.tsv
    └── get_data.R
