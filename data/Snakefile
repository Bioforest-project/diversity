configfile: "config/config.yml"

import pandas as pd
sites = pd.read_table(config["sites"])

rule all:
   input:
      expand("climate/{site}_climate.tsv",
              site=sites.site),
      expand("soil/{site}_soil.tsv",
              site=sites.site),
      expand("landscape/{site}_landscape.tsv",
              site=sites.site)

rule get_climate:
    input:
        config["sites"]
    output:
        "climate/{site}_climate.tsv"
    log:
        "logs/{site}_climate.log"
    benchmark:
        "benchmarks/{site}_climate.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "config/logging_diversity.yml"
    params:
        site="{site}"
    script:
      "scripts/get_climate.py"
      
rule get_soil:
    input:
        config["sites"]
    output:
        "soil/{site}_soil.tsv"
    log:
        "logs/{site}_soil.log"
    benchmark:
        "benchmarks/{site}_soil.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "config/logging_diversity.yml"
    params:
        site="{site}"
    script:
      "scripts/get_climate.py"
      
rule get_landscape:
    input:
        config["sites"]
    output:
        "landscape/{site}_landscape.tsv"
    log:
        "logs/{site}_landscape.log"
    benchmark:
        "benchmarks/{site}_landscape.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "config/logging_diversity.yml"
    params:
        site="{site}",
        radius=config["landscape_radius"]
    script:
      "scripts/get_landscape.py"
      
