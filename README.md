# diversity
Oct 16, 2024

[![](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![lint](https://github.com/Bioforest-project/diversity/workflows/lint/badge.svg)](https://github.com/Bioforest-project/diversity/actions?query=workflow%3Alint)

**diversity** is a sub-project of the **BioForest** project aimed at
studying the Interactions between tree biodiversity, forest dynamics and
climate in managed tropical forests with a pan-tropical approach.
Specifically, **diversity** focus on the **effects of logging
disturbance on tropical forest diversity - a global assessment** with
the following question: What is the effect of logging intensity on
biodiversity, and how does it change through time?

## Usage

All **diversity** analyses rely on the quarto documents (`files.qmd`)
that can be run with R and associated environment defined with
[renv](https://github.com/Bioforest-project/species#0).

## Project

**species** includes:

- Analyse of the data with associated [documentation and
  figures](https://bioforest-project.github.io/LoggingDiversity/):
  - Reproductive analyses in `files.qmd`
  - Resulting pages in `docs/`
  - Document structure definition in `_quarto.yml`
- Data in `data/` with:
  - All raw data in `raw_data/`

  - All derived data in `derived_sata/`
- Intermediary files in `outputs/`
- Figures in `figures/`
- Models in `models/`
- R environment definition with
  [renv](https://rstudio.github.io/renv/articles/renv.html) in `renv/`
  and `renv/lock`
- R files (`.Rbuildignore` , `.Rdata` , `.Rprofile` , `.Rhistory`)
- Git and GitHub files (`.gitignore` , `.github/`)
- Project documentation (`README.qmd` , `README.md` , `NEWS.md`,
  `LICENSE` )

## Contribution

You can contribute to the project by forking the repository on github
and cloning the fork to your machine using several options, including
GitHub desktop GUI. Further informations on contribution are detailed in
the online document:
<https://bioforest-project.github.io/data_preparation/98_contributing.html>.

## Help

Please preferentially create an issue on GitHub for any questions, bugs
or help needed regarding **diversity**:
<https://github.com/Bioforest-project/LoggingDiversity/issues> . You may
however reach us by mail with people from the core group (see below).

## Core group

- Sylvain Schmitt (sylvain.schmitt@cirad.fr)
- Mithila Unkule (mithila.unkule@fondationbiodiversite.fr)
- Genoveva Gatti (genogatti@gmail.com)
- David Burslem (d.burslem@abdn.ac.uk)
- Andes Hamuraby Rozak (andes.hamuraby.rozak@brin.go.id)
- Natalia Bedrij (nabedrij@gmail.com)
- Verginia Wortel (wortelv@gmail.com)

The whole group consist of participants to the [Bioforest
project](https://www.fondationbiodiversite.fr/la-frb-en-action/programmes-et-projets/le-cesab/bioforest/).

![](https://www.fondationbiodiversite.fr/wp-content/uploads/2023/10/bioforest-ws1_web.jpeg)
