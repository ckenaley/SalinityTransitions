
## Supplementary Content for O’Connor and Kenaley

### Description of the data and file structure

The scripts and data files within this repository were used in the
analysis of O’Connor and Kenaley (202x). We undertook stochastic
character mapping and multistate speciation and extinction modeling of
salinity habitat and diadromy. This repository contains the scripts to
perform the analysis and the data that were produced.

Scripts include:

- `Salinity_Spatial.R`: GBIF querying and assignment of records to
  marine, fresh, or estuarine locations
- `Salinity_ML_analysis.R`: Random forest and gradient boosted machine
  learning analysis of halotolerance
- `Salinity_SIMMAP.R`: Stochastic character mapping
- `Salinity_MuSSE.R`: Multi-State speciation and extinction modeling
- `Salinity_Figures.R`: Code that produced (rough drafts of) figures in
  the paper

Data produced by these scripts and geospatial data are placed withing
the `data` directory.

Table S1 is rendered by the RMD `Table_S1.RMD`

### Running the analysis

Each script loads required data saved in RDS (all contained in the
`data` directory), so one can begin anywhere. For our analysis, we
executed the scripts in the order above.

### Access information

These data and code are available on [Dryad](10.5061/dryad.bzkh189g4)
and [github](https://github.com/ckenaley/SalinityTransitions)
