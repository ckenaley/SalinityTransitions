
## Supplementary Content for O’Connor and Kenaley

### Description of the data and file structure

The scripts and data files within this repository were used in the
analysis of O’Connor and Kenaley (202x). We undertook stochastic
character mapping and multistate speciation and extinction modeling of
salinity habitat and diadromy. This repository contains the scripts to
perform the analysis and (if on Dryad) the data that were produced.

Scripts include:

- `Salinity_Spatial.R`: GBIF querying and assignment of records to
  marine, fresh, or estuarine locations
- `Salinity_ML_analysis.R`: Random forest and gradient boosted machine
  learning analysis of halotolerance
- `Salinity_SIMMAP.R`: Stochastic character mapping
- `Salinity_MuSSE.R`: Multi-State speciation and extinction modeling
- `Salinity_Figures.R`: Code that produced (rough drafts of) figures in
  the paper

Table S1 is rendered by the RMD `Table_S1.RMD`

### Running the analysis

For our analysis, we executed the scripts in the order above. They
require a number of data files, including shape files for spatial
analysis. These data can be downloaded from
[Dryad](https://doi.org/10.5061/dryad.bzkh189g4) (the “data.zip”
directory in that repository). We suggest downloading the data,
uncompressing them, and placing them in the same directory as the
scripts.

This data directory contains output from the analysis as well. One can,
of course, ignore these and just run the analysis; however, this is time
consuming, on the order of many hours. Therefore, we make the outputs of
each script available so that one can run each independently.

### Access information

These data and code are available on
[Dryad](https://doi.org/10.5061/dryad.bzkh189g4) and
[github](https://github.com/ckenaley/SalinityTransitions)

### Contact

For questions about the analysis, contact the corresponding author at
kenaley at bc.edu.
