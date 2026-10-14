# Snow parameter estimation with a Large Sample Emulator (LSE)

Code used to assign SUMMA snow parameters to the model elements (HRUs) of the
three study basins in Casson et al. (EGUsphere preprint,
[doi:10.5194/egusphere-2025-6066](https://doi.org/10.5194/egusphere-2025-6066)),
following Sect. 4.5 and Appendix C of the manuscript.

Eleven SUMMA snow parameters are estimated. Rather than calibrating each basin
directly, a Large Sample Emulator is trained on snow-station simulations and
then used as an inexpensive surrogate to search the parameter space at each
HRU, conditional on that HRU's local climate and physiographic attributes.
The target metric is the Taylor Skill Score (TSS; Taylor, 2001).

| Basin    | HRUs |
|----------|------|
| Bow      | 1181 |
| Chena    |  329 |
| Tuolumne |  127 |

## Emulator

The emulator is trained outside this folder, following Tang et al. (2025), 
on 162 snow stations selected by K-means clustering of climate attributes 
from 770 candidate stations (NorSWE and SNOTEL), forced with ERA5. 
Each station contributes 1000 SUMMA runs, built up over an iterative sequence 
starting from 500 parameter sets (162 000 SUMMA runs in total). The trained 
random forests are read from `emulators/RF_step_<i>_<metric>.RData` as `emulatorSummary$emulator`; 
the emulator used here corresponds to step `iStep` of that sequence.

Predictors are the eleven snow parameters, normalised to [0, 1], together with
26 site attributes (elevation plus climate attributes). The emulator therefore
returns the expected TSS for a candidate parameter set at a given site.

## Steps

### 1. Site attributes — `Rfiles/01_ObtainAtts.R`

Computes the 26 attributes for every HRU of one basin from the GPEP ensemble
mean forcing (`summa_mean_forcing.nc`) and the SUMMA attributes file
(`attributes.nc`), over water years 2011-2018 (2010-10-01 to 2018-09-30).
Attributes cover precipitation, temperature, PET and aridity, precipitation
seasonality, snow fraction, wind, shortwave radiation and specific humidity,
each as an annual mean and as winter (DJF) and spring (MAM) means.

Set `basin`, `dirForcing`, `dirSetting` and `dirSave` at the top of the script,
and run one basin at a time. Attributes for each HRU are computed in parallel
with `mclapply`; `fn_getMeanForcing.R` holds the per-HRU calculation.

Output: `<basin>_avgForcing.RData` and the matching `|`-delimited `.txt`, plus
`CheckHruIds.RData` for confirming that the HRU order of the forcing and
attributes files agree.

### 2. Parameter search — `Rfiles/02_ObtainParams.R`

Loads the emulator and, for each HRU, searches the emulated response surface
with a genetic algorithm (`GA::ga`, real-valued, population `gaPop`, `gaIters`
generations, initialised from a Latin hypercube). Only one search is carried
out — no further SUMMA iterations — and the parameter set with the best
emulated TSS is assigned to the HRU. `fn_paramToMetric.R` evaluates a candidate
set through the emulator; parameter ranges and transformations are read from
`usefulFiles/snow_parameters_<paramVersion>.xlsx`, and the [0, 1] values
returned by the GA are back-transformed to physical units.

Run through the SLURM scripts in `Tasks/`, which call:

```bash
Rscript 02_ObtainParams.R ${dirYamlFile} ${yamlFileName} ${nCores} ${basin} ${targetMetric}
```

Output, per HRU, under `output/<basin>/hruID_<id>/`:

- `RF_BestParamSet_targetMetric_<metric>.txt` — the parameter set passed to SUMMA
- `RF_allParams_targetMetric<metric>.RData` — the best `nNewParams` candidate sets
- `RF_paramSearch_targetMetric_<metric>.RData` — the GA population and fitness log

An HRU is skipped if its `BestParamSet` file already exists, so interrupted
jobs can be resubmitted as-is.

### 3. Emulator error figure — `plots/01_calibration_error.R`

Produces the Appendix C figure: emulated against SUMMA TSS for the 162 training
stations, and the error across iteration steps for the 608 stations not used in
training.

## Configuration

`Rfiles/config_LSE.yaml` holds the settings shared by step 2:

| Key | Meaning |
|---|---|
| `paramVersion` | version of the snow parameter table in `usefulFiles/` |
| `versionAtts`  | version of the attribute table read from `dirAtts` |
| `iStep`        | emulator iteration step to use |
| `nNewParams`   | candidate parameter sets kept from the GA population |
| `gaPop`, `gaIters` | GA population size and number of generations |
| `dir*`         | locations of the emulators, functions, attributes, parameter table and output |

## Requirements

Run on Anvil with:

```bash
module load r/4.4.1 netcdf-c/4.7.4 netcdf-fortran/4.5.3
```

R packages: `ncdf4`, `stringr`, `lubridate`, `dplyr`, `readxl`, `yaml`,
`parallel`, `doParallel`, `lhs`, `GA`, `ranger`, `xgboost`, `airGR`, `SoilHyP`
(steps 1-2) and `ggplot2`, `sf`, `hexbin`, `patchwork`, `cowplot` (step 3).

## Contact
dave.casson@ucalgary.ca
nicolas.vasquez@ucalgary.ca
