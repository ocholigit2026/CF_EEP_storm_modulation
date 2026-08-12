# CF_EEP_storm_modulation

MATLAB analysis code for *Storm-strength modulation of solar wind coupling
function-energetic electron precipitation relationship* (Ocholi et al., 2026).

## Structure

```
data/
  load_and_prepare_data.m      Loads and time-aligns solar wind, coupling
                                function, and POES/MetOp MEPED EEP flux data

storms/
  identify_storms.m            Storm identification, epoch extraction, and
                                weak/moderate/intense classification
  select_peak_events.m         Helper: enforces minimum separation between
                                candidate storm peaks

analysis/
  generate_table1_storm_distribution.m   Table 1: storm counts by regime
  generate_table2_regression_results.m   Table 2: storm-regime regression

figures/
  generate_figure1_lagged_correlation.m       Figure 1
  generate_figure2_integration_window.m       Figure 2
  generate_figure3_normalized_binned_response.m  Figure 3
  generate_figure4_native_binned_response.m      Figure 4
```

## Usage

Each script in `analysis/` and `figures/` is self-contained: it calls
`load_and_prepare_data` and `identify_storms`, then performs its analysis.
Run any script directly from within its folder (each adds `../data` and
`../storms` to the path).

## Required input files

Place these in the working directory (or update the path passed to
`load_and_prepare_data`):

- `fluxdata_2009_2019H.mat` - POES/MetOp MEPED particle flux (Evans & Greer, 2004)
- `swdataTT_2009_2019H_filtered.mat` - OMNIWeb solar wind and Dst data
- `energy_couple_data_2009_2019.mat` - derived coupling function time series

## Storm identification parameters

- Dst threshold: -30 nT
- Minimum storm separation: 96 hours
- Epoch window: 48 hours before to 72 hours after storm center
- Regime classification (by minimum Dst): weak (-30 to -50), moderate
  (-50 to -100), intense (\<= -100)
