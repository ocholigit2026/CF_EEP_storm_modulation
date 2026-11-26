# CF_EEP_storm_modulation

MATLAB analysis code for *Storm-strength modulation of solar wind coupling
function-energetic electron precipitation relationship* (Ocholi et al., 2026).

## Structure

```
 functions/
  load_and_prepare_data.m      Loads and time-aligns solar wind, coupling
                                function, and POES/MetOp MEPED EEP flux data
  identify_storms_new.m        Storm identification, epoch extraction, and
                                weak/moderate/intense classification
  select_peak_events.m         Helper: enforces minimum separation between
                                candidate storm peaks
  extract_flux_timeseries_daily.m    Extract satellite data from defined L-shell
                                       range
  CFs_daily_mean.m              Calculate daily averages of coupling functions from
                                  hourly data.   
     
 data/
  generate_table1_regression_results.m   Table 1: storm-regime regression

 figures/
  generate_figure1_lagged_correlation.m       Figure 1
  generate_figure2_integration_window.m       Figure 2
  generate_figure3_native_binned_response.m      Figure 3
  generate_figure4_normalized_binned_response.m  Figure 4
  generate_figure1_daily.m      Figure S1
  generate_figure2_daily.m      Figure S2

```

## Usage

Each script in `data/` and `figures/` is self-contained: it calls
`load_and_prepare_data_new` and `identify_storms_new`, then performs its analysis.
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
