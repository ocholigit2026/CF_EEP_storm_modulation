function D = load_and_prepare_data_new()
% LOAD_AND_PREPARE_DATA
% Load solar wind, coupling functions and EEP flux data.
% Produce both hourly and daily time-aligned datasets.

%% Explicit paths to data files

fluxFile = 'C:\Users\mamud\DATA HOUSE\poes\raw_data\MEPED_particle_reproduced_optimized.mat';

swFile = 'C:\Users\mamud\DATA HOUSE\OMNIWEB DATA\hourly\swdataTT_2009_2019H_filtered.mat';

coupleFile = 'C:\Users\mamud\DATA HOUSE\energy_couple_fns\energy_couple_data_2009_2019.mat';


%% Load variables

load(fluxFile, 'fluxdata_2009_2019H_4D_new');
load(swFile, 'swdata_2009_2019H');
load(coupleFile, 'coupling_fns_tt');


%% Synchronize solar wind and coupling-function timetables

cpl_swdata = synchronize(swdata_2009_2019H, coupling_fns_tt);

time = cpl_swdata.Time;


%% Extract flux channels (0-degree telescope)

mep0e1_3D = squeeze(fluxdata_2009_2019H_4D_new(7,:,:,:));
mep0e2_3D = squeeze(fluxdata_2009_2019H_4D_new(8,:,:,:));
mep0e3_3D = squeeze(fluxdata_2009_2019H_4D_new(9,:,:,:));


%% Generate hourly and daily EEP indices

[D.flux30_hourly, D.flux30_daily, D.daily_time, D.n30] = ...
    extract_flux_timeseries_daily(mep0e1_3D, time);

[D.flux100_hourly, D.flux100_daily, ~, D.n100] = ...
    extract_flux_timeseries_daily(mep0e2_3D, time);

[D.flux300_hourly, D.flux300_daily, ~, D.n300] = ...
    extract_flux_timeseries_daily(mep0e3_3D, time);


%% Extract hourly solar wind and coupling functions

D.time   = time;
D.Ey     = cpl_swdata.vBs;
D.Newell = cpl_swdata.dphidt;
D.KL     = cpl_swdata.Ekl;
D.EPS    = cpl_swdata.eps;
D.Dst    = cpl_swdata.Dst;


%% Calculate daily coupling-function and Dst averages

[D.Ey_daily, ~, D.nEy] = CFs_daily_mean(D.Ey, time);

[D.Newell_daily, ~, D.nNewell] = CFs_daily_mean(D.Newell, time);

[D.KL_daily, ~, D.nKL] = CFs_daily_mean(D.KL, time);

[D.EPS_daily, ~, D.nEPS] = CFs_daily_mean(D.EPS, time);

[D.Dst_daily, ~, D.nDst] = CFs_daily_mean(D.Dst, time);


%% Calculate daily EEP availability is already provided by flux function

% D.daily_time is the common daily timestamp vector generated from time.


%% Remove rows with missing values from HOURLY dataset only

valid = ~isnan(D.flux30_hourly) & ...
        ~isnan(D.flux100_hourly) & ...
        ~isnan(D.flux300_hourly) & ...
        ~isnan(D.KL) & ...
        ~isnan(D.EPS) & ...
        ~isnan(D.Newell) & ...
        ~isnan(D.Ey);

hourly_fields = {'time','Ey','Newell','KL','EPS','Dst', ...
                 'flux30_hourly','flux100_hourly','flux300_hourly'};

for i = 1:numel(hourly_fields)

    field = hourly_fields{i};

    if strcmp(field,'time')
        D.time = D.time(valid);
    else
        D.(field) = D.(field)(valid);
    end

end

end
