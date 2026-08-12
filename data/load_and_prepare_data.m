function D = load_and_prepare_data(dataDir)
% LOAD_AND_PREPARE_DATA Load solar wind, coupling function, and EEP flux
% data and merge into a single time-aligned dataset.
%
%   D = LOAD_AND_PREPARE_DATA(dataDir) loads the particle flux, solar
%   wind, and coupling function .mat files from dataDir (default:
%   current folder) and returns a struct D of time-aligned vectors,
%   with rows containing any missing coupling function or flux value
%   removed.
%
% Required files in dataDir:
%   fluxdata_2009_2019H.mat            (POES/MetOp MEPED particle flux)
%   swdataTT_2009_2019H_filtered.mat   (OMNIWeb solar wind / Dst data)
%   energy_couple_data_2009_2019.mat   (derived coupling function time series)
%
% Fields of D:
%   time                      - datetime vector
%   Ey, Newell, KL, EPS       - coupling function time series
%   Dst                       - Dst geomagnetic index
%   flux30, flux100, flux300  - EEP flux (0-deg telescope, MLT- and
%                                L=4-6-averaged), for >30/>100/>300 keV

if nargin < 1
    dataDir = pwd;
end

load(fullfile(dataDir, 'fluxdata_2009_2019H.mat'));
load(fullfile(dataDir, 'swdataTT_2009_2019H_filtered.mat'));
load(fullfile(dataDir, 'energy_couple_data_2009_2019.mat'));

cpl_swdata = synchronize(swdata_2009_2019H, coupling_fns_tt);

% Particle data is stored as [variable x time x MLT x L]. Variable order:
% lval(1), folon(2), folat(3), mlt(4), pas0(5), pas90(6),
% mep0e1(7), mep0e2(8), mep0e3(9), mep90e1(10), mep90e2(11), mep90e3(12)
% Only the 0-deg telescope channels (7-9) are used in this analysis.
mep0e1_3D = squeeze(fluxdata_2009_2019H_4D(7, :, :, :));
mep0e2_3D = squeeze(fluxdata_2009_2019H_4D(8, :, :, :));
mep0e3_3D = squeeze(fluxdata_2009_2019H_4D(9, :, :, :));

cpl_swdata.flux30  = extract_flux_timeseries(mep0e1_3D);
cpl_swdata.flux100 = extract_flux_timeseries(mep0e2_3D);
cpl_swdata.flux300 = extract_flux_timeseries(mep0e3_3D);

D.time    = cpl_swdata.Time;
D.Ey      = cpl_swdata.vBs;
D.Newell  = cpl_swdata.dphidt;
D.KL      = cpl_swdata.Ekl;
D.EPS     = cpl_swdata.eps;
D.Dst     = cpl_swdata.Dst;
D.flux30  = cpl_swdata.flux30;
D.flux100 = cpl_swdata.flux100;
D.flux300 = cpl_swdata.flux300;

% Remove records with a missing coupling function or flux value
valid = ~isnan(D.flux30) & ~isnan(D.flux100) & ~isnan(D.flux300) & ...
        ~isnan(D.KL) & ~isnan(D.EPS) & ~isnan(D.Newell) & ~isnan(D.Ey);

fn = fieldnames(D);
for i = 1:numel(fn)
    D.(fn{i}) = D.(fn{i})(valid);
end

end

function flux_ts = extract_flux_timeseries(flux_3D)
% Average a [time x MLT x L] flux array over all MLT sectors, then over
% the L = 4-6 shell range (array indices 13:21 in the source grid).
flux_ts = squeeze(mean(permute(flux_3D, [2 3 1]), 'omitmissing'));
flux_ts = squeeze(mean(flux_ts(13:21, :), 'omitmissing'))';
end
