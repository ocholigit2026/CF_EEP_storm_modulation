function S = identify_storms(D, varargin)
% IDENTIFY_STORMS
% Identify storms from hourly Dst and extract hourly and daily
% storm-centred EEP flux and coupling-function data.
%
% Daily aggregation:
% Day -1: -36 to -13 hours
% Day  0: -12 to +11 hours
% Day +1: +12 to +35 hours
% Day +2: +36 to +59 hours
% Day +3: +60 to +83 hours
%
% Each daily interval contains 24 hourly observations.
% Minimum valid hourly observations per daily mean = 12.
%
% Original hourly outputs are preserved.

p = inputParser;

addParameter(p, 'Threshold', -30);
addParameter(p, 'MinSeparation', 96);
addParameter(p, 'PreWindow', 48);
addParameter(p, 'PostWindow', 72);
addParameter(p, 'MinDailyHours', 12);

parse(p, varargin{:});
opt = p.Results;

%% Extract hourly variables

Dst = D.Dst(:);
time = D.time(:);

flux30  = get_flux_field(D, 'flux30');
flux100 = get_flux_field(D, 'flux100');
flux300 = get_flux_field(D, 'flux300');

EPS    = D.EPS(:);
KL     = D.KL(:);
Ey     = D.Ey(:);
Newell = D.Newell(:);

%% Check time alignment

n = length(Dst);

if any([length(time), length(flux30), length(flux100), ...
        length(flux300), length(EPS), length(KL), ...
        length(Ey), length(Newell)] ~= n)

    error('The hourly data arrays and time vector must have equal lengths.');
end

%% Identify storm centres

[~, locs] = findpeaks(-Dst);

stormCandidates = locs(Dst(locs) <= opt.Threshold);

stormCenters = select_peak_events( ...
    stormCandidates, -Dst, opt.MinSeparation);

pre  = opt.PreWindow;
post = opt.PostWindow;

T = pre + post + 1;
N = length(stormCenters);

%% Allocate hourly storm arrays

Dst_st     = NaN(N,T);

flux30_st  = NaN(N,T);
flux100_st = NaN(N,T);
flux300_st = NaN(N,T);

eps_st     = NaN(N,T);
KL_st      = NaN(N,T);
Ey_st      = NaN(N,T);
Newell_st  = NaN(N,T);

stormSize = NaN(N,1);

%% Daily interval definitions

% Indices relative to the 121-hour storm window.
% The storm minimum is at index 49.

dailyBins = [
     13  36
    37  60
    61  84
    73  96
    85 108
];

nDaily = size(dailyBins,1);

DailyNames = {'Day -1','Day 0', ...
              'Day +1','Day +2', 'Day +3'};

% Daily arrays for individual storm events.
flux30_daily  = NaN(N,nDaily);
flux100_daily = NaN(N,nDaily);
flux300_daily = NaN(N,nDaily);

EPS_daily    = NaN(N,nDaily);
KL_daily     = NaN(N,nDaily);
Ey_daily     = NaN(N,nDaily);
Newell_daily = NaN(N,nDaily);

%% Extract hourly storm intervals and calculate daily means

for i = 1:N

    center = stormCenters(i);

    start_idx = center - pre;
    end_idx   = center + post;

    % Exclude storms too close to the data boundaries.
    if start_idx < 1 || end_idx > n
        continue
    end

    % Hourly storm intervals
    Dst_st(i,:)     = Dst(start_idx:end_idx)';

    flux30_st(i,:)  = flux30(start_idx:end_idx)';
    flux100_st(i,:) = flux100(start_idx:end_idx)';
    flux300_st(i,:) = flux300(start_idx:end_idx)';

    eps_st(i,:)     = EPS(start_idx:end_idx)';
    KL_st(i,:)      = KL(start_idx:end_idx)';
    Ey_st(i,:)      = Ey(start_idx:end_idx)';
    Newell_st(i,:)  = Newell(start_idx:end_idx)';

    % Classify storm using minimum Dst in the extracted interval.
    Dst_min = min(Dst_st(i,:), [], 'omitnan');

    if Dst_min <= -100
        stormSize(i) = 3;

    elseif Dst_min <= -50 && Dst_min > -100
        stormSize(i) = 2;

    elseif Dst_min <= opt.Threshold && Dst_min > -50
        stormSize(i) = 1;
    end

    % Calculate daily means for each complete 24-hour interval.
    for d = 1:nDaily

        ii = dailyBins(d,1):dailyBins(d,2);

        flux30_daily(i,d) = daily_mean( ...
            flux30_st(i,ii), opt.MinDailyHours);

        flux100_daily(i,d) = daily_mean( ...
            flux100_st(i,ii), opt.MinDailyHours);

        flux300_daily(i,d) = daily_mean( ...
            flux300_st(i,ii), opt.MinDailyHours);

        EPS_daily(i,d) = daily_mean( ...
            eps_st(i,ii), opt.MinDailyHours);

        KL_daily(i,d) = daily_mean( ...
            KL_st(i,ii), opt.MinDailyHours);

        Ey_daily(i,d) = daily_mean( ...
            Ey_st(i,ii), opt.MinDailyHours);

        Newell_daily(i,d) = daily_mean( ...
            Newell_st(i,ii), opt.MinDailyHours);

    end
end

%% Group storms by intensity

idx = {
    stormSize == 1
    stormSize == 2
    stormSize == 3
};

S.stormTimes  = time(stormCenters);
S.stormSize   = stormSize;

S.StormNames = {'Weak','Moderate','Intense'};

S.EnergyNames = {'>30 keV','>100 keV','>300 keV'};

S.CFNames = {'E_y','\epsilon','Kan-Lee','Newell'};

%% Preserve original hourly outputs

fluxArrays = {flux30_st, flux100_st, flux300_st};

S.Flux = cell(3,3);

for e = 1:3
    for r = 1:3

        S.Flux{e,r} = reshape( ...
            fluxArrays{e}(idx{r},:), [], 1);

    end
end

cfArrays = {Ey_st, eps_st, KL_st, Newell_st};

S.CF = cell(4,3);

for c = 1:4
    for r = 1:3

        S.CF{c,r} = reshape( ...
            cfArrays{c}(idx{r},:), [], 1);

    end
end

%% Preserve individual storm intervals by regime

fluxArrays = {flux30_st, flux100_st, flux300_st};
cfArrays = {Ey_st, eps_st, KL_st, Newell_st};

S.FluxEvents = cell(3,3);
S.CFEvents = cell(4,3);

S.EventDates = cell(1,3);
S.EventIndices = cell(1,3);

for r = 1:3

    eventIdx = find(idx{r});

    S.EventIndices{r} = eventIdx;
    S.EventDates{r} = S.stormTimes(eventIdx);

    for e = 1:3

        S.FluxEvents{e,r} = fluxArrays{e}(eventIdx,:);

    end

    for c = 1:4

        S.CFEvents{c,r} = cfArrays{c}(eventIdx,:);

    end

end


%% Store daily storm-event arrays

fluxDailyArrays = {
    flux30_daily
    flux100_daily
    flux300_daily
};

cfDailyArrays = {
    Ey_daily
    EPS_daily
    KL_daily
    Newell_daily
};

S.FluxDailyEvents = cell(3,3);
S.CFDailyEvents = cell(4,3);

for e = 1:3
    for r = 1:3

        S.FluxDailyEvents{e,r} = ...
            fluxDailyArrays{e}(idx{r},:);

    end
end

for c = 1:4
    for r = 1:3

        S.CFDailyEvents{c,r} = ...
            cfDailyArrays{c}(idx{r},:);

    end
end

S.DailyNames = DailyNames;
S.DailyOffsets = [-36 -13; -12 11; 12 35; 36 59];

S.DailyBins = dailyBins;
S.MinDailyHours = opt.MinDailyHours;

end


%% Local helper: obtain flux field using either naming convention

function x = get_flux_field(D, baseName)

hourlyName = [baseName '_hourly'];

if isfield(D, hourlyName)
    x = D.(hourlyName)(:);

elseif isfield(D, baseName)
    x = D.(baseName)(:);

else
    error('Neither %s nor %s exists in D.', ...
        hourlyName, baseName);
end

end


%% Local helper: calculate daily mean

function m = daily_mean(x, minHours)

x = x(isfinite(x));

if numel(x) >= minHours
    m = mean(x);
else
    m = NaN;
end
end
