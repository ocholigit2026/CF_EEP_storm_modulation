function S = identify_storms(D, varargin)
% IDENTIFY_STORMS Detect geomagnetic storms from the Dst index and
% segment coupling-function / EEP-flux time series into storm-centered
% epochs, classified by storm intensity regime.
%
%   S = IDENTIFY_STORMS(D) identifies storm centers in D.Dst using a
%   threshold of Dst <= -30 nT with a minimum 96-hour separation between
%   storms, extracts a window from PRE hours before to POST hours after
%   each storm center, and classifies each storm as weak / moderate /
%   intense based on its minimum Dst.
%
% Name-value options:
%   'Threshold'     Dst threshold for storm identification (default -30)
%   'MinSeparation' Minimum hours between storm centers   (default 96)
%   'PreWindow'     Hours before storm center to include  (default 48)
%   'PostWindow'    Hours after storm center to include   (default 72)
%
% Fields of S:
%   stormTimes            - datetime of each storm center
%   stormSize             - 1=weak, 2=moderate, 3=intense (NaN=unclassified)
%   CF                    - 4x3 cell array {coupling function, storm regime}
%                            of vectorized, storm-segmented values.
%                            Rows: Ey, EPS, KL, Newell. Columns: weak,
%                            moderate, intense.
%   Flux                  - 3x3 cell array {energy channel, storm regime}
%                            of vectorized, storm-segmented flux values
%                            (native units, not log-transformed).
%                            Rows: >30, >100, >300 keV.
%   CFNames, EnergyNames, StormNames - labels for CF/Flux rows/columns

p = inputParser;
addParameter(p, 'Threshold', -30);
addParameter(p, 'MinSeparation', 96);
addParameter(p, 'PreWindow', 48);
addParameter(p, 'PostWindow', 72);
parse(p, varargin{:});
opt = p.Results;

% --- Step 1: find Dst minima at or below threshold ---
[~, locs] = findpeaks(-D.Dst);
stormCandidates = locs(D.Dst(locs) <= opt.Threshold);

% --- Step 2: enforce minimum separation between storm centers ---
stormCenters = select_peak_events(stormCandidates, -D.Dst, opt.MinSeparation);

pre  = opt.PreWindow;
post = opt.PostWindow;
T = pre + post + 1;
N = length(stormCenters);

Dst_st     = NaN(N, T);
flux30_st  = NaN(N, T);
flux100_st = NaN(N, T);
flux300_st = NaN(N, T);
eps_st     = NaN(N, T);
KL_st      = NaN(N, T);
Ey_st      = NaN(N, T);
Newell_st  = NaN(N, T);
stormSize  = NaN(N, 1);

for i = 1:N
    center = stormCenters(i);
    start_idx = center - pre;
    end_idx   = center + post;

    if start_idx < 1 || end_idx > length(D.Dst)
        continue  % skip storms too close to the data boundary
    end

    Dst_st(i,:)     = D.Dst(start_idx:end_idx);
    flux30_st(i,:)  = D.flux30(start_idx:end_idx);
    flux100_st(i,:) = D.flux100(start_idx:end_idx);
    flux300_st(i,:) = D.flux300(start_idx:end_idx);
    eps_st(i,:)     = D.EPS(start_idx:end_idx);
    KL_st(i,:)      = D.KL(start_idx:end_idx);
    Ey_st(i,:)      = D.Ey(start_idx:end_idx);
    Newell_st(i,:)  = D.Newell(start_idx:end_idx);

    Dst_min = min(Dst_st(i,:));
    if Dst_min <= -100
        stormSize(i) = 3;  % intense
    elseif Dst_min <= -50
        stormSize(i) = 2;  % moderate
    elseif Dst_min <= opt.Threshold
        stormSize(i) = 1;  % weak
    end
end

idx = {stormSize == 1, stormSize == 2, stormSize == 3};

S.stormTimes  = D.time(stormCenters);
S.stormSize   = stormSize;
S.StormNames  = {'Weak', 'Moderate', 'Intense'};
S.EnergyNames = {'>30 keV', '>100 keV', '>300 keV'};
S.CFNames     = {'E_y', '\epsilon', 'Kan-Lee', 'Newell'};

fluxArrays = {flux30_st, flux100_st, flux300_st};
S.Flux = cell(3, 3);
for e = 1:3
    for r = 1:3
        S.Flux{e, r} = reshape(fluxArrays{e}(idx{r}, :), [], 1);
    end
end

cfArrays = {Ey_st, eps_st, KL_st, Newell_st};
S.CF = cell(4, 3);
for c = 1:4
    for r = 1:3
        S.CF{c, r} = reshape(cfArrays{c}(idx{r}, :), [], 1);
    end
end

end
