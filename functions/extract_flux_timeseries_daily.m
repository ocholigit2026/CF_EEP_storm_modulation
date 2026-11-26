function [flux_hourly, flux_daily, daily_time, n_valid] = ...
    extract_flux_timeseries_daily(flux_3D, t_hr)

% Dimensions of flux_3D:
% 1 = time, 2 = MLT, 3 = L-shell
%
% Spatial averaging:
% First average over MLT, then over L-shell indices 13:21.

% Step 1: Average over MLT
flux_ts = squeeze(mean(permute(flux_3D, [2 3 1]), ...
    'omitmissing'));

% Step 2: Average over selected L-shell bins
flux_hourly = squeeze(mean(flux_ts(13:21, :), ...
    'omitmissing'))';

% Step 3: Assign each hourly observation to its calendar day
day_group = dateshift(t_hr(:), 'start', 'day');

[daily_time, ~, group] = unique(day_group);

% Step 4: Calculate daily means and count valid hourly values
flux_daily = NaN(length(daily_time), 1);
n_valid = zeros(length(daily_time), 1);

for i = 1:length(daily_time)

    idx = (group == i);

    valid_values = flux_hourly(idx);
    valid_values = valid_values(~ismissing(valid_values));

    n_valid(i) = numel(valid_values);

    % Accept daily mean only if at least 12 hourly values exist
    if n_valid(i) >= 12
        flux_daily(i) = mean(valid_values, 'omitmissing');
    end

end

end

