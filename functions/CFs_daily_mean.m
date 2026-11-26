function [daily_values, daily_time, n_valid] = CFs_daily_mean(hourly_values, time)

% Calculate daily averages from hourly data.
%
% Inputs:
% hourly_values : hourly coupling-function values
% time          : corresponding hourly datetime vector
%
% Outputs:
% daily_values  : daily mean values
% daily_time    : corresponding daily timestamps
% n_valid       : number of valid hourly observations per day
%
% A daily mean is accepted only when at least 12 valid
% hourly observations are available.

% Ensure column vectors
hourly_values = hourly_values(:);
time = time(:);

% Check that the inputs have matching lengths
if numel(hourly_values) ~= numel(time)
    error('Hourly values and time must have the same length.');
end

% Group timestamps by calendar day
day_group = dateshift(time, 'start', 'day');

% Identify unique days and their group indices
[daily_time, ~, group] = unique(day_group);

% Preallocate outputs
daily_values = NaN(length(daily_time), 1);
n_valid = zeros(length(daily_time), 1);

% Calculate daily means
for i = 1:length(daily_time)

    idx = (group == i);

    values = hourly_values(idx);
    values = values(~isnan(values));

    n_valid(i) = numel(values);

    if n_valid(i) >= 12
        daily_values(i) = mean(values);
    end

end

end

