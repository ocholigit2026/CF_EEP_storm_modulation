%% GENERATE_TABLE1_STORM_DISTRIBUTION
% Reproduces Table 1: the distribution of geomagnetic storms identified
% from Dst index data (2009-2019) by intensity regime (weak / moderate /
% intense), as reported in the manuscript.

addpath('../data', '../storms')

D = load_and_prepare_data();
S = identify_storms(D);

counts = [sum(S.stormSize == 1); ...
          sum(S.stormSize == 2); ...
          sum(S.stormSize == 3)];
total = sum(counts);

Table1 = table([counts; total], ...
    'VariableNames', {'NumberOfStorms'}, ...
    'RowNames', [S.StormNames, {'Total'}]);

disp('Table 1. Distribution of selected storms according to intensity.')
disp(Table1)
