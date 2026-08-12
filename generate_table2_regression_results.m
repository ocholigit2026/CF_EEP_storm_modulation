%% GENERATE_TABLE2_REGRESSION_RESULTS
% Reproduces Table 2: storm-regime regression analysis of normalized
% coupling function strength vs. log10(EEP flux), for each electron
% energy channel, storm regime, and coupling function.
%
%   Model: log10(Flux) = a + b * (z-scored coupling function)

addpath('../data', '../storms')

D = load_and_prepare_data();
S = identify_storms(D);

nEnergy = numel(S.EnergyNames);
nRegime = numel(S.StormNames);
nCF     = numel(S.CFNames);

Results = cell(nEnergy * nRegime * nCF, 11);
k = 1;

for iE = 1:nEnergy
    for iS = 1:nRegime
        y = log10(S.Flux{iE, iS});

        for iCF = 1:nCF
            x = S.CF{iCF, iS};
            x = (x - mean(x, 'omitnan')) ./ std(x, 'omitmissing');

            valid = isfinite(x) & isfinite(y);
            mdl = fitlm(x(valid), y(valid));

            a  = mdl.Coefficients.Estimate(1);
            b  = mdl.Coefficients.Estimate(2);
            CI = coefCI(mdl);

            Results(k, :) = { ...
                S.EnergyNames{iE}, ...
                S.StormNames{iS}, ...
                S.CFNames{iCF}, ...
                a, b, ...
                sprintf('[%.3f, %.3f]', CI(1,1), CI(1,2)), ...
                sprintf('[%.3f, %.3f]', CI(2,1), CI(2,2)), ...
                mdl.Rsquared.Ordinary, ...
                mdl.RMSE, ...
                mdl.Coefficients.pValue(2), ...
                sprintf('EEP = %.3f + %.3f x CF', a, b)};
            k = k + 1;
        end
    end
end

RegressionTable = cell2table(Results, 'VariableNames', ...
    {'Energy', 'Storm', 'CouplingFunction', 'Intercept_a', 'Slope_b', ...
     'CI95_a', 'CI95_b', 'Rsquared', 'RMSE', 'pValue', 'RegressionEquation'});

disp('Table 2. Storm-regime regression results.')
disp(RegressionTable)

% writetable(RegressionTable, 'Storm_Regression_Results.xlsx');  % optional export
