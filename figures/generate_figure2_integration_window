%% GENERATE_FIGURE2_INTEGRATION_WINDOW
% Figure 2: correlation between log10(EEP flux) and coupling functions
% integrated over progressively increasing time windows, identifying the
% characteristic integration window for each storm regime and energy
% channel.

addpath('../data', '../storms')

D = load_and_prepare_data();
S = identify_storms(D);

windows = 1:48;
nEnergy = numel(S.EnergyNames);
nRegime = numel(S.StormNames);
nCF     = numel(S.CFNames);

PeakWindow = zeros(nCF, nRegime, nEnergy);
PeakCorr   = zeros(nCF, nRegime, nEnergy);

figure('Color', 'w', 'Position', [100 50 1200 900])
for iE = 1:nEnergy
    for iS = 1:nRegime
        subplot(nEnergy, nRegime, (iE-1)*nRegime + iS)
        hold on

        y = log10(S.Flux{iE, iS});

        for iCF = 1:nCF
            CF = S.CF{iCF, iS};
            corrVals = nan(size(windows));
            for w = windows
                CFint = movsum(CF, [w 0], 'omitnan');
                corrVals(w) = corr(y, CFint, 'Rows', 'complete');
            end

            [PeakCorr(iCF, iS, iE), idxPeak] = max(corrVals);
            PeakWindow(iCF, iS, iE) = windows(idxPeak);

            plot(windows, corrVals, 'LineWidth', 2)
        end

        axis tight; grid on; box on
        ylim([0 0.52])

        if iE == 1
            title(S.StormNames{iS})
        end
        if iS == 1
            ylabel([S.EnergyNames{iE} ' Correlation'])
        end
        if iE == nEnergy
            xlabel('Integration Window [hours]')
        end
        if iE == 1 && iS == nRegime
            legend(S.CFNames, 'Location', 'best')
        end
    end
end
