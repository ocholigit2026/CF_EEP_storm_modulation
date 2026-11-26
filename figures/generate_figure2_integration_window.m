%% GENERATE_FIGURE2_INTEGRATION_WINDOW
% Figure 2: correlation between log10(EEP flux) and coupling functions
% integrated over progressively increasing time windows, identifying the
% characteristic integration window for each storm regime and energy
% channel.

addpath('C:\Users\mamud\Desktop\storm_modulation_scripts')

D = load_and_prepare_data_new();
S = identify_storms_new(D);

windows = 1:48;
nEnergy = numel(S.EnergyNames);
nRegime = numel(S.StormNames);
nCF     = numel(S.CFNames);

PeakWindow = nan(nCF, nRegime, nEnergy);
PeakCorr   = nan(nCF, nRegime, nEnergy);

Rint = NaN(4, 3, 3, length(windows));

figure('Color', 'w', 'Position', [100 50 1200 900])
for iE = 1:nEnergy
    for iS = 1:nRegime
        subplot(nEnergy, nRegime, (iE-1)*nRegime + iS)
        hold on

        Yevents = log10(S.FluxEvents{iE, iS});

        for iCF = 1:nCF
            Xevents = S.CFEvents{iCF, iS};
            for k = 1:length(windows)
                w = windows(k);
                Xpool = [];
                Ypool = [];
                for s = 1:size(Xevents, 1)
                    x = Xevents(s, :);
                    y = Yevents(s, :);
                    % Integrate within the individual storm
                    xint = movsum(x, [w-1 0], 'omitnan');
                    % Retain only complete integration windows
                    validCols = w:length(x);
                    xs = xint(validCols);
                    ys = y(validCols);
                    valid = isfinite(xs) & isfinite(ys);
                    Xpool = [Xpool; xs(valid)];
                    Ypool = [Ypool; ys(valid)];
                end
                if numel(Xpool) > 2
                    Rint(iCF, iE, iS, k) = corr(Xpool(:), Ypool(:), "type","Pearson");
                end
            end

            [PeakCorr(iCF, iE, iS), idxPeak] = max(Rint(iCF, iE, iS, :));
            PeakWindow(iCF, iE, iS) = windows(idxPeak);

            % corrVals = nan(size(windows));
            % for w = windows
            %     CFint = movsum(CF, [w 0], 'omitnan');
            %     corrVals(w) = corr(y, CFint, 'Rows', 'complete');
            % end

            % [PeakCorr(iCF, iS, iE), idxPeak] = max(corrVals);
            % PeakWindow(iCF, iS, iE) = windows(idxPeak);
            % 
            plot(windows, squeeze(Rint(iCF, iE, iS, :)), 'LineWidth', 2)
        end
        
        xticks([1 10 20 30 40 48])
        axis tight; grid on; box on
        ylim([0 0.8])

        if iE == 1
            title(S.StormNames{iS})
        end
        if iS == 1
            ylabel([S.EnergyNames{iE} ' Correlation'])
        end
        if iE == nEnergy
            xlabel('Integration Window [hours]')
        end
        if iE == 1 && iS == 1
            legend(S.CFNames, 'Location', 'best')
        end
    end
end

%% Peak window / correlation summary
for iCF = 1:nCF
    fprintf('\n=== Peak %s correlations ===\n', S.CFNames{iCF});
    for iE = 1:nEnergy
        fprintf('\n%s\n', S.EnergyNames{iE});
        for iS = 1:nRegime
            fprintf('%s: r = %.3f at window = %d h\n', S.StormNames{iS}, ...
                PeakCorr(iCF, iE, iS), PeakWindow(iCF, iE, iS));
        end
    end
end