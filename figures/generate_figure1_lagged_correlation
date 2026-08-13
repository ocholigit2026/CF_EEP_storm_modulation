%% GENERATE_FIGURE1_LAGGED_CORRELATION
% Figure 1: storm-regime lagged correlation between each coupling
% function and log10(EEP flux), identifying the lag of peak correlation
% for each electron energy channel and storm regime.

addpath('../data', '../storms')

D = load_and_prepare_data();
S = identify_storms(D);

maxLag = 48;
lags = 0:maxLag;
nEnergy = numel(S.EnergyNames);
nRegime = numel(S.StormNames);
nCF     = numel(S.CFNames);

corrCF  = nan(nCF, nEnergy, nRegime, length(lags));
peakLag = nan(nCF, nEnergy, nRegime);
peakR   = nan(nCF, nEnergy, nRegime);

for iE = 1:nEnergy
    for iS = 1:nRegime
        y_full = log10(S.Flux{iE, iS});

        for iCF = 1:nCF
            x  = S.CF{iCF, iS};
            xz = (x - mean(x, 'omitnan')) ./ std(x, 'omitmissing');

            for k = 1:length(lags)
                lag = lags(k);
                xShift = xz(1:end-lag);
                yShift = y_full(1+lag:end);
                valid = isfinite(xShift) & isfinite(yShift);
                if sum(valid) > 10
                    corrCF(iCF, iE, iS, k) = corr(xShift(valid), yShift(valid), ...
                        'Rows', 'complete');
                end
            end

            [peakR(iCF, iE, iS), idxPeak] = max(corrCF(iCF, iE, iS, :));
            peakLag(iCF, iE, iS) = lags(idxPeak);
        end
    end
end

%% Plot
figure('Color', 'w', 'Position', [50 50 1500 1100]);
for iE = 1:nEnergy
    for iS = 1:nRegime
        subplot(nEnergy, nRegime, (iE-1)*nRegime + iS)
        hold on
        for iCF = 1:nCF
            plot(lags, squeeze(corrCF(iCF, iE, iS, :)), 'LineWidth', 2)
        end
        yline(0, 'k--')
        grid on; box on; axis tight
        ylim([0 0.7])

        if iE == nEnergy
            xlabel('Lag [hours]')
        end
        if iS == 1
            ylabel([S.EnergyNames{iE} ' Correlation'])
        end
        if iE == 1
            title(S.StormNames{iS})
        end
        if iE == 1 && iS == 1
            legend(S.CFNames, 'Location', 'best')
        end
    end
end

%% Peak lag / correlation summary
for iCF = 1:nCF
    fprintf('\n=== Peak %s correlations ===\n', S.CFNames{iCF});
    for iE = 1:nEnergy
        fprintf('\n%s\n', S.EnergyNames{iE});
        for iS = 1:nRegime
            fprintf('%s: r = %.3f at lag = %d h\n', S.StormNames{iS}, ...
                peakR(iCF, iE, iS), peakLag(iCF, iE, iS));
        end
    end
end
