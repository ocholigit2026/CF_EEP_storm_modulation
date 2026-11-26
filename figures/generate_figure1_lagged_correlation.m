%% GENERATE_FIGURE1_LAGGED_CORRELATION
% Figure 1: storm-regime lagged correlation between each coupling
% function and log10(EEP flux), identifying the lag of peak correlation
% for each electron energy channel and storm regime.

addpath('C:\Users\mamud\Desktop\storm_modulation_scripts')

D = load_and_prepare_data_new();
S = identify_storms_new(D);

maxLag = 48;
lags = 0:maxLag;
nEnergy = numel(S.EnergyNames);
nRegime = numel(S.StormNames);
nCF     = numel(S.CFNames);

corrCF  = nan(nCF, nEnergy, nRegime, length(lags));
peakLag = nan(nCF, nEnergy, nRegime);
peakR   = nan(nCF, nEnergy, nRegime);
Rlag  = nan(nCF, nEnergy, nRegime, length(lags));

for iE = 1:nEnergy
    for iS = 1:nRegime
        yz = log10(S.FluxEvents{iE, iS});

        for iCF = 1:nCF
            x  = S.CFEvents{iCF, iS};
            xz = (x - mean(x, 'omitnan')) ./ std(x, 'omitmissing');

            for k = 1:length(lags)
                lag = lags(k);
                Xpool = [];
                Ypool = [];
                for s = 1:size(xz,1)
                    xx = xz(s, :);
                    yy = yz(s, :);
                    if lag == 0
                        xs = xx;
                        ys = yy;
                    else
                        xs = xx(1:end-lag);
                        ys = yy(1+lag:end);
                    end
                    valid = isfinite(xs) & isfinite(ys);
                    Xpool = [Xpool; xs(valid)];
                    Ypool = [Ypool; ys(valid)];
                end
                if numel(Xpool) > 2
                    Rlag(iCF, iE, iS, k) = corr(Xpool(:), Ypool(:), ...
                        'type', 'Pearson');
                end

            end

            [peakR(iCF, iE, iS), idxPeak] = max(Rlag(iCF, iE, iS, :));
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
            plot(lags, squeeze(Rlag(iCF, iE, iS, :)), 'LineWidth', 2)
        end
        yline(0, 'k--')
        grid on; box on; axis tight
        ylim([0 0.6])

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
