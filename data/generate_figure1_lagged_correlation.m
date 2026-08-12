%% GENERATE_FIGURE3_NORMALIZED_BINNED_RESPONSE
% Figure 3: storm-regime binned log10(EEP flux) response as a function
% of normalized (z-scored) coupling function strength, using
% equal-population quantile bins, shown separately for each electron
% energy channel.

addpath('../data', '../storms')

D = load_and_prepare_data();
S = identify_storms(D);

nEnergy = numel(S.EnergyNames);
nRegime = numel(S.StormNames);
nCF     = numel(S.CFNames);
nbins   = 12;

CFcolor = [0.0000 0.4470 0.7410   % Ey
           0.8500 0.3250 0.0980   % epsilon
           0.8500 0.6500 0.1300   % Kan-Lee
           0.4940 0.1840 0.5560]; % Newell

yLimits = {[2.9 6.0], [2.4 4.5], [2.0 3.5]};

figure('Color', 'w', 'Position', [100 50 1200 1000]);
hLeg = gobjects(nCF, 1);

for iE = 1:nEnergy
    for iS = 1:nRegime
        subplot(nEnergy, nRegime, (iE-1)*nRegime + iS)
        hold on

        for iCF = 1:nCF
            x = S.CF{iCF, iS};
            x = (x - mean(x, 'omitnan')) ./ std(x, 'omitmissing');
            y = log10(S.Flux{iE, iS});

            valid = isfinite(x) & isfinite(y);
            x = x(valid);
            y = y(valid);

            edges = unique(quantile(x, linspace(0, 1, nbins+1)));
            centres = nan(1, length(edges)-1);
            ymean   = nan(1, length(edges)-1);

            for b = 1:length(edges)-1
                if b < length(edges)-1
                    inBin = x >= edges(b) & x < edges(b+1);
                else
                    inBin = x >= edges(b) & x <= edges(b+1);
                end
                centres(b) = mean([edges(b), edges(b+1)]);
                if any(inBin)
                    ymean(b) = median(y(inBin), 'omitnan');
                end
            end

            p = plot(centres, ymean, '-o', 'Color', CFcolor(iCF,:), ...
                'MarkerFaceColor', CFcolor(iCF,:), 'LineWidth', 2, 'MarkerSize', 5);
            if iE == 1 && iS == 1
                hLeg(iCF) = p;
            end
        end

        grid on; box on
        ylim(yLimits{iE})

        if iE == 1
            title(S.StormNames{iS}, 'FontWeight', 'bold')
        end
        if iS == 1
            ylabel(['log_{10} Flux ', S.EnergyNames{iE}])
        end
        if iE == nEnergy && iS == 2
            xlabel('Relative Coupling strength')
        end
        set(gca, 'FontSize', 11)
    end
end

legend(hLeg, S.CFNames, 'Location', 'southeast', 'Box', 'on', 'FontSize', 9)
