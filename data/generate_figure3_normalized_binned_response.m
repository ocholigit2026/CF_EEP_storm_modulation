%% GENERATE_FIGURE4_NATIVE_BINNED_RESPONSE
% Figure 4: storm-regime binned log10(EEP flux) response as a function
% of each coupling function in its native (non-normalized) units, using
% common quantile bins across storm regimes, with interquartile-range
% shading and lines split by electron energy channel.

addpath('../data', '../storms')

D = load_and_prepare_data();
S = identify_storms(D);

nEnergy = numel(S.EnergyNames);
nRegime = numel(S.StormNames);
nCF     = numel(S.CFNames);
nBins   = 20;
minPoints = 5;

LineColours = lines(nEnergy);
hLine = gobjects(nEnergy, 1);

figure('Color', 'w', 'Position', [50 30 1300 1200]);

for iCF = 1:nCF
    % Common quantile bin edges built from all storm regimes combined
    Xall = vertcat(S.CF{iCF, :});
    Xall = Xall(isfinite(Xall));

    edges = unique(quantile(Xall, linspace(0, 1, nBins+1)));
    if length(edges) < nBins + 1
        edges = linspace(min(Xall), max(Xall), nBins+1);
    end
    nEdges = length(edges) - 1;

    for iS = 1:nRegime
        subplot(nCF, nRegime, (iCF-1)*nRegime + iS)
        hold on; box on

        for iE = 1:nEnergy
            x = S.CF{iCF, iS};
            y = log10(S.Flux{iE, iS});

            valid = isfinite(x) & isfinite(y);
            x = x(valid);
            y = y(valid);

            xMedian = nan(nEdges, 1);
            yMedian = nan(nEdges, 1);
            y25     = nan(nEdges, 1);
            y75     = nan(nEdges, 1);

            for b = 1:nEdges
                if b < nEdges
                    inBin = x >= edges(b) & x < edges(b+1);
                else
                    inBin = x >= edges(b) & x <= edges(b+1);
                end
                if sum(inBin) >= minPoints
                    xMedian(b) = median(x(inBin));
                    yMedian(b) = median(y(inBin));
                    y25(b) = prctile(y(inBin), 25);
                    y75(b) = prctile(y(inBin), 75);
                end
            end

            keep = isfinite(xMedian) & isfinite(yMedian);
            xMedian = xMedian(keep); yMedian = yMedian(keep);
            y25 = y25(keep); y75 = y75(keep);

            fill([xMedian; flipud(xMedian)], [y25; flipud(y75)], ...
                LineColours(iE,:), 'FaceAlpha', 0.18, 'EdgeColor', 'none');
            h = plot(xMedian, yMedian, 'Color', LineColours(iE,:), 'LineWidth', 2);
            if iCF == 1 && iS == 1
                hLine(iE) = h;
            end
        end

        grid on
        if iCF == 1
            title(sprintf('%s (N=%d)', S.StormNames{iS}, length(S.CF{iCF,iS})))
        end
        if iS == 1
            ylabel('log_{10} Flux', 'FontWeight', 'bold')
        end
        xlabel(S.CFNames{iCF}, 'FontWeight', 'bold')
        axis tight
        set(gca, 'FontSize', 10, 'LineWidth', 1)
    end
end

lgd = legend(hLine, S.EnergyNames, 'Orientation', 'horizontal');
lgd.Units = 'normalized';
lgd.Position = [0.36 0.955 0.28 0.03];
lgd.Box = 'off';
lgd.FontSize = 11;

% Sync y-axis limits across all subplots
ax = findall(gcf, 'Type', 'axes');
yl = cell2mat(get(ax, 'YLim'));
ylim_global = [min(yl(:,1)) max(yl(:,2))];
set(ax, 'YLim', ylim_global)
