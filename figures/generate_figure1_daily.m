%% GENERATE_FIGURE1_DAILY_LAGGED_CORRELATION
%
% Daily-resolution sensitivity analysis for Figure 1.
%
% Examines storm-regime lagged correlations between daily mean
% coupling functions and log10(daily mean EEP flux).
%
% Daily intervals:
% Day -1, Day 0, Day +1, Day +2, +3
%
% Lag is expressed in days.

clear;
clc;

addpath('C:\Users\mamud\Desktop\storm_modulation_scripts');

D = load_and_prepare_data_new();

S = identify_storms_new(D);

%% Settings

lags = 0:3;

nEnergy = numel(S.EnergyNames);
nRegime = numel(S.StormNames);
nCF = numel(S.CFNames);

corrCF = NaN(nCF,nEnergy,nRegime,length(lags));

peakLag = NaN(nCF,nEnergy,nRegime);
peakR = NaN(nCF,nEnergy,nRegime);

%% Calculate daily lagged correlations

for iE = 1:nEnergy

    for iS = 1:nRegime

        % Daily EEP arrays: rows represent individual storms.
        FluxEvents = S.FluxDailyEvents{iE,iS};

        % Log-transform the daily mean EEP flux.
        YEvents = log10(FluxEvents);

        for iCF = 1:nCF

            CFEvents = S.CFDailyEvents{iCF,iS};

            % Standardize coupling function within storm regime.
            xMean = mean(CFEvents(:),'omitnan');
            xStd  = std(CFEvents(:),'omitnan');

            if ~isfinite(xStd) || xStd == 0
                continue
            end

            XEvents = (CFEvents - xMean) ./ xStd;

            for k = 1:length(lags)

                lag = lags(k);

                xAll = [];
                yAll = [];

                % Calculate pairs within each storm separately.
                for iStorm = 1:size(XEvents,1)

                    x = XEvents(iStorm,:);
                    y = YEvents(iStorm,:);

                    xShift = x(1:end-lag);
                    yShift = y(1+lag:end);

                    valid = isfinite(xShift) & isfinite(yShift);

                    xAll = [xAll, xShift(valid)];
                    yAll = [yAll, yShift(valid)];

                end

                if numel(xAll) > 10

                    corrCF(iCF,iE,iS,k) = corr( ...
                        xAll(:),yAll(:),'Rows','complete');

                end
            end

            % Peak correlation and corresponding lag.
            currentCorr = squeeze(corrCF(iCF,iE,iS,:));

            if any(isfinite(currentCorr))

                [peakR(iCF,iE,iS),idxPeak] = ...
                    max(currentCorr);

                peakLag(iCF,iE,iS) = lags(idxPeak);

            end
        end
    end
end

%% Plot daily lagged correlations

figure('Color','w','Position',[50 50 1500 1100]);

for iE = 1:nEnergy

    for iS = 1:nRegime

        subplot(nEnergy,nRegime,(iE-1)*nRegime+iS);

        hold on;

        for iCF = 1:nCF

            plot(lags, ...
                squeeze(corrCF(iCF,iE,iS,:)), ...
                'LineWidth',2);

        end

        xticks([0 1 2 3])

        yline(0,'k--');

        grid on;
        box on;
        axis tight;

        ylim([0 0.9]);

        if iE == nEnergy
            xlabel('Lag [days]');
        end

        if iS == 1
            ylabel([S.EnergyNames{iE} ' Correlation']);
        end

        if iE == 1
            title(S.StormNames{iS});
        end

        if iE == 1 && iS == 1
            legend(S.CFNames,'Location','best');
        end

    end
end


% sgtitle('Daily-resolution storm-regime lagged correlations');

%% Print peak correlation summary

for iCF = 1:nCF

    fprintf('\n=== Daily peak %s correlations ===\n', ...
        S.CFNames{iCF});

    for iE = 1:nEnergy

        fprintf('\n%s\n',S.EnergyNames{iE});

        for iS = 1:nRegime

            fprintf('%s: r = %.3f at lag = %.0f days\n', ...
                S.StormNames{iS}, ...
                peakR(iCF,iE,iS), ...
                peakLag(iCF,iE,iS));

        end
    end
end
