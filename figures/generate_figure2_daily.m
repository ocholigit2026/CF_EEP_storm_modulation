%% GENERATE_FIGURE2_DAILY_INTEGRATION_WINDOW
%
% Daily-resolution sensitivity analysis for Figure 2.
%
% Examines correlations between log10(daily mean EEP flux)
% and daily coupling functions integrated over progressively
% increasing windows.
%
% Integration windows are expressed in days.

clear;
clc;

addpath('C:\Users\mamud\Desktop\storm_modulation_scripts');

D = load_and_prepare_data_new();

S = identify_storms_new(D);

%% Settings

windows = 1:4;

nEnergy = numel(S.EnergyNames);
nRegime = numel(S.StormNames);
nCF = numel(S.CFNames);

PeakWindow = NaN(nCF,nRegime,nEnergy);
PeakCorr = NaN(nCF,nRegime,nEnergy);

corrResults = NaN(nCF,nRegime,nEnergy,length(windows));

%% Calculate daily integrated correlations

figure('Color','w','Position',[100 50 1200 900]);

for iE = 1:nEnergy

    for iS = 1:nRegime

        subplot(nEnergy,nRegime,(iE-1)*nRegime+iS);

        hold on;

        % Daily EEP flux, log-transformed after averaging.
        FluxEvents = S.FluxDailyEvents{iE,iS};

        YEvents = log10(FluxEvents);

        for iCF = 1:nCF

            CFEvents = S.CFDailyEvents{iCF,iS};

            corrVals = NaN(size(windows));

            for iw = 1:length(windows)

                w = windows(iw);

                xAll = [];
                yAll = [];

                % Integrate separately within each storm.
                for iStorm = 1:size(CFEvents,1)

                    CF = CFEvents(iStorm,:);
                    y = YEvents(iStorm,:);

                    % Trailing integration over w daily values.
                    CFint = movsum(CF,[w-1 0],'omitnan');

                    valid = isfinite(CFint) & isfinite(y);

                    xAll = [xAll, CFint(valid)];
                    yAll = [yAll, y(valid)];

                end

                if numel(xAll) > 10

                    corrVals(iw) = corr( ...
                        xAll(:),yAll(:),'Rows','complete');

                end
            end

            corrResults(iCF,iS,iE,:) = corrVals;

            % Peak correlation and corresponding integration window.
            if any(isfinite(corrVals))

                [PeakCorr(iCF,iS,iE),idxPeak] = max(corrVals);

                PeakWindow(iCF,iS,iE) = windows(idxPeak);

            end

            plot(windows,corrVals,'LineWidth',2);

        end

        xticks([1 2 3 4])

        axis tight;
        grid on;
        box on;

        ylim([0 0.9]);

        if iE == 1
            title(S.StormNames{iS});
        end

        if iS == 1
            ylabel([S.EnergyNames{iE} ' Correlation']);
        end

        if iE == nEnergy
            xlabel('Integration Window [days]');
        end

        if iE == 1 && iS == 1
            legend(S.CFNames,'Location','best');
        end

    end
end

% sgtitle('Daily-resolution storm-regime integration-window correlations');

%% Print peak correlation summary

for iCF = 1:nCF

    fprintf('\n=== Daily integration results: %s ===\n', ...
        S.CFNames{iCF});

    for iE = 1:nEnergy

        fprintf('\n%s\n',S.EnergyNames{iE});

        for iS = 1:nRegime

            fprintf('%s: r = %.3f at window = %.0f days\n', ...
                S.StormNames{iS}, ...
                PeakCorr(iCF,iS,iE), ...
                PeakWindow(iCF,iS,iE));

        end
    end
end
