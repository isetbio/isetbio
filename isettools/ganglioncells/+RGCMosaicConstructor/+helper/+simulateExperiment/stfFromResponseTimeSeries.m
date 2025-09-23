function [theSTFamplitudeSpectrum, theSTFphaseSpectrum, hFig] = stfFromResponseTimeSeries(...
	spatialFrequenciesCPD, theSTFresponses, spatialPhasesDegs, temporalSupportSeconds, varargin)

        % Parse input
        p = inputParser;
        p.addParameter('visualizeSinusoidalFits', false, @islogical);
        p.parse(varargin{:});
        visualizeSinusoidalFits = p.Results.visualizeSinusoidalFits;

	% Determine temporal frequency
        totalTime = temporalSupportSeconds(end)-temporalSupportSeconds(1);
        temporalFrequencyHz = 1/totalTime;
        timeSupportHR = temporalSupportSeconds(1):0.01:temporalSupportSeconds(end);

        theSTFamplitudeSpectrum = zeros(1, numel(spatialFrequenciesCPD));
        theSTFphaseSpectrum = zeros(1, numel(spatialFrequenciesCPD));

        hFig = [];

        if (visualizeSinusoidalFits)
                hFig = figure(555); clf;
                set(hFig, 'Position', [10 10 1250 1150], 'Color', [1 1 1]);
                pause(1);
                cols = 4;
                rows = round(max([20 numel(spatialFrequenciesCPD)])/cols);

                subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', rows, ...
                   'colsNum', cols, ...
                   'heightMargin',  0.04, ...
                   'widthMargin',    0.06, ...
                   'leftMargin',     0.04, ...
                   'rightMargin',    0.00, ...
                   'bottomMargin',   0.04, ...
                   'topMargin',      0.02);

                minR = 0.2*max(abs(theSTFresponses(:)));

                for iSF = 1:numel(spatialFrequenciesCPD)
                        % Retrieve the response time series for for all frames of this stimulus
                        theResponseTimeSeries = squeeze(theSTFresponses(iSF,:));

                        % Fit a sinusoid to extract the magnitude and phase of the response
                        [theFittedResponseHR, fitParams] = RGCMosaicConstructor.helper.fit.sinusoidToResponseTimeSeries(...
                                temporalSupportSeconds, theResponseTimeSeries, temporalFrequencyHz, timeSupportHR);
                        theSTFamplitudeSpectrum(1,iSF) = fitParams(1);
                        theSTFphaseSpectrum(1,iSF) = fitParams(2);

                        if (numel(temporalSupportSeconds) < 32)
                                markerSize = 200;
                                MarkerFaceAlpha =  0.8;
                                MarkerEdgeAlpha = 1.0;
                        else
                                markerSize = 100;
                                MarkerFaceAlpha =  0.5;
                                MarkerEdgeAlpha = 0.5;
                        end

                        theRow = floor((iSF-1)/cols)+1;
                        theCol = mod(iSF-1, cols)+1;
                        ax = subplot('Position', subplotPosVectors(theRow,theCol).v);
                        scatter(ax, temporalSupportSeconds*1e3, theResponseTimeSeries, markerSize, ...
                                'MarkerFaceAlpha', MarkerFaceAlpha , 'MarkerEdgeAlpha', MarkerEdgeAlpha, 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.1 0.1 0.1]);
                        hold(ax, 'on');
                        plot(ax,timeSupportHR*1e3, theFittedResponseHR, 'r-', 'LineWidth', 1.5);
                        maxR = 1.1*max([minR max(abs(theResponseTimeSeries(:))) max(abs(theFittedResponseHR(:)))]);

                        yTickLabels = -10*minR:minR:10*minR;
                        set(ax, 'XTick', 0:100:timeSupportHR(end)*1e3, 'YTick', yTickLabels, 'YTickLabels', sprintf('%2.1f\n', yTickLabels));
                        grid(ax, 'on'); box(ax, 'on');
                        set(ax, 'XColor', [0.4 0.4 0.4], 'YColor', [0.4 0.4 0.4]);
                        set(ax, 'FontSize', 14)
                        set(ax, 'YLim', maxR *[-1 1]);
                        if (theRow < rows)
                                set(ax, 'XTickLabel', {});
                        else
                                xtickangle(ax, 90);
                                xlabel(ax, 'time (msec)')
                        end

                        title(ax,sprintf('%2.2f c/deg', spatialFrequenciesCPD(iSF)));
                        drawnow
                end
        else
                parfor iSF = 1:numel(spatialFrequenciesCPD)
                        % Retrieve the response time series for for all frames of this stimulus
                        theResponseTimeSeries = squeeze(theSTFresponses(iSF,:));

                        % Fit a sinusoid to extract the magnitude and phase of the response
                        [theFittedResponseHR, fitParams] = RGCMosaicConstructor.helper.fit.sinusoidToResponseTimeSeries(...
                                temporalSupportSeconds, theResponseTimeSeries, temporalFrequencyHz, timeSupportHR);
                        theSTFamplitudeSpectrum(1,iSF) = fitParams(1);
                        theSTFphaseSpectrum(1,iSF) = fitParams(2);
                end
        end

end