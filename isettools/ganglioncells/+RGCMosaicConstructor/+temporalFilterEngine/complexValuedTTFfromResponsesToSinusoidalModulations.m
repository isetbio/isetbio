%
% RGCMosaicConstructor.temporalFilterEngine.complexValuedTTFfromResponsesToSinusoidalModulations
%

function [theComplexTTF, temporalFrequencySupportHz] = complexValuedTTFfromResponsesToSinusoidalModulations(...
        temporalFrequenciesExamined, ...
        temporalSupportSecondsAllConditions, ...
        responsesAllConditions, ...
        allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries, ...
        phaseComputationParams, ...
        visualizeSinusoidalFits, ...
        stimulusShape, ...
        theColor)


    theTTFamplitude = zeros(1, numel(temporalFrequenciesExamined));
    theTTFphaseDegs = zeros(1, numel(temporalFrequenciesExamined));

    


    for iTF = 1:numel(temporalFrequenciesExamined)

        % Retrieve the photocurrent response data for this TF
        temporalSupportSecondsForThisTF = squeeze(temporalSupportSecondsAllConditions(iTF,:));
        responseForThisTF = squeeze(responsesAllConditions(iTF,:));
        temporalFrequency = temporalFrequenciesExamined(iTF);

        % find valid time bins for this TF (we have stored 1 period of the
        % respons, so each TF response has a different length)
        nanIndices = find(isnan(responseForThisTF));
        if (isempty(nanIndices))
            theTimeBins = 1:size(responsesAllConditions,2);
        else
            theTimeBins = 1:(nanIndices(1)-1);
        end

        temporalSupportSecondsForThisTF = temporalSupportSecondsForThisTF(theTimeBins);
        responseForThisTF = responseForThisTF(theTimeBins);
        dt = temporalSupportSecondsForThisTF(2)-temporalSupportSecondsForThisTF(1);

        if (~isempty(phaseComputationParams))
            responseDelayBins = round((phaseComputationParams.delayMilliSeconds/1000)/dt);
            responseForThisTF = circshift(responseForThisTF, responseDelayBins);
        end

        assert(isempty(find(isnan(responseForThisTF(:)))), 'did not remove all nan part of the response');
        assert(isempty(find(isnan(temporalSupportSecondsForThisTF(:)))), 'did not remove all nan part of the response temporal support');

        
        temporalSupport1Msec = 0:1/1000:(temporalSupportSecondsForThisTF(end)+dt);

        [theFittedResponse, fittedParams, meanResponse] = RGCMosaicConstructor.helper.fit.sinusoidToResponseTimeSeries(...
            temporalSupportSecondsForThisTF, ...
            responseForThisTF, ...
            temporalFrequency, ...
            temporalSupport1Msec, ...
            'allowOffset', allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries);

        theTTFamplitude(iTF) = fittedParams(1);
        theTTFphaseDegs(iTF) = fittedParams(2);

     

        if (visualizeSinusoidalFits)
            hFig = figure(2); clf;
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
            theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
            ax = theAxes{1,1};

            if (temporalSupportSecondsForThisTF(end) > 1)
                dTick = 0.4;
            elseif (temporalSupportSecondsForThisTF(end) > 0.5)
                dTick = 0.2;
            elseif (temporalSupportSecondsForThisTF(end) > 0.25)
                dTick = 0.1;
            elseif (temporalSupportSecondsForThisTF(end) > 0.125)
                dTick = 0.05;
            elseif (temporalSupportSecondsForThisTF(end) > 0.05)
                dTick = 0.025;
            elseif (temporalSupportSecondsForThisTF(end) > 0.03)
                dTick = 0.01;
            elseif (temporalSupportSecondsForThisTF(end) > 0.01)
                dTick = 0.005;
            else
                dTick = 0.002;
            end
            xTicks = (0:dTick:2.0)*1e3;
            xTickLabels = sprintf('%.0f\n', xTicks);

            XLims = [0 temporalSupportSecondsForThisTF(end)]*1e3;
            YLims = [-1 1];

            dT = temporalSupportSecondsForThisTF(2)-temporalSupportSecondsForThisTF(1);
            skip = max([1 round(1/temporalFrequency/16/dT)]);

            plot(temporalSupport1Msec*1e3, temporalSupport1Msec*0, 'k-', 'LineWidth', 1.5);
            hold(ax, 'on')
            scatter(ax,temporalSupportSecondsForThisTF(1:skip:end)*1e3, responseForThisTF(1:skip:end)-meanResponse, 300, ...
                'MarkerFaceColor', theColor, 'MarkerEdgeColor', theColor * 0.5, 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);
            
            plot(ax, temporalSupport1Msec*1e3, theFittedResponse-meanResponse, 'k-', 'LineWidth', 4);
            plot(ax, temporalSupport1Msec*1e3, theFittedResponse-meanResponse, 'y-', 'LineWidth', 2);
            yTickLabels = {'-1.0' '' '-.50' '' '0' '' '+.50' '' '+1.0'};
            set(ax, 'XLim', XLims, 'YLim', YLims, 'XTick', xTicks, 'YTick', -1:0.25:1, 'YTickLabel', yTickLabels, 'XTickLabel', xTickLabels);
            grid(ax, 'on');
            xlabel(ax,'time (msec)');
            ylabel(ax,'pAmps');
            title(ax, sprintf('%s, TF = %2.2fHz', stimulusShape, temporalFrequency));

            % Finalize figure using the Publication-Ready format
            PublicationReadyPlotLib.applyFormat(ax,ff);
            PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

            % Export
            visualizationPDFfileName = sprintf('%s_TF_%2.2fHz', stimulusShape, temporalFrequency);
            exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
            pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
            theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

            % Generate the path if we need to
            RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);

            thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
            NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');
        end % visualizeSinusoidalFits
    end % iTF


    % Phase in radians
    theTTFphaseRadians = theTTFphaseDegs/180*pi;


    % Complex TTF from amplitude and phase
    theComplexTTF = theTTFamplitude .* exp(-1i * theTTFphaseRadians);


    % Adjust TTF for missing 0 Hz point.
    [theComplexTTF, temporalFrequencySupportHz] = RGCMosaicConstructor.temporalFilterEngine.adjustTTFtoDealWithMissingTTFsampleAt0Hz(...
        theComplexTTF, temporalFrequenciesExamined);

    % Clean up the TTF by removing phase outliers
    theComplexTTF = RGCMosaicConstructor.temporalFilterEngine.removePhaseOutliersFromTTF(...
        temporalFrequencySupportHz, theComplexTTF);

end
