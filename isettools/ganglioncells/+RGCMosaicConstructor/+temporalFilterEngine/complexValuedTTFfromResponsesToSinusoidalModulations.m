%
% RGCMosaicConstructor.temporalFilterEngine.complexValuedTTFfromResponsesToSinusoidalModulations
%

function [theComplexTTF, temporalFrequencySupportHz] = complexValuedTTFfromResponsesToSinusoidalModulations(...
        temporalFrequenciesExamined, ...
        temporalSupportSecondsAllConditions, ...
        responsesAllConditions, ...
        allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries, ...
        visualizeSinusoidalFits, ...
        stimulusShape, ...
        theSinudoidalFitColor)


    theTTFamplitude = zeros(1, numel(temporalFrequenciesExamined));
    theTTFphaseDegs = zeros(1, numel(temporalFrequenciesExamined));

    
    temporalSupport1Msec = cell(1, numel(temporalFrequenciesExamined));
    validTemporalSupportSeconds = cell(1, numel(temporalFrequenciesExamined));
    validResponses = cell(1, numel(temporalFrequenciesExamined));

    theFittedResponses = cell(1, numel(temporalFrequenciesExamined));
    fittedParams = [];
    meanResponses = zeros(1, numel(temporalFrequenciesExamined));

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


        assert(isempty(find(isnan(responseForThisTF(:)))), 'did not remove all nan part of the response');
        assert(isempty(find(isnan(temporalSupportSecondsForThisTF(:)))), 'did not remove all nan part of the response temporal support');

        
        temporalSupport1Msec{iTF} = 0:1/1000:(temporalSupportSecondsForThisTF(end)+dt);

        validTemporalSupportSeconds{iTF} = temporalSupportSecondsForThisTF;
        validResponses{iTF} = responseForThisTF;

        [theFittedResponses{iTF}, fittedParams(iTF,:), meanResponses(iTF)] = RGCMosaicConstructor.helper.fit.sinusoidToResponseTimeSeries(...
            temporalSupportSecondsForThisTF, ...
            responseForThisTF, ...
            temporalFrequency, ...
            temporalSupport1Msec{iTF}, ...
            'allowOffset', allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries);

        theTTFamplitude(iTF) = fittedParams(iTF,1);
        theTTFphaseDegs(iTF) = fittedParams(iTF,2);

    end % iTF

    if (visualizeSinusoidalFits)

        for iTF = 1:numel(temporalFrequenciesExamined)
        
           

            hFig = figure(2); clf;
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure', ...
                'darkScheme', true);


            theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
            ax = theAxes{1,1};

            temporalSupportSecondsForThisTF = validTemporalSupportSeconds{iTF};
            responseForThisTF = validResponses{iTF};
            temporalFrequency = temporalFrequenciesExamined(iTF);

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
            YLims = 2*[-1 1];

            dT = temporalSupportSecondsForThisTF(2)-temporalSupportSecondsForThisTF(1);
            skip = max([1 round(1/temporalFrequency/16/dT)]);

            plot(temporalSupport1Msec{iTF}*1e3, temporalSupport1Msec{iTF}*0, '-', 'Color', ff.labelColor, 'LineWidth', 1.5);
            hold(ax, 'on')
 

            % We subtract the mean response, which is stronger at high TFs.
            % We believe this is so because there is not enough 'warm up' time in the
            % high TFs to reach zero baseline
            meanResponse = meanResponses(iTF);

            scatter(ax,temporalSupportSecondsForThisTF(1:skip:end)*1e3, responseForThisTF(1:skip:end)-meanResponse, 18*18, ...
                'MarkerFaceColor', theSinudoidalFitColor, 'MarkerEdgeColor', theSinudoidalFitColor* 0.5, 'MarkerFaceAlpha', 1.0, 'LineWidth', 1.5);
            
            plot(ax, temporalSupport1Msec{iTF}*1e3, theFittedResponses{iTF}-meanResponse, 'k-', 'LineWidth', 4);
            plot(ax, temporalSupport1Msec{iTF}*1e3, theFittedResponses{iTF}-meanResponse, 'w-', 'LineWidth', 2);
            yTickLabels = {'-2' '' '-1' '' '0' '' '+1' '' '+2'};
            set(ax, 'XLim', XLims, 'YLim', YLims, 'XTick', xTicks, 'YTick', -2:0.5:2, 'YTickLabel', yTickLabels, 'XTickLabel', xTickLabels);
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
        
        end % iTF
    end % visualizeSinusoidalFits


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
