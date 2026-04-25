%
%RGCMosaicAnalyzer.compute.MRGCtemporalFiltersFromPhotocurrentsBasedTTF
%
%
function MRGCtemporalFiltersFromPhotocurrentsBasedTTF(...
    temporalFilterSynthesisMethod, ...
    targetCellImpulseResponseSource, ...
    stimulusShape, ...
    allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries, ...
    theTargetRGCwithIndex, ...
    theMRGCMosaicTTFResponsesFullFileName, ...
    theAnalyzedTTFsFullFileName, ...
    visualizeSinusoidalFits)


    % Load the measured TTFs
    load(theMRGCMosaicTTFResponsesFullFileName, ...
        'theMRGCMosaic', 'stimParams', 'TTFparamsStruct', ...
        'computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex', ...
        'theMRGCMosaicTTFresponsesAllConditions', ...
        'theMRGCMosaicTemporalSupportSecondsAllConditions');

    % Assert that the sources are configured correctly
    assert(strcmp(stimParams.stimulusShape, TTFparamsStruct.stimulusShape), ...
        'Stimulus shapes do not agree in ''stimulusParams'' and ''TTFparamsStruct''.');

    assert( ...
        (contains(targetCellImpulseResponseSource, 'center') && (strcmp(stimParams.stimulusShape, 'spot'))) || ...
        (contains(targetCellImpulseResponseSource, 'surround') && (strcmp(stimParams.stimulusShape, 'annulus'))), ...
         'Stimulus shape (''%s'') does not agree with the source of the targetCellImpulse response (''%s'').', ...
         stimParams.stimulusShape, targetCellImpulseResponseSource)


    temporalFrequenciesExamined = TTFparamsStruct.tfSupport;
    temporalFrequenciesExamined = temporalFrequenciesExamined(temporalFrequenciesExamined< 120);


    theTTFamplitude = zeros(1, numel(temporalFrequenciesExamined));
    theTTFphaseDegs = zeros(1, numel(temporalFrequenciesExamined));

    for iTF = 1:numel(temporalFrequenciesExamined)

        % Retrieve the photocurrent response data for this TF
        theTemporalSupportSecondsForThisTF = squeeze(theMRGCMosaicTemporalSupportSecondsAllConditions(iTF,:));
        theMRGCMosaicPhotocurrentResponsesForThisTF = squeeze(theMRGCMosaicTTFresponsesAllConditions(iTF,:,computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex));

        nanIndices = find(isnan(squeeze(theMRGCMosaicTTFresponsesAllConditions(iTF,:,computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex))));
        if (isempty(nanIndices))
            theTimeBins = 1:size(theMRGCMosaicTTFresponsesAllConditions,2);
        else
            theTimeBins = 1:(nanIndices(1)-1);
        end

        theTemporalSupportSecondsForThisTF = theTemporalSupportSecondsForThisTF(theTimeBins);
        theMRGCMosaicPhotocurrentResponsesForThisTF = theMRGCMosaicPhotocurrentResponsesForThisTF(theTimeBins);

        assert(isempty(find(isnan(theMRGCMosaicPhotocurrentResponsesForThisTF(:)))), 'did not remove all nan part of the response');
        assert(isempty(find(isnan(theTemporalSupportSecondsForThisTF(:)))), 'did not remove all nan part of the response temporal support');

        dt = theTemporalSupportSecondsForThisTF(2)-theTemporalSupportSecondsForThisTF(1);
        temporalSupport1Msec = 0:1/1000:(theTemporalSupportSecondsForThisTF(end)+dt);

        [theFittedResponse, fittedParams] = RGCMosaicConstructor.helper.fit.sinusoidToResponseTimeSeries(...
            theTemporalSupportSecondsForThisTF, ...
            theMRGCMosaicPhotocurrentResponsesForThisTF, ...
            temporalFrequenciesExamined(iTF), ...
            temporalSupport1Msec, ...
            'allowOffset', allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries);

        theTTFamplitude(iTF) = fittedParams(1);
        theTTFphaseDegs(iTF) = fittedParams(2);

        if (visualizeSinusoidalFits)
            hFig = figure(1); clf;
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
            theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
            ax = theAxes{1,1};

            if (theTemporalSupportSecondsForThisTF(end) > 1)
                dTick = 0.4;
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.5)
                dTick = 0.2;
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.3)
                dTick = 0.1;
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.1)
                dTick = 0.05;
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.06)
                dTick = 0.02;
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.03)
                dTick = 0.01;
            else
                dTick = 0.005;
            end
            xTicks = 0:dTick:2.0;
            if (theTemporalSupportSecondsForThisTF(end) > 0.5)
                xTickLabels = sprintf('%.1f\n', xTicks);
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.1)
                xTickLabels = sprintf('%.2f\n', xTicks);
            else
                xTickLabels = sprintf('%.3f\n', xTicks);
            end

            XLims = [0 theTemporalSupportSecondsForThisTF(end)];
            YLims = [-2 2];

            plot(ax,theTemporalSupportSecondsForThisTF, theMRGCMosaicPhotocurrentResponsesForThisTF, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
            hold(ax, 'on')
            plot(ax, temporalSupport1Msec, theFittedResponse, 'k-', 'LineWidth', 1.5);
            set(ax, 'XLim', XLims, 'YLim', YLims, 'XTick', xTicks, 'YTick', [-20:1:2], 'XTickLabel', xTickLabels);
            grid(ax, 'on');
            xlabel(ax,'time (sec)');
            ylabel(ax,'pAmps');
            title(ax, sprintf('%s, TF = %2.2fHz', stimulusShape, temporalFrequenciesExamined(iTF)));

            % Finalize figure using the Publication-Ready format
            PublicationReadyPlotLib.applyFormat(ax,ff);
            PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

            % Export
            visualizationPDFfileName = sprintf('%s_TF_%2.2fHz', stimulusShape, temporalFrequenciesExamined(iTF));
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
    thePhotocurrentBasedMRGCcellTTF = theTTFamplitude .* exp(-1i * (theTTFphaseRadians));


    % Adjust TTF for missing 0 Hz point
    minimumDelaySecondsForEstimationOfBaseline = 0.8;
    [thePhotocurrentBasedMRGCcellTTF, temporalFrequenciesExamined] = RGCMosaicConstructor.temporalFilterEngine.adjustTTFtoDealWithMissingTTFsampleAt0Hz(...
        thePhotocurrentBasedMRGCcellTTF, temporalFrequenciesExamined, minimumDelaySecondsForEstimationOfBaseline);

    % Add a delay of 0 msec to make the photocurrent impulse response (as esimated from its TTF) is causal
    delaySeconds = 0/1000;
    omega = 2 * pi * temporalFrequenciesExamined;
    thePhotocurrentBasedMRGCcellTTF = exp(-1i * omega * delaySeconds) .* thePhotocurrentBasedMRGCcellTTF;


    verifyOffsetCorrection = true;
    if (verifyOffsetCorrection)
        % Verify the the impulse response has a baseline very close to 0
        thePhotocurrentBasedImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                    thePhotocurrentBasedMRGCcellTTF , temporalFrequenciesExamined);
    
        figure(1);
        plot(thePhotocurrentBasedImpulseResponseData.temporalSupportSeconds, thePhotocurrentBasedImpulseResponseData.amplitude, 'ko');
        hold on;
        plot(thePhotocurrentBasedImpulseResponseData.temporalSupportSeconds, thePhotocurrentBasedImpulseResponseData.temporalSupportSeconds*0, 'r-');
        disp('Is baseline near 0? Hit enter to proceed') 
        pause
    end



    % Compute theTargetCascadedFilterTTF 
    switch (targetCellImpulseResponseSource)
        case 'Benardete&Kaplan 1997, Figure 6 (ON), center'
            params = RGCmodels.BenardeteKaplan1997.figure6CenterSurroundFilterParams('ON');
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.centerIR.pVector, temporalFrequenciesExamined);
    
        case 'Benardete&Kaplan 1997, Figure 6 (OFF), center'
            params = RGCmodels.BenardeteKaplan1997.figure6CenterSurroundFilterParams('OFF');
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.centerIR.pVector, temporalFrequenciesExamined);

        case 'Benardete&Kaplan 1997, Figure 7, center'
            params = RGCmodels.BenardeteKaplan1997.figure7CenterSurroundFilterParams();
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.centerIR.pVector, temporalFrequenciesExamined);

        case 'Benardete&Kaplan 1997, Figure 6 (ON), surround'
            params = RGCmodels.BenardeteKaplan1997.figure6CenterSurroundFilterParams('ON');
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.surroundIR.pVector, temporalFrequenciesExamined);


        case 'Benardete&Kaplan 1997, Figure 6 (OFF), surround'
            params = RGCmodels.BenardeteKaplan1997.figure6CenterSurroundFilterParams('OFF');
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.surroundIR.pVector, temporalFrequenciesExamined);

        case 'Benardete&Kaplan 1997, Figure 7, surround'
            params = RGCmodels.BenardeteKaplan1997.figure7CenterSurroundFilterParams();
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.surroundIR.pVector, temporalFrequenciesExamined);

        otherwise
            error('Unknown source for target impulse response: ''%s''.', targetCellImpulseResponseSource);
    end



    % Derive theInnerRetinaTTF based on the theTargetCascadedFilterTTF and thePhotocurrentBasedMRGCcellTTF
    switch (temporalFilterSynthesisMethod)
        case 'direct division of TTFs'
            theInnerRetinaTTF = theTargetCascadedFilterTTF ./ thePhotocurrentBasedMRGCcellTTF;

        case 'delay + highpass filter'
           
            % Filter to fit
            filterType =  'differenceOfLowPassFilters'; % 'dampedOscillationFiter'; %'delayLeadLagFilter'; %'delayHighPassFilter'; % 'asymmetricBandPassFilter'; %'delayLeadLagFilter';
            
            % Frequency weighting
            minFrequencyHzWithUnitWeight = 3.0;
            frequencyWeights = linspace(0,1, numel(temporalFrequenciesExamined));
            frequencyWeights(temporalFrequenciesExamined>=minFrequencyHzWithUnitWeight) = frequencyWeights(temporalFrequenciesExamined==minFrequencyHzWithUnitWeight);
            frequencyWeights = frequencyWeights/max(frequencyWeights);

            solverType = 'multi-start'; % 'fmincon'; %'global-search' ; % 'fmincon'; %'multi-start'; 'fmincon'; %'multi-start'; %'global-search';
            multiStartsNum = 32;
            useParallel = ~true;

            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.innerRetinaFilterBetweenPhotocurrentBasedAndTargetTTF(...
                temporalFrequenciesExamined, theTargetCascadedFilterTTF, thePhotocurrentBasedMRGCcellTTF, ...
                frequencyWeights, filterType, solverType, multiStartsNum, useParallel);

        otherwise
            error('Unknown temporal filter synthesis method: ''%s''.', temporalFilterSynthesisMethod);
    end


   

    % The photocurrents based mRGC impulse response 
    thePhotocurrentBasedMRGCcellImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        thePhotocurrentBasedMRGCcellTTF, temporalFrequenciesExamined);

    % Compute the achieved cascaded filter TTF
    theAchievedCascadedFilterMRGCcellTTF = theInnerRetinaTTF .* thePhotocurrentBasedMRGCcellTTF;

    % The achieved cascaded filter mRGC impulse response 
    theAchievedCascadedFilterMRGCcellImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        theAchievedCascadedFilterMRGCcellTTF, temporalFrequenciesExamined);

    scalingFactor = max(abs(thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude))/max(abs(theAchievedCascadedFilterMRGCcellImpulseResponseData .amplitude))
    
    % Scale the innerRetinaTTF by it
    theInnerRetinaTTF = theInnerRetinaTTF * scalingFactor;
   

    % Compute theInnerRetinaImpulseResponse
    theInnerRetinaImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        theInnerRetinaTTF, temporalFrequenciesExamined);

    % Compute the achieved cascaded filter TTF
    theAchievedCascadedFilterTTF = theInnerRetinaTTF .* thePhotocurrentBasedMRGCcellTTF;

    % Compute the achieved impulse response
    theAchievedCascadedFilterImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        theAchievedCascadedFilterTTF, temporalFrequenciesExamined);


    % Compute theTargetCascadedFiltermpulseResponse
    theTargetCascadedFilterImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        theTargetCascadedFilterTTF, temporalFrequenciesExamined);


   
    theMacaqueImpulseResponseData = [];
    theMacaqueTTFData = [];

    switch (targetCellImpulseResponseSource)
        case 'Benardete&Kaplan 1997, Figure 6 (ON), center'
            theMacaqueImpulseResponseData = RGCmodels.BenardeteKaplan1997.digitizedData.ONcenterDiskImpulseResponseFromFigure6;

        case 'Benardete&Kaplan 1997, Figure 6 (ON), surround'
            theMacaqueImpulseResponseData = RGCmodels.BenardeteKaplan1997.digitizedData.ONcenterAnnulusImpulseResponseFromFigure6;
    end

    % Plot the impulse response functions
    % Plot the TTF
    hFig = figure(60); clf;
    
    plotComboData(hFig, ...
        temporalFrequenciesExamined, ...
        thePhotocurrentBasedMRGCcellTTF/max(abs(thePhotocurrentBasedMRGCcellTTF)), ...
        theTargetCascadedFilterTTF / max(abs(theTargetCascadedFilterTTF)), ...
        theAchievedCascadedFilterTTF/max(abs(theAchievedCascadedFilterTTF)), ...
        theInnerRetinaTTF/max(abs(theInnerRetinaTTF)), ...
        targetCellImpulseResponseSource, ...
        [0 1.01], [0.3 200], [0.3 1 3 10 30 100], 'log', 'frequency (Hz)', ...
        'TTFs', ...
        'withMacaqueData', theMacaqueTTFData);


    hFig = figure(61); clf;
    plotComboData(hFig, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.temporalSupportSeconds*1e3, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude, ...
        theTargetCascadedFilterImpulseResponseData.amplitude, ...
        theAchievedCascadedFilterImpulseResponseData.amplitude, ...
        theInnerRetinaImpulseResponseData.amplitude, ...
        targetCellImpulseResponseSource, ...
        1.01*[-1 1], [0 200], 0:20:200, 'linear', 'time (msec)', ...
        'IRFs', ...
        'withMacaqueData', theMacaqueImpulseResponseData);


    
end

function plotComboData(hFig, ...
    theAxisData, ...
    thePhotocurrentBasedMRGCcellData, ...
    theTargetCascadedFilterData, ...
    theAchievedCascadedFilterData, ...
    theInnerRetinaData, ...
    targetCellImpulseResponseSource, ...
    yAxisLims, xAxisLims, xAxisTicks, xAxisScale, xAxisLabel, ...
    pdfFileNamePostFix, varargin)

    p = inputParser;
    p.addParameter('withMacaqueData', [], @(x)(isempty(x) || (isstruct(x))));
    p.parse(varargin{:});
    theMacaqueData = p.Results.withMacaqueData;
    

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard very wide figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    allPlotHandles = [];
    allLegends = {};

    hold(ax, 'on');

    % The photocurrents-based Data
    if (strcmp(pdfFileNamePostFix, 'TTFs'))
        theAmplitude = abs(thePhotocurrentBasedMRGCcellData);
    else
        theAmplitude = thePhotocurrentBasedMRGCcellData;
    end
    maxPhotocurrentsBasedAmplitude = max(abs(theAmplitude));

    p1 = plot(ax, theAxisData, theAmplitude/maxPhotocurrentsBasedAmplitude, 'ro-', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
    allPlotHandles(numel(allPlotHandles)+1) = p1;
    allLegends{numel(allLegends)+1} = 'photocurrents-derived';

    if (~isempty(theTargetCascadedFilterData))
        % The target cascaded filter TTF
        if (strcmp(pdfFileNamePostFix, 'TTFs'))
            theAmplitude = abs(theTargetCascadedFilterData);
        else
            theAmplitude = theTargetCascadedFilterData;
        end
        maxTargetAmplitude = max(abs(theAmplitude));
    
        plot(ax, theAxisData, theAmplitude/maxTargetAmplitude, 'k-', ...
                 'LineWidth', 4, ...
                 'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
        p2 = plot(ax,theAxisData, theAmplitude/maxTargetAmplitude, 'c-', ...
                 'LineWidth', 2, ...
                 'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
        allPlotHandles(numel(allPlotHandles)+1) = p2;
        allLegends{numel(allLegends)+1} = sprintf('%s (fitted model)', targetCellImpulseResponseSource);
    else
        maxTargetAmplitude = [];
    end



    % The achieved cascaded filter TTF
    if (strcmp(pdfFileNamePostFix, 'TTFs'))
        theAmplitude = abs(theAchievedCascadedFilterData);
    else
        theAmplitude = theAchievedCascadedFilterData;
    end

    plot(ax, theAxisData, theAmplitude/maxPhotocurrentsBasedAmplitude, 'b-', ...
             'LineWidth', 4)
    p5 = plot(ax,theAxisData, theAmplitude/maxPhotocurrentsBasedAmplitude, 'r-', ...
             'LineWidth', 3, ...
             'Color', [0.0 1 1]);

    allPlotHandles(numel(allPlotHandles)+1) = p5;
    allLegends{numel(allLegends)+1} = 'ISETBio model (photocurrent x inner retina filter cascade)';


    % The derived inner retina filter TTF
    if (strcmp(pdfFileNamePostFix, 'TTFs'))
        theAmplitude = abs(theInnerRetinaData);
    else
        theAmplitude = theInnerRetinaData;
    end
    maxAmplitude = max(abs(theAmplitude));

    plot(ax,theAxisData, theAmplitude/maxAmplitude, 'k-', ...
              'LineWidth', 4, ...
              'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
    p4 = plot(ax,theAxisData, theAmplitude/maxAmplitude, 'y-', ...
              'LineWidth', 2, ...
              'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);

    allPlotHandles(numel(allPlotHandles)+1) = p4;
    allLegends{numel(allLegends)+1} = 'inner retina filter';


    if (~isempty(theMacaqueData))
        if (strcmp(pdfFileNamePostFix, 'TTFs'))
            theAxisData = [];
            theAmplitude = [];
        else
            theAxisData = theMacaqueData.temporalSupportSeconds * 1e3;
            theAmplitude = theMacaqueData.amplitude;
        end
        if (isempty(maxTargetAmplitude))
            maxTargetAmplitude = max(abs(theAmplitude));
        end

        p7 = plot(ax,theAxisData, theAmplitude/maxTargetAmplitude, 'bv', ...
                  'LineWidth', 2, ...
                  'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.9]);
    
        allPlotHandles(numel(allPlotHandles)+1) = p7;
        allLegends{numel(allLegends)+1} = sprintf('%s (raw data)', targetCellImpulseResponseSource);
    end


    legend(ax, allPlotHandles, allLegends, 'Location', 'NorthEast', 'NumColumns', 1);
    xlabel(ax, xAxisLabel)
    set(ax, 'XLim', xAxisLims, 'XTick', xAxisTicks, 'YLim', yAxisLims, 'XScale', xAxisScale);
    grid(ax, 'on')
    ff.legendBox = 'on';

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    %PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

    % Export
    visualizationPDFfileName = sprintf('%s_%s', targetCellImpulseResponseSource, pdfFileNamePostFix);
    exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
    pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

    thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');
 
end
