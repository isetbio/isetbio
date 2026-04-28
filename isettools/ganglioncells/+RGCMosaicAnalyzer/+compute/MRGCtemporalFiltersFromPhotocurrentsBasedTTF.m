%
%RGCMosaicAnalyzer.compute.MRGCtemporalFiltersFromPhotocurrentsBasedTTF
%
%
function MRGCtemporalFiltersFromPhotocurrentsBasedTTF(...
    innerRetinaFilterDerivationParams, ...
    targetCellImpulseResponseSource, ...
    stimulusShape, ...
    allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries, ...
    theTargetRGCindex, ...
    theMRGCMosaicTTFResponsesFullFileName, ...
    theAnalyzedTTFsFullFileName, ...
    visualizeSinusoidalFits)


    % Load the measured TTF responses
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

    % Assert that the specified RGCindex matches computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex
    assert(theTargetRGCindex == computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex, ...
        sprintf('Different targetRGCindex specified (%d)than what the one for which the TTFresponses were computed for (%d).\nThis will result in a TTF computed a stimulus not exactly centered on RGC %d\n', ...
        theTargetRGCindex, computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex, theTargetRGCindex));
   
    % Compute the complex temporal transfer function (TTF) theTargetRGCindex
    minimumDelaySecondsForEstimationOfBaseline = 800/1000;

    % Compute thePhotocurrentBasedMRGCcellTTF, for a zero-baseline,
    % centered & Tukey windowed impulse response
    [thePhotocurrentBasedMRGCcellTTF, temporalFrequencySupportHz] = computeComplexValuedTTFforSingleMRGC(...
        TTFparamsStruct, ...
        theMRGCMosaicTemporalSupportSecondsAllConditions, ...
        theMRGCMosaicTTFresponsesAllConditions, ...
        theTargetRGCindex, ...
        allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries, ...
        visualizeSinusoidalFits);


    % Compute the TTF after adjusting for a zero baseline IR
    verifyOffsetCorrection = ~true;
    centerIR = true;
    [thePhotocurrentBasedMRGCcellTTF, delaySecondsForPhotocurrentIR] = RGCMosaicConstructor.temporalFilterEngine.zeroBaselineWindowedCenteredTTF(...
        thePhotocurrentBasedMRGCcellTTF, temporalFrequencySupportHz, minimumDelaySecondsForEstimationOfBaseline, centerIR, verifyOffsetCorrection);

    % Compute theTargetCascadedFilterTTF from the Benardete&Kaplan work
    switch (targetCellImpulseResponseSource)
        case 'Benardete&Kaplan 1997, Figure 6 (ON), center'
            params = RGCmodels.BenardeteKaplan1997.figure6CenterSurroundFilterParams('ON');
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.centerIR.pVector, temporalFrequencySupportHz);
    
        case 'Benardete&Kaplan 1997, Figure 6 (OFF), center'
            params = RGCmodels.BenardeteKaplan1997.figure6CenterSurroundFilterParams('OFF');
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.centerIR.pVector, temporalFrequencySupportHz);

        case 'Benardete&Kaplan 1997, Figure 7, center'
            params = RGCmodels.BenardeteKaplan1997.figure7CenterSurroundFilterParams();
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.centerIR.pVector, temporalFrequencySupportHz);

        case 'Benardete&Kaplan 1997, Figure 6 (ON), surround'
            params = RGCmodels.BenardeteKaplan1997.figure6CenterSurroundFilterParams('ON');
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.surroundIR.pVector, temporalFrequencySupportHz);


        case 'Benardete&Kaplan 1997, Figure 6 (OFF), surround'
            params = RGCmodels.BenardeteKaplan1997.figure6CenterSurroundFilterParams('OFF');
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.surroundIR.pVector, temporalFrequencySupportHz);

        case 'Benardete&Kaplan 1997, Figure 7, surround'
            params = RGCmodels.BenardeteKaplan1997.figure7CenterSurroundFilterParams();
            theTargetCascadedFilterTTF = RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF(...
                params.surroundIR.pVector, temporalFrequencySupportHz);

        otherwise
            error('Unknown source for target impulse response: ''%s''.', targetCellImpulseResponseSource);
    end


    % Compute the TTF after adjusting for a zero baseline IR
    verifyOffsetCorrection = ~true;
    centerIR = false;
    [theTargetCascadedFilterTTF, delaySecondsForTargetIR] = RGCMosaicConstructor.temporalFilterEngine.zeroBaselineWindowedCenteredTTF(...
        theTargetCascadedFilterTTF, temporalFrequencySupportHz, minimumDelaySecondsForEstimationOfBaseline, centerIR, verifyOffsetCorrection);

    % Normalize the TTFs
    thePhotocurrentBasedMRGCcellTTF = thePhotocurrentBasedMRGCcellTTF / max(abs(thePhotocurrentBasedMRGCcellTTF(:)));
    theTargetCascadedFilterTTF = theTargetCascadedFilterTTF / max(abs(theTargetCascadedFilterTTF(:)));


    % ONLY for debuggin purposes
    verifyOffsetCorrection = ~true;
    if (verifyOffsetCorrection)

        theTargetCascadedFilterImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                    theTargetCascadedFilterTTF, temporalFrequencySupportHz);


        thePhotocurrentBasedImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                    thePhotocurrentBasedMRGCcellTTF, temporalFrequencySupportHz);


        assert(numel(theTargetCascadedFilterImpulseResponseData.temporalSupportSeconds) == numel(thePhotocurrentBasedImpulseResponseData.temporalSupportSeconds), ...
            'unequal temporal support lengths');

        figure(1);
        p1 = plot(thePhotocurrentBasedImpulseResponseData.temporalSupportSeconds, thePhotocurrentBasedImpulseResponseData.amplitude/max(thePhotocurrentBasedImpulseResponseData.amplitude), 'r-');
        hold on;
        p2 = plot(theTargetCascadedFilterImpulseResponseData.temporalSupportSeconds, theTargetCascadedFilterImpulseResponseData.amplitude/max(thePhotocurrentBasedImpulseResponseData.amplitude), 'b-');
        plot(thePhotocurrentBasedImpulseResponseData.temporalSupportSeconds, thePhotocurrentBasedImpulseResponseData.temporalSupportSeconds*0, 'k-');
        legend([p1 p2], {'photocurrents-based IR', 'target IR'})
        set(gca, 'FontSize', 20)
        pause
    end



    

    % Derive theInnerRetinaTTF based on the theTargetCascadedFilterTTF and thePhotocurrentBasedMRGCcellTTF
    switch (innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod)
        case 'direct division of TTFs'

            idx = find(abs(thePhotocurrentBasedMRGCcellTTF)>10*eps);
            theInnerRetinaTTF = theTargetCascadedFilterTTF*0;
            theInnerRetinaTTF(idx) = theTargetCascadedFilterTTF(idx)./thePhotocurrentBasedMRGCcellTTF(idx);

            % Undo the delays we introduced to center the photocurrent and target TTFs
            delaySeconds = -delaySecondsForPhotocurrentIR;
            omega = 2 * pi * temporalFrequencySupportHz;
            thePhotocurrentBasedMRGCcellTTF = exp(-1i * omega * delaySeconds) .* thePhotocurrentBasedMRGCcellTTF;
        
            delaySeconds = -delaySecondsForTargetIR;
            theTargetCascadedFilterTTF  = exp(-1i * omega * delaySeconds) .* theTargetCascadedFilterTTF;
        
            % Undo the corresponding delay in the inner retinal TTF
            delaySeconds = -(delaySecondsForTargetIR-delaySecondsForPhotocurrentIR);
            theInnerRetinaTTF = exp(-1i * omega * delaySeconds) .* theInnerRetinaTTF;


            modelParams = [];
            frequencyWeights = temporalFrequencySupportHz*0+1;
            innerRetinaFilterDerivationParams.minFrequencyHzWithNonZeroWeight = temporalFrequencySupportHz(1);
            innerRetinaFilterDerivationParams.maxFrequencyHzWithNonZeroWeight = temporalFrequencySupportHz(end);

        case {'differenceOfLowPassFilters', 'differenceOfLowPassFilters2', 'sumOfLowPassFilters', 'dampedOscillationFilter', 'dampedOscillationLowPassCascadeFilter', 'delayLeadLagFilter', 'delayLeadLagFilter2', 'delayHighPassFilter', 'asymmetricBandPassFilter'}
           
            idx = find(abs(thePhotocurrentBasedMRGCcellTTF)>10*eps);
            theIdealInnerRetinaTTF = theTargetCascadedFilterTTF*0;
            theIdealInnerRetinaTTF(idx) = theTargetCascadedFilterTTF(idx)./thePhotocurrentBasedMRGCcellTTF(idx);

            % Undo the delays we introduced to center the photocurrent and target TTFs
            delaySeconds = -(delaySecondsForTargetIR-delaySecondsForPhotocurrentIR);
            omega = 2 * pi * temporalFrequencySupportHz;
            theIdealInnerRetinaTTF = exp(-1i * omega * delaySeconds) .* theIdealInnerRetinaTTF;

            % Undo the delays we introduced to center the photocurrent and target TTFs
            delaySeconds = -delaySecondsForPhotocurrentIR;
            omega = 2 * pi * temporalFrequencySupportHz;
            thePhotocurrentBasedMRGCcellTTF = exp(-1i * omega * delaySeconds) .* thePhotocurrentBasedMRGCcellTTF;
        

            % Frequency weighting
            [~,bin1] = min(abs(temporalFrequencySupportHz - innerRetinaFilterDerivationParams.minFrequencyHzWithNonZeroWeight));
            [~,bin2] = min(abs(temporalFrequencySupportHz - innerRetinaFilterDerivationParams.maxFrequencyHzWithNonZeroWeight));

            % Weights for residuals in frequency domain
            frequencyWeights = temporalFrequencySupportHz*0+1;
            frequencyWeights(1:bin1) = 0.0;
            frequencyWeights(bin2:end) = 0.0;
            frequencyWeights(bin1:bin2) = 1; % (1-((0:(bin2-bin1))*(1/(bin2-bin1)))).^1.5;
           
            % Time range for comptuting temporal residuals 
            temporalWeightingLimitsSeconds = [innerRetinaFilterDerivationParams.minTimeDelaySecondsWithUnitWeight innerRetinaFilterDerivationParams.maxTimeDelaySecondsWithUnitWeight];
 

           

            [theInnerRetinaTTF, modelParams] = RGCMosaicConstructor.temporalFilterEngine.deriveInnerRetinaFilterBetweenPhotocurrentBasedAndTargetTTF(...
                temporalFrequencySupportHz, theTargetCascadedFilterTTF, thePhotocurrentBasedMRGCcellTTF, theIdealInnerRetinaTTF, ...
                minimumDelaySecondsForEstimationOfBaseline, ...
                frequencyWeights, temporalWeightingLimitsSeconds, ...
                innerRetinaFilterDerivationParams.timeDomainResidualWeighting, ...
                innerRetinaFilterDerivationParams.residualIsBasedOnTTFofCascadedPhotocurrentInnerRetinaFilter, ...
                innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod, ...
                innerRetinaFilterDerivationParams.amplitudeSpectrumVsComplexSpectrumBias, ...
                innerRetinaFilterDerivationParams.solverType, ...
                innerRetinaFilterDerivationParams.multiStartsNum, ...
                innerRetinaFilterDerivationParams.useParallel);



        otherwise
            error('Unknown temporal filter synthesis method: ''%s''.', temporalFilterSynthesisMethod);
    end

    


    fprintf('Saving derived inner retina filter TTF to %s\n', theAnalyzedTTFsFullFileName);

    innerRetinaFilterDataStruct = struct(...
        'temporalFilterSynthesisMethod', innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod, ...
        'targetCellImpulseResponseSource', targetCellImpulseResponseSource, ...
        'temporalFrequencySupportHz', temporalFrequencySupportHz, ...
        'targetMacaqueTTF', theTargetCascadedFilterTTF, ...
        'achievedTargetMRGCcellTTF', thePhotocurrentBasedMRGCcellTTF .* theInnerRetinaTTF, ...
        'photocurrentBasedTTF', thePhotocurrentBasedMRGCcellTTF, ...
        'derivedInnerRetinaTTF', theInnerRetinaTTF, ...
        'fittedModelParams', modelParams, ...
        'frequencyWeights', frequencyWeights);

    save(theAnalyzedTTFsFullFileName, ...
        'theMRGCMosaic', 'stimParams', 'TTFparamsStruct', ...
        'innerRetinaFilterDerivationParams', ...
        'innerRetinaFilterDataStruct');
   

    hFig = figure(9876); clf;
    set(hFig, 'Position', [10 10 2000 800], 'Name', sprintf('%s - %s', targetCellImpulseResponseSource, stimulusShape))

    % The amplitude spectra of the target and the achieved TTFs
    ax = subplot('Position', [0.05 0.05 0.25 0.9]);
    p1 = RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
        ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.targetMacaqueTTF, 'o', ...
        true, false, [0 0 0], ...
        '');
    hold(ax, 'on');
    p2 = RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
        ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF, '-', ...
        true, false, [1 0 0], ...
        '');
    plot(ax, innerRetinaFilterDerivationParams.minFrequencyHzWithNonZeroWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);
    plot(ax, innerRetinaFilterDerivationParams.maxFrequencyHzWithNonZeroWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);
    legend(ax, [p1 p2], {'target' 'achieved'});


    % The phase spectra of the target and the achieved TTFs
    ax = subplot('Position', [0.4 0.05 0.25 0.9]);
    p1 = RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
        ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.targetMacaqueTTF, 'o', ...
        true, false, [0 0 0], ...
        '');
    hold(ax, 'on');
    p2 = RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
        ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF, '-', ...
        true, false, [1 0 0], ...
        '');
    plot(ax, innerRetinaFilterDerivationParams.minFrequencyHzWithNonZeroWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);
    plot(ax, innerRetinaFilterDerivationParams.maxFrequencyHzWithNonZeroWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);
    legend(ax, [p1 p2], {'target' 'achieved'});


    
    % The IRs of the target and the achieved TTFs
    theTargetMacaqueIR = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        innerRetinaFilterDataStruct.targetMacaqueTTF, ...
        innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        'causal', false); 

    % The achieved cascaded filter mRGC impulse response 
    theAchievedTargetMRGCcellImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF, ...
        innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        'causal', false); 


    % The photocurrents based mRGC impulse response 
    thePhotocurrentBasedMRGCcellImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        innerRetinaFilterDataStruct.photocurrentBasedTTF, temporalFrequencySupportHz, ...
        'causal', false); 

    % The derived inner retina impulse response
    theDerivedInnerRetinaImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTF, temporalFrequencySupportHz, ...
        'causal', false); 

    ax = subplot('Position', [0.7 0.05 0.25 0.9]);
    p1 = plot(ax, theTargetMacaqueIR.temporalSupportSeconds*1e3, theTargetMacaqueIR.amplitude, 'o');
    hold(ax, 'on');
    p2 = plot(ax, theAchievedTargetMRGCcellImpulseResponseData.temporalSupportSeconds*1e3, theAchievedTargetMRGCcellImpulseResponseData.amplitude, 'r-');
    set(ax, 'XLim', [0 300]);
    legend(ax, [p1 p2], {'target' 'achieved'});
    pause


   
    theMacaqueImpulseResponseData = [];
    theMacaqueTTFData = [];

    switch (targetCellImpulseResponseSource)
        case 'Benardete&Kaplan 1997, Figure 6 (ON), center'
            theMacaqueImpulseResponseData = RGCmodels.BenardeteKaplan1997.digitizedData.ONcenterDiskImpulseResponseFromFigure6;
            theMacaqueImpulseResponseData.amplitude = theMacaqueImpulseResponseData.amplitude/max(abs(theMacaqueImpulseResponseData.amplitude));

        case 'Benardete&Kaplan 1997, Figure 6 (ON), surround'
            theMacaqueImpulseResponseData = RGCmodels.BenardeteKaplan1997.digitizedData.ONcenterAnnulusImpulseResponseFromFigure6;
            theMacaqueImpulseResponseData.amplitude = theMacaqueImpulseResponseData.amplitude/max(abs(theMacaqueImpulseResponseData.amplitude));
    end

    % Plot the impulse response functions
    % Plot the TTF
    hFig = figure(60); clf;
    
    % The spectra
    plotComboData(hFig, ...
        temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.photocurrentBasedTTF/max(abs(innerRetinaFilterDataStruct.photocurrentBasedTTF)), ...
        innerRetinaFilterDataStruct.targetMacaqueTTF / max(abs(innerRetinaFilterDataStruct.targetMacaqueTTF)), ...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF/max(abs(innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF)), ...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTF/max(abs(innerRetinaFilterDataStruct.derivedInnerRetinaTTF)), ...
        targetCellImpulseResponseSource, ...
        [0 1.01], [0.3 200], [0.3 1 3 10 30 100], 'log', 'frequency (Hz)', ...
        'TTFs', ...
        'withMacaqueData', theMacaqueTTFData);


    hFig = figure(61); clf;
    plotComboData(hFig, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.temporalSupportSeconds*1e3, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude/max(abs(thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude)), ...
        theTargetMacaqueIR.amplitude/max(abs(theTargetMacaqueIR.amplitude)), ...
        theAchievedTargetMRGCcellImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
        theDerivedInnerRetinaImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
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

    p5 = plot(ax, theAxisData, theAmplitude/maxPhotocurrentsBasedAmplitude, 'k--', ...
             'LineWidth', 4);

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



function [theComplexTTF, temporalFrequencySupportHz] = computeComplexValuedTTFforSingleMRGC(...
        TTFparamsStruct, ...
        theMRGCMosaicTemporalSupportSecondsAllConditions, ...
        theMRGCMosaicTTFresponsesAllConditions, ...
        theRGCindex, ...
        allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries, ...
        visualizeSinusoidalFits)

    temporalFrequenciesExamined = TTFparamsStruct.tfSupport;

    theTTFamplitude = zeros(1, numel(temporalFrequenciesExamined));
    theTTFphaseDegs = zeros(1, numel(temporalFrequenciesExamined));

    for iTF = 1:numel(temporalFrequenciesExamined)

        % Retrieve the photocurrent response data for this TF
        theTemporalSupportSecondsForThisTF = squeeze(theMRGCMosaicTemporalSupportSecondsAllConditions(iTF,:));
        theMRGCMosaicPhotocurrentResponsesForThisTF = squeeze(theMRGCMosaicTTFresponsesAllConditions(iTF,:,theRGCindex));

        % find valid time bins for this TF (we have stored 1 period of the
        % respons, so each TF response has a different length)
        nanIndices = find(isnan(theMRGCMosaicPhotocurrentResponsesForThisTF));
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
    theComplexTTF = theTTFamplitude .* exp(-1i * (theTTFphaseRadians));

    % Adjust TTF for missing 0 Hz point.
    [theComplexTTF, temporalFrequencySupportHz] = RGCMosaicConstructor.temporalFilterEngine.adjustTTFtoDealWithMissingTTFsampleAt0Hz(...
        theComplexTTF, temporalFrequenciesExamined);

    % Clean up the TTF by removing phase outliers
    theComplexTTF = RGCMosaicConstructor.temporalFilterEngine.removePhaseOutliersFromTTF(...
        temporalFrequencySupportHz, theComplexTTF);


end
