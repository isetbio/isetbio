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
    delayMsecondsForPhaseComputation, ...
    visualizeSinusoidalFits, ...
    onlyVisualizePreviouslySynthesizedFilters)


     % Derive theInnerRetinaTTF based on the theTargetCascadedFilterTTF and thePhotocurrentBasedMRGCcellTTF
    switch (innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod)
        case 'direct division of TTFs'

            theAnalyzedTTFsFullFileName = strrep(theAnalyzedTTFsFullFileName, '???', 'direct');

        case {'differenceOfLowPassFilters', 'dampedOscillationFilter', 'dampedOscillationLowPassCascadeFilter', 'delayLeadLagFilter', 'delayHighPassFilter'}

            % Add the filter synthesis method name to theAnalyzedTTFsFullFileName
            switch (innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod)
                 case 'dampedOscillationFilter'
                     shortTemporalFilterSynthesisMethodName = 'dOsc';
                 case 'dampedOscillationLowPassCascadeFilter'
                     shortTemporalFilterSynthesisMethodName = 'dOscLP';
                 case 'differenceOfLowPassFilters'
                     shortTemporalFilterSynthesisMethodName = 'diffLP';
                 case 'delayLeadLagFilter'
                     shortTemporalFilterSynthesisMethodName = 'ld-lag';
                 case 'delayHighPassFilter'
                     shortTemporalFilterSynthesisMethodName = 'hp-lp';
                 otherwise
                        error('No short name for  synthesis method: ''%s''.', temporalFilterSynthesisMethod);
             end

             theAnalyzedTTFsFullFileName = strrep(theAnalyzedTTFsFullFileName, '???', shortTemporalFilterSynthesisMethodName);

        otherwise
            error('Unknown temporal filter synthesis method: ''%s''.', innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod);
    end


    if (onlyVisualizePreviouslySynthesizedFilters)
        visualizeInnerRetinaFilterDerivationResults(theAnalyzedTTFsFullFileName);
        return;
    end

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

    % How to compute the phase of the complex TTF
    if (isempty(delayMsecondsForPhaseComputation))
        phaseComputationParams = [];
    else
        if ( ...
                ( contains(targetCellImpulseResponseSource, 'ON') && contains(targetCellImpulseResponseSource, 'surround') ) || ...
                ( contains(targetCellImpulseResponseSource, 'OFF') && contains(targetCellImpulseResponseSource, 'center')  ) ...
           )
            phaseDegsAtZeroHz = 180;
        elseif( ...
                ( contains(targetCellImpulseResponseSource, 'OFF') && contains(targetCellImpulseResponseSource, 'surround') ) || ...
                ( contains(targetCellImpulseResponseSource, 'ON') && contains(targetCellImpulseResponseSource, 'center')  ) ...
           )
            phaseDegsAtZeroHz = 0;
        else
            error('Dont know how to interpted the targer cell data for phase computation: ''%s''.', targetCellImpulseResponseSource);
        end

        phaseComputationParams = struct(...
            'delayMilliSeconds', delayMsecondsForPhaseComputation, ...
            'phaseDegsAtZeroHz', phaseDegsAtZeroHz ...
            );
    end


    if (strcmp(stimParams.stimulusShape, 'annulus'))
        theColor = [50 110 180]/255;
    else
        theColor = [240 50 80]/255;
    end
    theColor = theColor/max(theColor);
    theSinudoidalFitColor = [1 0.5 0];

    % Compute thePhotocurrentBasedMRGCcellTTF

    [thePhotocurrentBasedMRGCcellTTF, temporalFrequencySupportHz] = ...
        RGCMosaicConstructor.temporalFilterEngine.complexValuedTTFfromResponsesToSinusoidalModulations(...
            TTFparamsStruct.tfSupport, ...
            theMRGCMosaicTemporalSupportSecondsAllConditions, ...
            squeeze(theMRGCMosaicTTFresponsesAllConditions(:,:,theTargetRGCindex)), ...
            allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries, ...
            phaseComputationParams, ...
            visualizeSinusoidalFits, ...
            stimParams.stimulusShape, ...
            theColor, theSinudoidalFitColor);
    

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

        figure(1); clf;
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

            modelParams = [];
            frequencyWeights = temporalFrequencySupportHz*0+1;
            innerRetinaFilterDerivationParams.minFrequencyHzWithNonZeroWeight = temporalFrequencySupportHz(1);
            innerRetinaFilterDerivationParams.maxFrequencyHzWithNonZeroWeight = temporalFrequencySupportHz(end);


        case {'differenceOfLowPassFilters', 'dampedOscillationFilter', 'dampedOscillationLowPassCascadeFilter', 'delayLeadLagFilter', 'delayHighPassFilter'}
           
            idx = find(abs(thePhotocurrentBasedMRGCcellTTF)>10*eps);
            theIdealInnerRetinaTTF = theTargetCascadedFilterTTF*0;
            theIdealInnerRetinaTTF(idx) = theTargetCascadedFilterTTF(idx)./thePhotocurrentBasedMRGCcellTTF(idx);

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
            error('Unknown temporal filter synthesis method: ''%s''.', innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod);
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
        'theMRGCMosaic', 'theTargetRGCindex', 'targetCellImpulseResponseSource', ...
        'stimParams', 'TTFparamsStruct', ...
        'innerRetinaFilterDerivationParams', ...
        'innerRetinaFilterDataStruct');


    hFig = figure(9876); clf;
    set(hFig, 'Position', [10 10 2000 800], 'Name', sprintf('%s - %s', targetCellImpulseResponseSource, stimParams.stimulusShape))

    % The amplitude spectra of the target and the achieved TTFs
    ax = subplot('Position', [0.05 0.05 0.25 0.9]);
    p1 = RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
        ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.targetMacaqueTTF, 'o', ...
        true, false, [0 0 0], 1.0, ...
        '');
    hold(ax, 'on');
    p2 = RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
        ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF, '-', ...
        true, false, [1 0 0], 1.0, ...
        '');
    plot(ax, innerRetinaFilterDerivationParams.minFrequencyHzWithNonZeroWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);
    plot(ax, innerRetinaFilterDerivationParams.maxFrequencyHzWithNonZeroWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);
    legend(ax, [p1 p2], {'target' 'achieved'});


    % The phase spectra of the target and the achieved TTFs
    ax = subplot('Position', [0.4 0.05 0.25 0.9]);
    p1 = RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
        ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.targetMacaqueTTF, 'o', ...
        true, false, [0 0 0], 1.0, ...
        '');
    hold(ax, 'on');
    p2 = RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
        ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF, '-', ...
        true, false, [1 0 0], 1.0, ...
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
        innerRetinaFilterDataStruct.photocurrentBasedTTF, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        'causal', false); 

    % The derived inner retina impulse response
    theDerivedInnerRetinaImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTF, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        'causal', false); 

    ax = subplot('Position', [0.7 0.05 0.25 0.9]);
    p1 = plot(ax, theTargetMacaqueIR.temporalSupportSeconds*1e3, theTargetMacaqueIR.amplitude, 'o');
    hold(ax, 'on');
    p2 = plot(ax, theAchievedTargetMRGCcellImpulseResponseData.temporalSupportSeconds*1e3, theAchievedTargetMRGCcellImpulseResponseData.amplitude, 'r-');
    set(ax, 'XLim', [0 300]);
    legend(ax, [p1 p2], {'target' 'achieved'});


end


function visualizeInnerRetinaFilterDerivationResults(theAnalyzedTTFsFullFileName)

    load(theAnalyzedTTFsFullFileName, ...
        'theMRGCMosaic', 'theTargetRGCindex', 'targetCellImpulseResponseSource',...
        'stimParams', 'TTFparamsStruct', ...
        'innerRetinaFilterDerivationParams', ...
        'innerRetinaFilterDataStruct');

    
    if (strcmp(stimParams.stimulusShape, 'annulus'))
        theColor = [50 110 180]/255;
    else
        theColor = [240 50 80]/255;
    end
    theColor = theColor/max(theColor);
    theAlpha = 1.0;

%
% THE PHOTOCURRENTS-BASED TTF
%

    % Plot the amplitude spectrum of the photocurrents-based synthetic mRGC TTF
    hFig = figure(1001); clf;
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard tall figure', ...
        'darkScheme', true);
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
           ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, innerRetinaFilterDataStruct.photocurrentBasedTTF, 'o', ...
           false, false, theColor, theAlpha, ...
          '');
    yTicks = get(ax, 'YTick');
    yTickLabels = get(ax, 'YTickLabel');
    PublicationReadyPlotLib.applyFormat(ax,ff);
    set(ax, 'YTick', yTicks);
    set(ax, 'YTickLabel', yTickLabels);

    % Export
    visualizationPDFfileName = sprintf('mRGC_%d_photocurrentsBasedTTFamplitudeSpectrum_%s', theTargetRGCindex, stimParams.stimulusShape);
    exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
    pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

    % Generate the path if we need to
    RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);

    thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');


    % Plot the phase spectrum of the photocurrents-based synthetic mRGC TTF
    hFig = figure(1002); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
           ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, innerRetinaFilterDataStruct.photocurrentBasedTTF, 'o', ...
           false, false, theColor, theAlpha, ...
          '');

    PublicationReadyPlotLib.applyFormat(ax,ff);

    % Export
    visualizationPDFfileName = sprintf('mRGC_%d_photocurrentsBasedTTFphaseSpectrum_%s', theTargetRGCindex, stimParams.stimulusShape);
    exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
    pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

    % Generate the path if we need to
    RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);

    thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');



%
% THE TARGET TTF
%

    % Plot the amplitude spectrum of the target macaque mRGC TTF
    hFig = figure(1003); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
           ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, innerRetinaFilterDataStruct.targetMacaqueTTF, 'o', ...
           false, false, theColor, theAlpha, ...
          '');

    PublicationReadyPlotLib.applyFormat(ax,ff);

    % Export
    visualizationPDFfileName = sprintf('targetMacaque_%s_TTFamplitudeSpectrum_%s', targetCellImpulseResponseSource, stimParams.stimulusShape);
    exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
    pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

    % Generate the path if we need to
    RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);

    thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');


    % Plot the phase spectrum of the target macaque mRGC TTF
    hFig = figure(1004); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
           ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, innerRetinaFilterDataStruct.targetMacaqueTTF, 'o', ...
           false, false, theColor, theAlpha, ...
          '');

    PublicationReadyPlotLib.applyFormat(ax,ff);

    % Export
    visualizationPDFfileName = sprintf('targetMacaque_%s_TTFphaseSpectrum_%s', targetCellImpulseResponseSource, stimParams.stimulusShape);
    exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
    pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

    % Generate the path if we need to
    RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);

    thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');




%
% THE DERIVED INNER RETINA TTF
%

    % Plot the amplitude spectrum of the target macaque mRGC TTF
    hFig = figure(1003); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
           ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, innerRetinaFilterDataStruct.derivedInnerRetinaTTF, 'o', ...
           false, false, theColor, theAlpha, ...
          '');

    PublicationReadyPlotLib.applyFormat(ax,ff);

    % Export
    visualizationPDFfileName = sprintf('mRGC_%d_InnerRetina_TTFamplitudeSpectrum_%s_forTargetMacaque_%s_via_%s', theTargetRGCindex,  stimParams.stimulusShape, targetCellImpulseResponseSource, innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod);
    exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
    pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

    % Generate the path if we need to
    RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);

    thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');


    % Plot the phase spectrum of the target macaque mRGC TTF
    hFig = figure(1004); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
           ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, innerRetinaFilterDataStruct.derivedInnerRetinaTTF, 'o', ...
           false, false, theColor, theAlpha, ...
          '');

    PublicationReadyPlotLib.applyFormat(ax,ff);

    % Export
    visualizationPDFfileName = sprintf('mRGC_%d_InnerRetina_TTFphaseSpectrum_%s_forTargetMacaque_%s_via_%s', theTargetRGCindex,  stimParams.stimulusShape, targetCellImpulseResponseSource, innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod);
    
    exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
    pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

    % Generate the path if we need to
    RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);

    thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');



   
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



    % Photocurrents and derived inner retina data

    includePhotocurrentsBasedData = true;
    includeDerivedInnerRetinaData = true;
    includeTargetData = false;
    includeAchievedData = false;
    includeOriginalMacaqueMeasurementsData = false;



    % Plot the impulse response functions
    % Plot the TTF
    hFig = figure(60); clf;
    
    % The TTFs
    plotComboData(hFig, ...
        innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.photocurrentBasedTTF/max(abs(innerRetinaFilterDataStruct.photocurrentBasedTTF)), ...
        innerRetinaFilterDataStruct.targetMacaqueTTF / max(abs(innerRetinaFilterDataStruct.targetMacaqueTTF)), ...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF/max(abs(innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF)), ...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTF/max(abs(innerRetinaFilterDataStruct.derivedInnerRetinaTTF)), ...
        targetCellImpulseResponseSource, ...
        [0 1.01], [0.3 200], [0.3 1 3 10 30 100], 'log', 'frequency (Hz)', theColor, ...
        'TTFs-derived', ...
        'withMacaqueData', theMacaqueTTFData, ...
        'includePhotocurrentsBasedData', includePhotocurrentsBasedData, ...
        'includeTargetData', includeTargetData, ...
        'includeAchievedData', includeAchievedData, ...
        'includeDerivedInnerRetinaData', includeDerivedInnerRetinaData, ...
        'includeOriginalMacaqueMeasurementsData', includeOriginalMacaqueMeasurementsData);



    % The impulse responses
    hFig = figure(61); clf;
    plotComboData(hFig, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.temporalSupportSeconds*1e3, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude/max(abs(thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude)), ...
        theTargetMacaqueIR.amplitude/max(abs(theTargetMacaqueIR.amplitude)), ...
        theAchievedTargetMRGCcellImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
        theDerivedInnerRetinaImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
        targetCellImpulseResponseSource, ...
        1.01*[-1 1], [0 200], 0:25:250, 'linear', 'time (msec)', theColor, ...
        'IRFs-derived', ...
        'withMacaqueData', theMacaqueImpulseResponseData, ...
        'includePhotocurrentsBasedData', includePhotocurrentsBasedData, ...
        'includeTargetData', includeTargetData, ...
        'includeAchievedData', includeAchievedData, ...
        'includeDerivedInnerRetinaData', includeDerivedInnerRetinaData, ...
        'includeOriginalMacaqueMeasurementsData', includeOriginalMacaqueMeasurementsData);



    % target and achieved data
   includePhotocurrentsBasedData = false;
   includeDerivedInnerRetinaData = false;
   includeTargetData = true;
   includeAchievedData = true;
   includeOriginalMacaqueMeasurementsData = true;


    % Plot the impulse response functions
    % Plot the TTF
    hFig = figure(70); clf;
    
    % The TTFs
    plotComboData(hFig, ...
        innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.photocurrentBasedTTF/max(abs(innerRetinaFilterDataStruct.photocurrentBasedTTF)), ...
        innerRetinaFilterDataStruct.targetMacaqueTTF / max(abs(innerRetinaFilterDataStruct.targetMacaqueTTF)), ...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF/max(abs(innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF)), ...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTF/max(abs(innerRetinaFilterDataStruct.derivedInnerRetinaTTF)), ...
        targetCellImpulseResponseSource, ...
        [0 1.01], [0.3 200], [0.3 1 3 10 30 100], 'log', 'frequency (Hz)', theColor, ...
        'TTFs-achieved', ...
        'withMacaqueData', theMacaqueTTFData, ...
        'includePhotocurrentsBasedData', includePhotocurrentsBasedData, ...
        'includeTargetData', includeTargetData, ...
        'includeAchievedData', includeAchievedData, ...
        'includeDerivedInnerRetinaData', includeDerivedInnerRetinaData, ...
        'includeOriginalMacaqueMeasurementsData', includeOriginalMacaqueMeasurementsData);



    % The impulse responses
    hFig = figure(71); clf;
    plotComboData(hFig, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.temporalSupportSeconds*1e3, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude/max(abs(thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude)), ...
        theTargetMacaqueIR.amplitude/max(abs(theTargetMacaqueIR.amplitude)), ...
        theAchievedTargetMRGCcellImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
        theDerivedInnerRetinaImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
        targetCellImpulseResponseSource, ...
        1.01*[-1 1], [0 200], 0:25:250, 'linear', 'time (msec)', theColor, ...
        'IRFs-achieved', ...
        'withMacaqueData', theMacaqueImpulseResponseData, ...
        'includePhotocurrentsBasedData', includePhotocurrentsBasedData, ...
        'includeTargetData', includeTargetData, ...
        'includeAchievedData', includeAchievedData, ...
        'includeDerivedInnerRetinaData', includeDerivedInnerRetinaData, ...
        'includeOriginalMacaqueMeasurementsData', includeOriginalMacaqueMeasurementsData);


end


function plotComboData(hFig, ...
    theAxisData, ...
    thePhotocurrentBasedMRGCcellData, ...
    theTargetCascadedFilterData, ...
    theAchievedCascadedFilterData, ...
    theInnerRetinaData, ...
    targetCellImpulseResponseSource, ...
    yAxisLims, xAxisLims, xAxisTicks, xAxisScale, xAxisLabel, ...
    theDerivedInnerRetinaColor, pdfFileNamePostFix, varargin)

    p = inputParser;
    p.addParameter('withMacaqueData', [], @(x)(isempty(x) || (isstruct(x))));
    p.addParameter('includePhotocurrentsBasedData', false, @islogical);
    p.addParameter('includeTargetData', false, @islogical);
    p.addParameter('includeAchievedData', false, @islogical);
    p.addParameter('includeDerivedInnerRetinaData', false, @islogical);
    p.addParameter('includeOriginalMacaqueMeasurementsData', false, @islogical);

    p.parse(varargin{:});
    theMacaqueData = p.Results.withMacaqueData;
    includePhotocurrentsBasedData = p.Results.includePhotocurrentsBasedData;
    includeTargetData = p.Results.includeTargetData;
    includeAchievedData = p.Results.includeAchievedData;
    includeDerivedInnerRetinaData = p.Results.includeDerivedInnerRetinaData;
    includeOriginalMacaqueMeasurementsData = p.Results.includeOriginalMacaqueMeasurementsData;

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard tall figure', ...
        'darkScheme', true);
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    allPlotHandles = [];
    allLegends = {};

    hold(ax, 'on');

    thePhotocurrentColor = [225 120 45]/225;

    maxPhotocurrentsBasedAmplitude = [];
    if (includePhotocurrentsBasedData)
        % The photocurrents-based Data
        if (contains(pdfFileNamePostFix, 'TTFs'))
            theAmplitude = abs(thePhotocurrentBasedMRGCcellData);
        else
            theAmplitude = thePhotocurrentBasedMRGCcellData;
        end
        maxPhotocurrentsBasedAmplitude = max(abs(theAmplitude));
    
        p1 = scatter(ax, theAxisData, theAmplitude/maxPhotocurrentsBasedAmplitude, 12*12, ...
                'LineWidth', 1.5, ...
                'MarkerFaceColor', thePhotocurrentColor, ...
                'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', [1 1 1], 'MarkerEdgeAlpha', 0.5 );
        allPlotHandles(numel(allPlotHandles)+1) = p1;
        allLegends{numel(allLegends)+1} = 'photocurrents-derived';
    end

    

    if (includeAchievedData)
        % The achieved cascaded filter TTF
        if (contains(pdfFileNamePostFix, 'TTFs'))
            theAmplitude = abs(theAchievedCascadedFilterData);
        else
            theAmplitude = theAchievedCascadedFilterData;
        end
        maxAmplitude = max(abs(theAmplitude));

        if (~isempty(maxPhotocurrentsBasedAmplitude))
            theAmplitude = theAmplitude/maxPhotocurrentsBasedAmplitude;
        else
            theAmplitude = theAmplitude/maxAmplitude;
        end


        p3 = scatter(ax, theAxisData, theAmplitude, 12*12, ...
            'LineWidth', 1.5, ...
            'MarkerFaceColor', theDerivedInnerRetinaColor, ...
            'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', [1 1 1], 'MarkerEdgeAlpha', 0.5 );
        
    
        allPlotHandles(numel(allPlotHandles)+1) = p3;
        allLegends{numel(allLegends)+1} = 'photocurrent X inner retina filter cascade';
    end

    
    maxTargetAmplitude = [];
    if (includeTargetData)
        if (~isempty(theTargetCascadedFilterData))
            % The target cascaded filter TTF
            if (contains(pdfFileNamePostFix, 'TTFs'))
                theAmplitude = abs(theTargetCascadedFilterData);
            else
                theAmplitude = theTargetCascadedFilterData;
            end
            maxTargetAmplitude = max(abs(theAmplitude));
        
            plot(ax, theAxisData, theAmplitude/maxTargetAmplitude, 'k-', ...
                     'LineWidth', 4, ...
                     'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
            p2 = plot(ax,theAxisData, theAmplitude/maxTargetAmplitude, 'w-', ...
                     'LineWidth', 2, ...
                     'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
            allPlotHandles(numel(allPlotHandles)+1) = p2;
            allLegends{numel(allLegends)+1} = sprintf('%s (model)', strrep(targetCellImpulseResponseSource, 'Benardete&Kaplan 1997', 'B&K''97'));
        end
    end


    if (includeDerivedInnerRetinaData)
        % The derived inner retina filter TTF
        if (contains(pdfFileNamePostFix, 'TTFs'))
            theAmplitude = abs(theInnerRetinaData);
        else
            theAmplitude = theInnerRetinaData;
        end
        maxAmplitude = max(abs(theAmplitude));
        theAmplitude = theAmplitude/maxAmplitude;

        if (~contains(pdfFileNamePostFix, 'TTFs'))
            plot(ax,theAxisData, theAmplitude, '-', ...
                      'Color', theDerivedInnerRetinaColor, 'LineWidth', 3);
            plot(ax,theAxisData, theAmplitude, '-', ...
                      'Color', theDerivedInnerRetinaColor, 'LineWidth', 1.5);
        end
        p4 = scatter(ax, theAxisData, theAmplitude, 12*12, ...
            'LineWidth', 1.5, ...
            'MarkerFaceColor', theDerivedInnerRetinaColor, ...
            'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', [1 1 1], 'MarkerEdgeAlpha', 0.5 );
        
    
        allPlotHandles(numel(allPlotHandles)+1) = p4;
        allLegends{numel(allLegends)+1} = 'derived inner retina filter';
    end
    

    if (includeOriginalMacaqueMeasurementsData)
        if (~isempty(theMacaqueData))
            if (contains(pdfFileNamePostFix, 'TTFs'))
                theAxisData = [];
                theAmplitude = [];
            else
                theAxisData = theMacaqueData.temporalSupportSeconds * 1e3;
                theAmplitude = theMacaqueData.amplitude;
            end
            if (isempty(maxTargetAmplitude))
                maxTargetAmplitude = max(abs(theAmplitude));
            end
    
            p5 = plot(ax,theAxisData, theAmplitude/maxTargetAmplitude, 'wv', ...
                      'LineWidth', 2, ...
                      'MarkerSize', 12, 'MarkerFaceColor', [0 0 0]);
        
            allPlotHandles(numel(allPlotHandles)+1) = p5;
            allLegends{numel(allLegends)+1} = sprintf('%s (raw data)', strrep(targetCellImpulseResponseSource, 'Benardete&Kaplan 1997', 'B&K''97'));
        end
    end

    

    if (contains(pdfFileNamePostFix, 'IRFs'))
        if (contains(targetCellImpulseResponseSource, 'ON'))
            legend(ax, allPlotHandles, allLegends, 'Location', 'SouthEast', 'NumColumns', 1);
        else
            legend(ax, allPlotHandles, allLegends, 'Location', 'NorthEast', 'NumColumns', 1);
        end
    else
        if (contains(pdfFileNamePostFix, 'achieved'))
            legend(ax, allPlotHandles, allLegends, 'Location', 'SouthWest', 'NumColumns', 1);
        else
            legend(ax, allPlotHandles, allLegends, 'Location', 'West', 'NumColumns', 1);
        end
    end

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

