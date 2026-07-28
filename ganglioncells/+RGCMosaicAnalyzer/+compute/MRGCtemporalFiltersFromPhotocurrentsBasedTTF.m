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
    visualizeSinusoidalFits, ...
    onlyVisualizePreviouslySynthesizedFilters, varargin)

    p = inputParser;
    p.addParameter('importCenterTTFtoBackOut', false, @islogical);
    p.addParameter('importSurroundTTFtoBackOut', false, @islogical);
    p.addParameter('visualizeCenterSurroundContributions', true, @islogical);

    % Execute the parser
    p.parse(varargin{:});
    importCenterTTFtoBackOut = p.Results.importCenterTTFtoBackOut;
    importSurroundTTFtoBackOut = p.Results.importSurroundTTFtoBackOut;
    visualizeCenterSurroundContributions = p.Results.visualizeCenterSurroundContributions;

     % Derive theInnerRetinaTTF based on the theTargetCascadedFilterTTF and thePhotocurrentBasedMRGCcellTTF
    switch (innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod)
        case 'direct division of TTFs'

            theAnalyzedTTFsFullFileName = strrep(theAnalyzedTTFsFullFileName, '???', 'direct');

            if ((importSurroundTTFtoBackOut)&&(strcmp(stimulusShape, 'spot')))
                 theAnnulusAnalyzedTTFsFullFileName = strrep(strrep(theAnalyzedTTFsFullFileName, 'spot', 'annulus'), '_cnt_', '_srnd_');
          
                 % Read the theAnalyzedTTFsFullFileName for the annulus
                 d = load(theAnnulusAnalyzedTTFsFullFileName, 'innerRetinaFilterDataStruct');
                 assert(isfield(d.innerRetinaFilterDataStruct, 'surroundTTFtoBackOutInCenterTTFcomputation'), 'The ''annulus'' TTF must be analyzed before the ''spot'' TTF.');
                 importedSurroundTTFtoBackOut = d.innerRetinaFilterDataStruct.surroundTTFtoBackOutInCenterTTFcomputation;
            else
                 importedSurroundTTFtoBackOut = [];
            end


            if ((importCenterTTFtoBackOut)&&(strcmp(stimulusShape, 'annulus')))
                 theCenterAnalyzedTTFsFullFileName = strrep(strrep(theAnalyzedTTFsFullFileName, 'annulus', 'spot'), '_srnd_', '_cnt_');
          
                 % Read the theAnalyzedTTFsFullFileName for the annulus
                 d = load(theCenterAnalyzedTTFsFullFileName, 'innerRetinaFilterDataStruct');
                 assert(isfield(d.innerRetinaFilterDataStruct, 'centerTTFtoBackOutInCenterTTFcomputation'), 'The ''spot'' TTF must be analyzed before the ''annulus'' TTF.');
                 importedCenterTTFtoBackOut = d.innerRetinaFilterDataStruct.centerTTFtoBackOutInSurroundTTFcomputation;
            else
                importedCenterTTFtoBackOut = [];
            end



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
        if (strcmp(stimulusShape, 'annulus'))
            theColor = [50 110 180]/255;
            thePhotocurrentColor = [239 142 72]/255;
        else
            theColor = [240 50 80]/255;
            thePhotocurrentColor = [239 142 72]/255;
        end

        visualizeInnerRetinaFilterDerivationResults(theAnalyzedTTFsFullFileName, ...
            visualizeCenterSurroundContributions, thePhotocurrentColor, theColor);
        return;
    end

    % Load the measured TTF responses
    load(theMRGCMosaicTTFResponsesFullFileName, ...
        'theMRGCMosaic', 'stimParams', 'TTFparamsStruct', ...
        'computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex', ...
        'nonLinearitiesList', ...
        'computeMRGCMosaicResponsesAlsoWithCenterMechanismDeactivated', ...
        'computeMRGCMosaicResponsesAlsoWithSurroundMechanismDeactivated', ...
        'theMRGCMosaicTTFresponsesAllConditions', ...
        'theMRGCMosaicTemporalSupportSecondsAllConditions', ...
        'theMRGCMosaicTTFresponsesAllConditionsCenterDeactivated', ...
        'theMRGCMosaicTTFresponsesAllConditionsSurroundDeactivated');

    assert(strcmp(stimParams.stimulusShape, stimulusShape), ...
        'Stimulus shapes do not agree in ''stimParams'' and passed stimulusShape param.');

    % Assert that the sources are configured correctly
    assert(strcmp(stimParams.stimulusShape, TTFparamsStruct.stimulusShape), ...
        'Stimulus shapes do not agree in ''stimParams'' and ''TTFparamsStruct''.');

    assert( ...
        (contains(targetCellImpulseResponseSource, 'center') && (strcmp(stimParams.stimulusShape, 'spot'))) || ...
        (contains(targetCellImpulseResponseSource, 'surround') && (strcmp(stimParams.stimulusShape, 'annulus'))), ...
         'Stimulus shape (''%s'') does not agree with the source of the targetCellImpulse response (''%s'').', ...
         stimParams.stimulusShape, targetCellImpulseResponseSource)

    % Assert that the specified RGCindex matches computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex
    assert(theTargetRGCindex == computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex, ...
        sprintf('Different targetRGCindex specified (%d)than what the one for which the TTFresponses were computed for (%d).\nThis will result in a TTF computed a stimulus not exactly centered on RGC %d\n', ...
        theTargetRGCindex, computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex, theTargetRGCindex));



    
    theSinudoidalFitColor = [1 0.5 0];

    % Compute thePhotocurrentBasedMRGCcellTTF (composite center+surround)
    [thePhotocurrentBasedMRGCcellTTFnonNormalized, temporalFrequencySupportHz] = ...
        RGCMosaicConstructor.temporalFilterEngine.complexValuedTTFfromResponsesToSinusoidalModulations(...
            TTFparamsStruct.tfSupport, ...
            theMRGCMosaicTemporalSupportSecondsAllConditions, ...
            squeeze(theMRGCMosaicTTFresponsesAllConditions(:,:,theTargetRGCindex)), ...
            allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries, ...
            visualizeSinusoidalFits, ...
            stimParams.stimulusShape, ...
            theSinudoidalFitColor);
    
    if (computeMRGCMosaicResponsesAlsoWithCenterMechanismDeactivated)
        % Compute thePhotocurrentBasedMRGCcellTTF (center deactivated)
        thePhotocurrentBasedMRGCcellCenterDeactivatedTTFnonNormalized = ...
        RGCMosaicConstructor.temporalFilterEngine.complexValuedTTFfromResponsesToSinusoidalModulations(...
            TTFparamsStruct.tfSupport, ...
            theMRGCMosaicTemporalSupportSecondsAllConditions, ...
            squeeze(theMRGCMosaicTTFresponsesAllConditionsCenterDeactivated(:,:,theTargetRGCindex)), ...
            allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries, ...
            visualizeSinusoidalFits, ...
            stimParams.stimulusShape, ...
            theSinudoidalFitColor);
    else
        thePhotocurrentBasedMRGCcellCenterDeactivatedTTFnonNormalized = [];
    end

    if (computeMRGCMosaicResponsesAlsoWithSurroundMechanismDeactivated)
        % Compute thePhotocurrentBasedMRGCcellTTF (surround deactivated)
        thePhotocurrentBasedMRGCcellSurroundDeactivatedTTFnonNormalized = ...
        RGCMosaicConstructor.temporalFilterEngine.complexValuedTTFfromResponsesToSinusoidalModulations(...
            TTFparamsStruct.tfSupport, ...
            theMRGCMosaicTemporalSupportSecondsAllConditions, ...
            squeeze(theMRGCMosaicTTFresponsesAllConditionsSurroundDeactivated(:,:,theTargetRGCindex)), ...
            allowNonZeroBaselineInSineWaveFitsToResponseTimeSeries, ...
            visualizeSinusoidalFits, ...
            stimParams.stimulusShape, ...
            theSinudoidalFitColor);
    else
        thePhotocurrentBasedMRGCcellSurroundDeactivatedTTFnonNormalized = [];
    end



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


    % Normalize the TTFs to unit magnitude
    thePhotocurrentBasedMRGCcellTTF = thePhotocurrentBasedMRGCcellTTFnonNormalized / max(abs(thePhotocurrentBasedMRGCcellTTFnonNormalized(:)));
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




    % Derive theInnerRetinaTTF based on the theTargetCascadedFilterTTF and
    % the various photocurrents-based mRGC TTFs 

    exportedSurroundTTFtoBackOut = [];
    switch (innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod)
        case 'direct division of TTFs'

            % First thePhotocurrentBasedMRGCcellTTF (composite center+surround)
            theInnerRetinaTTF = directlyDeconvolve(theTargetCascadedFilterTTF, ...
                thePhotocurrentBasedMRGCcellTTF, []);
           
           
            if (computeMRGCMosaicResponsesAlsoWithCenterMechanismDeactivated)
                 % Only the surround mechanism is active
                
                theInnerRetinaTTFcenterDeactivated = directlyDeconvolve(theTargetCascadedFilterTTF, ...
                    thePhotocurrentBasedMRGCcellCenterDeactivatedTTFnonNormalized, []);

                if (strcmp(stimParams.stimulusShape, 'annulus'))
                    % See Equation 5 of "VSS2026: Deriving the inner retinal filter TTFs for the center
                    %                             and the surround mechanism"
                    % THIS IS WRONG. NEEDS UPDATING. Needs to be a
                    % partially deconvolved version of the macaque TTF 
                    % BUT BASED ON CONTRIBUTIONS, BETTER NOT TO DO IT.
                    exportedSurroundTTFtoBackOut =  thePhotocurrentBasedMRGCcellSurroundDeactivatedTTFnonNormalized .* theInnerRetinaTTFcenterDeactivated;
                end

            else
                theInnerRetinaTTFcenterDeactivated = [];
            end


            if (computeMRGCMosaicResponsesAlsoWithSurroundMechanismDeactivated)
                % Only center mechanism is active

                switch stimParams.stimulusShape
                    case 'annulus'
                        theInnerRetinaTTFsurroundDeactivated = directlyDeconvolve(theTargetCascadedFilterTTF, ...
                            thePhotocurrentBasedMRGCcellSurroundDeactivatedTTFnonNormalized, []);
                        
                    case 'spot'
                         if (isempty(importedSurroundTTFtoBackOut))
                             theInnerRetinaTTFsurroundDeactivated = directlyDeconvolve(theTargetCascadedFilterTTF, ...
                                thePhotocurrentBasedMRGCcellSurroundDeactivatedTTFnonNormalized, []);
                         else
                             % This should give us the center mechanism TTF only
                             theInnerRetinaTTFsurroundDeactivated = directlyDeconvolve(theTargetCascadedFilterTTF, ...
                                thePhotocurrentBasedMRGCcellSurroundDeactivatedTTFnonNormalized, importedSurroundTTFtoBackOut);
                         end

                    otherwise
                        error('stimulusShape (''%s'') is incorect. Must be either ''annulus'' or ''spot''.', stimParams.stimulusShape)
                end % switch

            else
                theInnerRetinaTTFsurroundDeactivated = [];
            end

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

            
    % Compute the finite time inner retina filter
    [theInnerRetinaTTFwithFiniteTimeIR, theInnerRetinaFiniteTimeTemporalImpulseResponseData] = ...
        RGCMosaicConstructor.temporalFilterEngine.finiteTimeIRbasedTTF(...
            temporalFrequencySupportHz, ...
            theInnerRetinaTTF, ...
            innerRetinaFilterDerivationParams.finiteTimeIRdurationSeconds, ...
            innerRetinaFilterDerivationParams.leftFiniteDurationWindowDurationSeconds, ...
            innerRetinaFilterDerivationParams.rightFiniteDurationWindowDurationSeconds);


    fprintf('Saving derived inner retina filter TTF to %s\n', theAnalyzedTTFsFullFileName);

    innerRetinaFilterDataStruct = struct(...
        'temporalFilterSynthesisMethod', innerRetinaFilterDerivationParams.temporalFilterSynthesisMethod, ...
        'targetCellImpulseResponseSource', targetCellImpulseResponseSource, ...
        'temporalFrequencySupportHz', temporalFrequencySupportHz, ...
        'targetMacaqueTTF', theTargetCascadedFilterTTF, ...
        'achievedTargetMRGCcellTTF', thePhotocurrentBasedMRGCcellTTF .* theInnerRetinaTTF, ...
        'achievedTargetMRGCcellTTFwithFiniteTimeIR', thePhotocurrentBasedMRGCcellTTF .* theInnerRetinaTTFwithFiniteTimeIR, ...
        'photocurrentBasedTTF', thePhotocurrentBasedMRGCcellTTF, ...
        'photocurrentBasedTTFnonNormalized', thePhotocurrentBasedMRGCcellTTFnonNormalized, ...
        'photocurrentBasedSurroundDeactivatedTTFnonNormalized', thePhotocurrentBasedMRGCcellSurroundDeactivatedTTFnonNormalized, ...
        'photocurrentBasedCenterDeactivatedTTFnonNormalized', thePhotocurrentBasedMRGCcellCenterDeactivatedTTFnonNormalized, ...
        'derivedInnerRetinaFiniteTimeTemporalImpulseResponseData', theInnerRetinaFiniteTimeTemporalImpulseResponseData, ...
        'derivedInnerRetinaTTF', theInnerRetinaTTF, ...
        'derivedInnerRetinaTTFfiniteTimeIR', theInnerRetinaTTFwithFiniteTimeIR, ...
        'derivedInnerRetinaTTFcenterDeactivated', theInnerRetinaTTFcenterDeactivated, ...
        'derivedInnerRetinaTTFsurroundDeactivated', theInnerRetinaTTFsurroundDeactivated, ...
        'surroundTTFtoBackOutInCenterTTFcomputation', exportedSurroundTTFtoBackOut, ...
        'fittedModelParams', modelParams, ...
        'frequencyWeights', frequencyWeights ...
        );


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


function theInnerRetinaTTF = directlyDeconvolve(theMacaqueTTF, thePhotocurrentBaseSyntheticCellTTF, photocurrentBasedTTFtoBeBackedOut)

    if (isempty(photocurrentBasedTTFtoBeBackedOut))
        idx = find(abs(thePhotocurrentBaseSyntheticCellTTF)>10*eps);
        theInnerRetinaTTF = theMacaqueTTF*0;
        theInnerRetinaTTF(idx) = theMacaqueTTF(idx)./thePhotocurrentBaseSyntheticCellTTF(idx);
    else
        % See Equation 5 of "VSS2026: Deriving the inner retinal filter TTFs for the center
        %                             and the surround mechanism"
        idx = find(abs(thePhotocurrentBaseSyntheticCellTTF)>10*eps);
        theInnerRetinaTTF = theMacaqueTTF*0;
        theInnerRetinaTTF(idx) = (theMacaqueTTF(idx) - photocurrentBasedTTFtoBeBackedOut(idx))./thePhotocurrentBaseSyntheticCellTTF(idx);
    end


end


function visualizeInnerRetinaFilterDerivationResults(theAnalyzedTTFsFullFileName, ...
    visualizeCenterSurroundContributions, thePhotocurrentColor, theColor)

    targetString = '_MRGCMosaic';
    idx = strfind(theAnalyzedTTFsFullFileName, targetString);
    theInnerRetinaFilterFiniteImpulseResponseFileName = sprintf('%s.mat',theAnalyzedTTFsFullFileName(1:idx-1));

    load(theAnalyzedTTFsFullFileName, ...
        'theMRGCMosaic', 'theTargetRGCindex', 'targetCellImpulseResponseSource',...
        'stimParams', 'TTFparamsStruct', ...
        'innerRetinaFilterDerivationParams', ...
        'innerRetinaFilterDataStruct');


    if (visualizeCenterSurroundContributions)
        hFig = figure(888);
        set(hFig, 'Name', stimParams.stimulusShape)
        subplot(1,2,1);
        p1 = plot(innerRetinaFilterDataStruct.temporalFrequencySupportHz, abs(innerRetinaFilterDataStruct.photocurrentBasedTTFnonNormalized), 'k-', 'LineWidth', 3); hold on
        p2 = plot(innerRetinaFilterDataStruct.temporalFrequencySupportHz, abs(innerRetinaFilterDataStruct.photocurrentBasedSurroundDeactivatedTTFnonNormalized), 'r-', 'LineWidth', 2);
        p3 = plot(innerRetinaFilterDataStruct.temporalFrequencySupportHz, abs(innerRetinaFilterDataStruct.photocurrentBasedCenterDeactivatedTTFnonNormalized), 'b-', 'LineWidth', 1);
        legend([p1 p2 p3], {'center+surround', 'only center', 'only surround'});
        title(sprintf('photocurrent responses to %s', stimParams.stimulusShape));
    
        subplot(1,2,2);
        plot(innerRetinaFilterDataStruct.temporalFrequencySupportHz, abs(innerRetinaFilterDataStruct.targetMacaqueTTF), 'k--');
        legend(targetCellImpulseResponseSource)
    end



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
           false, false, thePhotocurrentColor, theAlpha, ...
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
           false, false, thePhotocurrentColor, theAlpha, ...
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



    % Target + photocurrent magnitude magnitude (overlay in same plot)
    hFig = figure(2003); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};


    % The photocurrent phase plot
    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
           ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, innerRetinaFilterDataStruct.photocurrentBasedTTF, 'o', ...
           false, false, thePhotocurrentColor, theAlpha, ...
           '');
     % Add the magnitude spectrum of the target macaque mRGC TTF
    hold(ax, 'on')
    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
           ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, innerRetinaFilterDataStruct.targetMacaqueTTF, 'o', ...
           false, false, theColor, theAlpha, ...
          '');

    PublicationReadyPlotLib.applyFormat(ax,ff);

    % Export
    visualizationPDFfileName = sprintf('targetMacaque_%s_TTFmagnitudepectrum_%s_withPhotocurrentOverlay', targetCellImpulseResponseSource, stimParams.stimulusShape);
    exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
    pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

    % Generate the path if we need to
    RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);

    thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');




    % Target + photocurrent phase magnitude (overlay in same plot)
    hFig = figure(2004); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};


    % The photocurrent phase plot
    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
           ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, innerRetinaFilterDataStruct.photocurrentBasedTTF, 'o', ...
           false, false, thePhotocurrentColor, theAlpha, ...
           '');

    % Add the phase spectrum of the target macaque mRGC TTF
    hold(ax, 'on')
    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
           ax, innerRetinaFilterDataStruct.temporalFrequencySupportHz, innerRetinaFilterDataStruct.targetMacaqueTTF, 'o', ...
           false, false, theColor, theAlpha, ...
          '');

    PublicationReadyPlotLib.applyFormat(ax,ff);

    % Export
    visualizationPDFfileName = sprintf('targetMacaque_%s_TTFphaseSpectrum_%s_withPhotocurrentOverlay', targetCellImpulseResponseSource, stimParams.stimulusShape);
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

    % Plot the amplitude spectrum of the derived inner retina filter
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


    % Plot the phase spectrum of the derived inner retina filter
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

    % The achieved cascaded filter (with finite IR) mRGC impulse response 
    theAchievedTargetMRGCcellFiniteTimeImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTFwithFiniteTimeIR, ...
        innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        'causal', false); 


    % The photocurrents based mRGC impulse response 
    thePhotocurrentBasedMRGCcellImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        innerRetinaFilterDataStruct.photocurrentBasedTTF, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        'causal', false); 

    % The derived inner retina impulse response (infinite time)
    theDerivedInnerRetinaImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTF, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        'causal', false); 


    % The derived inner retina impulse response (finite time)
    theDerivedInnerRetinaFiniteTimeImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTFfiniteTimeIR, innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        'causal', false); 



    % Only macaque data

    includePhotocurrentsBasedData = ~true;
    includeDerivedInnerRetinaData = ~true;
    includeDerivedInnerRetinaDataFiniteTimeIR = ~true;
    includeTargetData = true;
    includeAchievedData = false;
    includeAchievedDataWithFiniteTimeIR = false;
    includeOriginalMacaqueMeasurementsData = ~true;


     % The TTFs
    hFig = figure(30); clf;
    
   
    plotComboData(hFig, ...
        innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.photocurrentBasedTTF/max(abs(innerRetinaFilterDataStruct.photocurrentBasedTTF)), ...
        innerRetinaFilterDataStruct.targetMacaqueTTF / max(abs(innerRetinaFilterDataStruct.targetMacaqueTTF)), ...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF/max(abs(innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF)), ...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTFwithFiniteTimeIR/max(innerRetinaFilterDataStruct.achievedTargetMRGCcellTTFwithFiniteTimeIR), ...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTF/max(abs(innerRetinaFilterDataStruct.derivedInnerRetinaTTF)), ...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTFfiniteTimeIR/max(abs(innerRetinaFilterDataStruct.derivedInnerRetinaTTFfiniteTimeIR)), ...
        targetCellImpulseResponseSource, ...
        [0 1.01], [0.3 200], [0.3 1 3 10 30 100], 'log', 'frequency (Hz)', theColor, ...
        'TTFs-macaque data', ...
        'withMacaqueData', theMacaqueTTFData, ...
        'includePhotocurrentsBasedData', includePhotocurrentsBasedData, ...
        'includeTargetData', includeTargetData, ...
        'includeAchievedData', includeAchievedData, ...
        'includeAchievedDataWithFiniteTimeIR', includeAchievedDataWithFiniteTimeIR, ...
        'includeDerivedInnerRetinaData', includeDerivedInnerRetinaData, ...
        'includeDerivedInnerRetinaDataFiniteTimeIR', includeDerivedInnerRetinaDataFiniteTimeIR, ...
        'includeOriginalMacaqueMeasurementsData', includeOriginalMacaqueMeasurementsData);


     includePhotocurrentsBasedData = ~true;
    includeDerivedInnerRetinaData = ~true;
    includeDerivedInnerRetinaDataFiniteTimeIR = ~true;
    includeTargetData = true;
    includeAchievedData = false;
    includeAchievedDataWithFiniteTimeIR = false;
    includeOriginalMacaqueMeasurementsData = true;


    % The impulse responses
    hFig = figure(31); clf;
    plotComboData(hFig, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.temporalSupportSeconds*1e3, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude/max(abs(thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude)), ...
        theTargetMacaqueIR.amplitude/max(abs(theTargetMacaqueIR.amplitude)), ...
        theAchievedTargetMRGCcellImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
        theAchievedTargetMRGCcellFiniteTimeImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellFiniteTimeImpulseResponseData.amplitude)), ...
        theDerivedInnerRetinaImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
        theDerivedInnerRetinaFiniteTimeImpulseResponseData.amplitude/max(abs(theDerivedInnerRetinaFiniteTimeImpulseResponseData.amplitude)), ...
        targetCellImpulseResponseSource, ...
        1.01*[-1 1], [0 200], 0:25:250, 'linear', 'time (msec)', theColor, ...
        'IRs-macaque data', ...
        'withMacaqueData', theMacaqueImpulseResponseData, ...
        'includePhotocurrentsBasedData', includePhotocurrentsBasedData, ...
        'includeTargetData', includeTargetData, ...
        'includeAchievedData', includeAchievedData, ...
        'includeAchievedDataWithFiniteTimeIR', includeAchievedDataWithFiniteTimeIR, ...
        'includeDerivedInnerRetinaData', includeDerivedInnerRetinaData, ...
        'includeDerivedInnerRetinaDataFiniteTimeIR', includeDerivedInnerRetinaDataFiniteTimeIR, ...
        'includeOriginalMacaqueMeasurementsData', includeOriginalMacaqueMeasurementsData);





    % Photocurrents and derived inner retina data

    includePhotocurrentsBasedData = true;
    includeDerivedInnerRetinaData = ~true;
    includeDerivedInnerRetinaDataFiniteTimeIR = true;
    includeTargetData = ~true;
    includeAchievedData = false;
    includeAchievedDataWithFiniteTimeIR = false;
    includeOriginalMacaqueMeasurementsData = false;



    % Plot the TTFs
    hFig = figure(60); clf;
    
    % The TTFs
    plotComboData(hFig, ...
        innerRetinaFilterDataStruct.temporalFrequencySupportHz, ...
        innerRetinaFilterDataStruct.photocurrentBasedTTF/max(abs(innerRetinaFilterDataStruct.photocurrentBasedTTF)), ...
        innerRetinaFilterDataStruct.targetMacaqueTTF / max(abs(innerRetinaFilterDataStruct.targetMacaqueTTF)), ...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF/max(abs(innerRetinaFilterDataStruct.achievedTargetMRGCcellTTF)), ...
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTFwithFiniteTimeIR/max(innerRetinaFilterDataStruct.achievedTargetMRGCcellTTFwithFiniteTimeIR), ...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTF/max(abs(innerRetinaFilterDataStruct.derivedInnerRetinaTTF)), ...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTFfiniteTimeIR/max(abs(innerRetinaFilterDataStruct.derivedInnerRetinaTTFfiniteTimeIR)), ...
        targetCellImpulseResponseSource, ...
        [0 1.01], [0.3 200], [0.3 1 3 10 30 100], 'log', 'frequency (Hz)', theColor, ...
        'TTFs-derived', ...
        'withMacaqueData', theMacaqueTTFData, ...
        'includePhotocurrentsBasedData', includePhotocurrentsBasedData, ...
        'includeTargetData', includeTargetData, ...
        'includeAchievedData', includeAchievedData, ...
        'includeAchievedDataWithFiniteTimeIR', includeAchievedDataWithFiniteTimeIR, ...
        'includeDerivedInnerRetinaData', includeDerivedInnerRetinaData, ...
        'includeDerivedInnerRetinaDataFiniteTimeIR', includeDerivedInnerRetinaDataFiniteTimeIR, ...
        'includeOriginalMacaqueMeasurementsData', includeOriginalMacaqueMeasurementsData);



    includePhotocurrentsBasedData = ~true;
    includeDerivedInnerRetinaData = ~true;
    includeDerivedInnerRetinaDataFiniteTimeIR = true;
    includeTargetData = ~true;
    includeAchievedData = false;
    includeAchievedDataWithFiniteTimeIR = false;
    includeOriginalMacaqueMeasurementsData = false;



    % The impulse responses
    hFig = figure(61); clf;
    plotComboData(hFig, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.temporalSupportSeconds*1e3, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude/max(abs(thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude)), ...
        theTargetMacaqueIR.amplitude/max(abs(theTargetMacaqueIR.amplitude)), ...
        theAchievedTargetMRGCcellImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
        theAchievedTargetMRGCcellFiniteTimeImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellFiniteTimeImpulseResponseData.amplitude)), ...
        theDerivedInnerRetinaImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
        theDerivedInnerRetinaFiniteTimeImpulseResponseData.amplitude/max(abs(theDerivedInnerRetinaFiniteTimeImpulseResponseData.amplitude)), ...
        targetCellImpulseResponseSource, ...
        1.01*[-1 1], [0 200], 0:25:250, 'linear', 'time (msec)', theColor, ...
        'IRFs-derived', ...
        'withMacaqueData', theMacaqueImpulseResponseData, ...
        'includePhotocurrentsBasedData', includePhotocurrentsBasedData, ...
        'includeTargetData', includeTargetData, ...
        'includeAchievedData', includeAchievedData, ...
        'includeAchievedDataWithFiniteTimeIR', includeAchievedDataWithFiniteTimeIR, ...
        'includeDerivedInnerRetinaData', includeDerivedInnerRetinaData, ...
        'includeDerivedInnerRetinaDataFiniteTimeIR', includeDerivedInnerRetinaDataFiniteTimeIR, ...
        'includeOriginalMacaqueMeasurementsData', includeOriginalMacaqueMeasurementsData);



   % target and achieved data
   includePhotocurrentsBasedData = false;
   includeDerivedInnerRetinaData = false;
   includeDerivedInnerRetinaDataFiniteTimeIR = false;
   includeTargetData = ~true;
   includeAchievedData = ~true;
   includeAchievedDataWithFiniteTimeIR = true;
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
        innerRetinaFilterDataStruct.achievedTargetMRGCcellTTFwithFiniteTimeIR/max(innerRetinaFilterDataStruct.achievedTargetMRGCcellTTFwithFiniteTimeIR), ...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTF/max(abs(innerRetinaFilterDataStruct.derivedInnerRetinaTTF)), ...
        innerRetinaFilterDataStruct.derivedInnerRetinaTTFfiniteTimeIR/max(abs(innerRetinaFilterDataStruct.derivedInnerRetinaTTFfiniteTimeIR)), ...
        targetCellImpulseResponseSource, ...
        [0 1.01], [0.3 200], [0.3 1 3 10 30 100], 'log', 'frequency (Hz)', theColor, ...
        'TTFs-achieved', ...
        'withMacaqueData', theMacaqueTTFData, ...
        'includePhotocurrentsBasedData', includePhotocurrentsBasedData, ...
        'includeTargetData', includeTargetData, ...
        'includeAchievedData', includeAchievedData, ...
        'includeAchievedDataWithFiniteTimeIR', includeAchievedDataWithFiniteTimeIR, ...
        'includeDerivedInnerRetinaData', includeDerivedInnerRetinaData, ...
        'includeDerivedInnerRetinaDataFiniteTimeIR', includeDerivedInnerRetinaDataFiniteTimeIR, ...
        'includeOriginalMacaqueMeasurementsData', includeOriginalMacaqueMeasurementsData);


    % The impulse responses
    hFig = figure(71); clf;
    plotComboData(hFig, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.temporalSupportSeconds*1e3, ...
        thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude/max(abs(thePhotocurrentBasedMRGCcellImpulseResponseData.amplitude)), ...
        theTargetMacaqueIR.amplitude/max(abs(theTargetMacaqueIR.amplitude)), ...
        theAchievedTargetMRGCcellImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
        theAchievedTargetMRGCcellFiniteTimeImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellFiniteTimeImpulseResponseData.amplitude)), ...
        theDerivedInnerRetinaImpulseResponseData.amplitude/max(abs(theAchievedTargetMRGCcellImpulseResponseData.amplitude)), ...
        theDerivedInnerRetinaFiniteTimeImpulseResponseData.amplitude/max(abs(theDerivedInnerRetinaFiniteTimeImpulseResponseData.amplitude)), ...
        targetCellImpulseResponseSource, ...
        1.01*[-1 1], [0 200], 0:25:250, 'linear', 'time (msec)', theColor, ...
        'IRFs-achieved', ...
        'withMacaqueData', theMacaqueImpulseResponseData, ...
        'includePhotocurrentsBasedData', includePhotocurrentsBasedData, ...
        'includeTargetData', includeTargetData, ...
        'includeAchievedData', includeAchievedData, ...
        'includeAchievedDataWithFiniteTimeIR', includeAchievedDataWithFiniteTimeIR, ...
        'includeDerivedInnerRetinaData', includeDerivedInnerRetinaData, ...
        'includeDerivedInnerRetinaDataFiniteTimeIR', includeDerivedInnerRetinaDataFiniteTimeIR, ...
        'includeOriginalMacaqueMeasurementsData', includeOriginalMacaqueMeasurementsData);



    % Trim data
    idx = find(theDerivedInnerRetinaFiniteTimeImpulseResponseData.temporalSupportSeconds <= ...
        innerRetinaFilterDerivationParams.finiteTimeIRdurationSeconds);
    theDerivedInnerRetinaFiniteTimeImpulseResponseData.temporalSupportSeconds = ...
        theDerivedInnerRetinaFiniteTimeImpulseResponseData.temporalSupportSeconds(idx);
    theDerivedInnerRetinaFiniteTimeImpulseResponseData.amplitude = ...
        theDerivedInnerRetinaFiniteTimeImpulseResponseData.amplitude(idx);
    
    figure(4444);
    plot(theDerivedInnerRetinaFiniteTimeImpulseResponseData.temporalSupportSeconds, theDerivedInnerRetinaFiniteTimeImpulseResponseData.amplitude, 'ko-');

    save(theInnerRetinaFilterFiniteImpulseResponseFileName, ...
        'theDerivedInnerRetinaFiniteTimeImpulseResponseData' ...
        );

    fprintf('\n******** \nDerived inner retina filter for %s saved in %s\n*********\n', ...
        stimParams.stimulusShape, theInnerRetinaFilterFiniteImpulseResponseFileName);
end


function plotComboData(hFig, ...
    theAxisData, ...
    thePhotocurrentBasedMRGCcellData, ...
    theTargetCascadedFilterData, ...
    theAchievedCascadedFilterData, ...
    theAchievedCascadedFilterFiniteTimeIRdata, ...
    theInnerRetinaData, ...
    theInnerRetinaDataFiniteTimeIR, ...
    targetCellImpulseResponseSource, ...
    yAxisLims, xAxisLims, xAxisTicks, xAxisScale, xAxisLabel, ...
    theColor, pdfFileNamePostFix, varargin)

    p = inputParser;
    p.addParameter('withMacaqueData', [], @(x)(isempty(x) || (isstruct(x))));
    p.addParameter('includePhotocurrentsBasedData', false, @islogical);
    p.addParameter('includeTargetData', false, @islogical);
    p.addParameter('includeAchievedData', false, @islogical);
    p.addParameter('includeAchievedDataWithFiniteTimeIR', false, @islogical);
    p.addParameter('includeDerivedInnerRetinaData', false, @islogical);
    p.addParameter('includeDerivedInnerRetinaDataFiniteTimeIR', false, @islogical);
    p.addParameter('includeOriginalMacaqueMeasurementsData', false, @islogical);

    p.parse(varargin{:});
    theMacaqueData = p.Results.withMacaqueData;
    includePhotocurrentsBasedData = p.Results.includePhotocurrentsBasedData;
    includeTargetData = p.Results.includeTargetData;
    includeAchievedData = p.Results.includeAchievedData;
    includeAchievedDataWithFiniteTimeIR = p.Results.includeAchievedDataWithFiniteTimeIR;
    includeDerivedInnerRetinaData = p.Results.includeDerivedInnerRetinaData;
    includeDerivedInnerRetinaDataFiniteTimeIR = p.Results.includeDerivedInnerRetinaDataFiniteTimeIR;
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
       
         if (contains(pdfFileNamePostFix, 'TTFs'))
                 yAxisLims = [0.001 1.0];
                 set(ax, 'XScale', 'log', 'XLim', [0.3 200], 'XTick', [0.1 0.3 1 3 10 30 100]);
                 set(ax, 'YScale', 'log', 'YLim', yAxisLims , 'YTick', [0.001 0.01 0.1 1], 'YTickLabel', {'.001' '.01', '.1', '1'})
         end

        allPlotHandles(numel(allPlotHandles)+1) = p1;
        allLegends{numel(allLegends)+1} = 'photocurrents-derived';
    end

   
    

    if (includeAchievedDataWithFiniteTimeIR)
        % The achieved cascaded filter TTF
        if (contains(pdfFileNamePostFix, 'TTFs'))
            theAmplitude = abs(theAchievedCascadedFilterFiniteTimeIRdata);
        else
            theAmplitude = theAchievedCascadedFilterFiniteTimeIRdata;
        end
        maxAmplitude = max(abs(theAmplitude));

        if (~isempty(maxPhotocurrentsBasedAmplitude))
            theAmplitude = theAmplitude/maxPhotocurrentsBasedAmplitude;
        else
            theAmplitude = theAmplitude/maxAmplitude;
        end


        p3a = scatter(ax, theAxisData, theAmplitude, 12*12, ...
            'LineWidth', 1.5, ...
            'MarkerFaceColor', theColor, ...
            'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', [1 1 1], 'MarkerEdgeAlpha', 0.5);
        
         if (contains(pdfFileNamePostFix, 'TTFs'))
                 yAxisLims = [0.001 1.0];
                 set(ax, 'XScale', 'log', 'XLim', [0.3 200], 'XTick', [0.1 0.3 1 3 10 30 100]);
                 set(ax, 'YScale', 'log', 'YLim', yAxisLims , 'YTick', [0.001 0.01 0.1 1], 'YTickLabel', {'.001' '.01', '.1', '1'})
          end

        allPlotHandles(numel(allPlotHandles)+1) = p3a;
        allLegends{numel(allLegends)+1} = 'photocurrent X inner retina filter cascade (finite time)';
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


        plot(ax, theAxisData, theAmplitude, '-', ...
            'LineWidth', 4.0, ...
            'Color', [1 1 1]);
        
        p3b = plot(ax, theAxisData, theAmplitude, '-', ...
            'LineWidth', 2, ...
            'Color', theColor);

         if (contains(pdfFileNamePostFix, 'TTFs'))
                 yAxisLims = [0.001 1.0];
                 set(ax, 'XScale', 'log', 'XLim', [0.3 200], 'XTick', [0.1 0.3 1 3 10 30 100]);
                 set(ax, 'YScale', 'log', 'YLim', yAxisLims , 'YTick', [0.001 0.01 0.1 1], 'YTickLabel', {'.001' '.01', '.1', '1'})
         end

        allPlotHandles(numel(allPlotHandles)+1) = p3b;
        allLegends{numel(allLegends)+1} = 'photocurrent X inner retina filter cascade (infinite time)';
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
        
            plot(ax, theAxisData, theAmplitude/maxTargetAmplitude, 'w-', ...
                     'LineWidth', 5, ...
                     'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
            p2 = plot(ax,theAxisData, theAmplitude/maxTargetAmplitude, '-', 'Color', theColor, ...
                     'LineWidth', 3, ...
                     'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);

             if (contains(pdfFileNamePostFix, 'TTFs'))
                 yAxisLims = [0.001 1.0];
                 set(ax, 'XScale', 'log', 'XLim', [0.3 200], 'XTick', [0.1 0.3 1 3 10 30 100]);
                 set(ax, 'YScale', 'log', 'YLim', yAxisLims , 'YTick', [0.001 0.01 0.1 1], 'YTickLabel', {'.001' '.01', '.1', '1'})
             end

            allPlotHandles(numel(allPlotHandles)+1) = p2;
            allLegends{numel(allLegends)+1} = sprintf('%s (model)', strrep(targetCellImpulseResponseSource, 'Benardete&Kaplan 1997', 'B&K''97'));
        end
    end


    if (includeDerivedInnerRetinaDataFiniteTimeIR)
        % The derived inner retina filter TTF
        if (contains(pdfFileNamePostFix, 'TTFs'))
            theAmplitude = abs(theInnerRetinaDataFiniteTimeIR);
        else
            theAmplitude = theInnerRetinaDataFiniteTimeIR;
        end
        maxAmplitude = max(abs(theAmplitude));
        theAmplitude = theAmplitude/maxAmplitude;

        if (~contains(pdfFileNamePostFix, 'TTFs'))
            plot(ax,theAxisData, theAmplitude, '-', ...
                      'Color', theColor, 'LineWidth', 3);
            plot(ax,theAxisData, theAmplitude, '-', ...
                      'Color', theColor, 'LineWidth', 1.5);
        end
        p4a = scatter(ax, theAxisData, theAmplitude, 12*12, ...
            'LineWidth', 1.5, ...
            'MarkerFaceColor', theColor, ...
            'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', [1 1 1], 'MarkerEdgeAlpha', 0.5 );
        
    
         if (contains(pdfFileNamePostFix, 'TTFs'))
                 yAxisLims = [0.001 1.0];
                 set(ax, 'XScale', 'log', 'XLim', [0.3 200], 'XTick', [0.1 0.3 1 3 10 30 100]);
                 set(ax, 'YScale', 'log', 'YLim', yAxisLims , 'YTick', [0.001 0.01 0.1 1], 'YTickLabel', {'.001' '.01', '.1', '1'})
         end

        allPlotHandles(numel(allPlotHandles)+1) = p4a;
        allLegends{numel(allLegends)+1} = 'derived inner retina filter (finite time)';

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

        plot(ax, theAxisData, theAmplitude, '-', ...
            'LineWidth', 4.0, ...
            'Color', [1 1 1]);
        
        p4b = plot(ax, theAxisData, theAmplitude, '-', ...
            'LineWidth', 2, ...
            'Color', theColor);


         if (contains(pdfFileNamePostFix, 'TTFs'))
                 yAxisLims = [0.001 1.0];
                 set(ax, 'XScale', 'log', 'XLim', [0.3 200], 'XTick', [0.1 0.3 1 3 10 30 100]);
                 set(ax, 'YScale', 'log', 'YLim', yAxisLims , 'YTick', [0.001 0.01 0.1 1], 'YTickLabel', {'.001' '.01', '.1', '1'})
         end
    
        allPlotHandles(numel(allPlotHandles)+1) = p4b;
        allLegends{numel(allLegends)+1} = 'derived inner retina filter (infinite time)';
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
        
             if (contains(pdfFileNamePostFix, 'TTFs'))
                 yAxisLims = [0.001 1.0];
                 set(ax, 'XScale', 'log', 'XLim', [0.3 200], 'XTick', [0.1 0.3 1 3 10 30 100]);
                 set(ax, 'YScale', 'log', 'YLim', yAxisLims , 'YTick', [0.001 0.01 0.1 1], 'YTickLabel', {'.001' '.01', '.1', '1'})
             end

            allPlotHandles(numel(allPlotHandles)+1) = p5;
            allLegends{numel(allLegends)+1} = sprintf('%s (raw data)', strrep(targetCellImpulseResponseSource, 'Benardete&Kaplan 1997', 'B&K''97'));
        end
    end

    

    if (contains(pdfFileNamePostFix, 'IRFs')) || (contains(pdfFileNamePostFix, 'IRs'))
        if (contains(pdfFileNamePostFix, 'achieved')) || (contains(pdfFileNamePostFix, 'macaque'))
            if ( ...
                    ((contains(targetCellImpulseResponseSource, 'ON')) && contains(targetCellImpulseResponseSource, 'center')) || ...
                    ((contains(targetCellImpulseResponseSource, 'OFF')) && contains(targetCellImpulseResponseSource, 'surround')) ...
                )
                legend(ax, allPlotHandles, allLegends, 'Location', 'SouthEast', 'NumColumns', 1);
            else
                legend(ax, allPlotHandles, allLegends, 'Location', 'NorthEast', 'NumColumns', 1);
            end
        elseif (contains(pdfFileNamePostFix, 'derived'))
            if (contains(targetCellImpulseResponseSource, 'ON'))
                legend(ax, allPlotHandles, allLegends, 'Location', 'SouthEast', 'NumColumns', 1);
            else
                legend(ax, allPlotHandles, allLegends, 'Location', 'NorthEast', 'NumColumns', 1);
            end
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

