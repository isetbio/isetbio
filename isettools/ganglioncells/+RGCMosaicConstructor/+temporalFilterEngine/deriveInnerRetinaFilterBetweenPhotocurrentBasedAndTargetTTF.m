%
% RGCMosaicConstructor.temporalFilterEngine.deriveInnerRetinaFilterBetweenPhotocurrentBasedAndTargetTTF
%

function [theInnerRetinaTTF, modelParams] = deriveInnerRetinaFilterBetweenPhotocurrentBasedAndTargetTTF(...
                temporalFrequencySupportHz, theTargetTTF, thePhotocurrentsBasedTTF, theIdealInnerRetinaTTF, ...
                frequencyWeights, temporalWeightingLimitsSeconds, timeDomainResidualWeighting, ...
                residualIsBasedOnTTFofCascadedPhotocurrentInnerRetinaFilter, ...
                temporalFilterSynthesisMethod, ...
                amplitudeSpectrumVsComplexSpectrumBias, ...
                solverType, multiStartsNum, useParallel)


    % Normalized target TTFs
    theTargetTTF = theTargetTTF / max(abs(theTargetTTF(:)));
    theIdealInnerRetinaTTF = theIdealInnerRetinaTTF / max(abs(theIdealInnerRetinaTTF(:)));

    theTargetImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                theTargetTTF, temporalFrequencySupportHz, ...
                'causal', false);

    theIdealInnerRetinaFilterResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                theIdealInnerRetinaTTF, temporalFrequencySupportHz, ...
                'causal', false);

    theIdealInnerRetinaFilterResponseData.amplitude = ...
        theIdealInnerRetinaFilterResponseData.amplitude/max(abs(theIdealInnerRetinaFilterResponseData.amplitude));


    mPos = max(theIdealInnerRetinaFilterResponseData.amplitude(:));
    mNeg = max(-theIdealInnerRetinaFilterResponseData.amplitude(:));
    if (mNeg > mPos)
        waveformPolarity = -1;
    else
        waveformPolarity = 1;
    end
    


    idx = find(frequencyWeights>0.0);
    minFrequencyToIncludeWithUnitWeight = temporalFrequencySupportHz(idx(1));
    maxFrequencyToIncludeWithUnitWeight = temporalFrequencySupportHz(idx(end));

    generateMagnitudePhaseSpectraPlots(temporalFrequencySupportHz, ...
        [minFrequencyToIncludeWithUnitWeight, maxFrequencyToIncludeWithUnitWeight],...
        theTargetTTF, thePhotocurrentsBasedTTF, theIdealInnerRetinaTTF);


    hFig = figure(1000); clf;
    set(hFig, 'Position', [10 10 1100 1000]);
    ax = subplot(1,1,1);

    switch (temporalFilterSynthesisMethod)
        case 'delayLeadLagFilter'
            [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = RGCMosaicConstructor.temporalFilterEngine.delayLeadLagFilter([],temporalFrequencySupportHz);

        case 'delayLeadLagFilter2'
            [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = RGCMosaicConstructor.temporalFilterEngine.delayLeadLagFilter2([],temporalFrequencySupportHz);

        case 'delayHighPassFilter'
            [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = RGCMosaicConstructor.temporalFilterEngine.delayHighPassFilter([],temporalFrequencySupportHz);

        case 'asymmetricBandPassFilter'
            [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = RGCMosaicConstructor.temporalFilterEngine.asymmetricBandPassFilter([],temporalFrequencySupportHz);

        case 'dampedOscillationFilter'
               [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter([], temporalFrequencySupportHz);

        case 'dampedOscillationLowPassCascadeFilter'
                [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = RGCMosaicConstructor.temporalFilterEngine.dampedOscillationLowPassCascadeFilter([], temporalFrequencySupportHz);


        case 'differenceOfLowPassFilters'
             [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters([], temporalFrequencySupportHz);

        case 'differenceOfLowPassFilters2'
             [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters2([], temporalFrequencySupportHz);


        case 'sumOfLowPassFilters'
             [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = RGCMosaicConstructor.temporalFilterEngine.sumOfLowPassFilters([], temporalFrequencySupportHz);


        otherwise
            error('Unknown filterTye: ''%s''.', temporalFilterSynthesisMethod);
    end

    modelParams.finalValues = modelParams.initialValues;
    modelParams.scaling = 'linear';

    saveIntermediateResults = false;
    showCurrentParamValuesPlot = true;
    showFitResult = true;

    % Figure showing the current fitted model params
    hFig = figure(2000); clf;
    set(hFig, 'Position', [1230 3750 800 950]);
    fittedModelParamsAxes = subplot(1,1,1);

    % Figure showing the current fit
    progressFigureHandle = figure(201); clf;
    set(progressFigureHandle, 'Position', [1 35 1200 1250]);


    modelParams.initialValues

    objectiveFunctionToMinimize = @(x)theObjectiveFunctionToMinimize(x, ...
        temporalFrequencySupportHz, ...
        theTargetTTF, thePhotocurrentsBasedTTF, frequencyWeights, ...
        temporalFilterSynthesisMethod, amplitudeSpectrumVsComplexSpectrumBias, ...
        fittedModelParamsAxes, progressFigureHandle, modelParams, ...
        showCurrentParamValuesPlot, showFitResult);

    

    residualsSequence = [];
    problem = createOptimProblem('fmincon',...
          'objective', objectiveFunctionToMinimize , ...
          'x0', modelParams.initialValues, ...
          'lb', modelParams.lowerBounds, ...
          'ub', modelParams.upperBounds, ...
          'options', optimoptions(...
            'fmincon',...
            'Display', 'none', ...
            'Algorithm', 'interior-point',... % 'sqp', ... % 'interior-point',...
            'GradObj', 'off', ...
            'DerivativeCheck', 'off', ...
            'MaxFunEvals', 10^5, ...
            'MaxIter', 10^4) ...
          );


     switch (solverType)
            case 'multi-start'
                % Setup the multi-start solver
                ms = MultiStart(...
                'Display', 'iter', ...
                'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
                'UseParallel', useParallel);
    
                % Run the multi-start
                modelParams.finalValues = run(ms, problem, multiStartsNum)
    
            case 'global-search'
    
                gs = GlobalSearch;
                gs = GlobalSearch(gs,'XTolerance',1e-3,'StartPointsToRun','bounds');
                modelParams.finalValues = run(gs,problem)
    
    
            case 'fmincon'
                modelParams.finalValues = fmincon(problem)
    
            otherwise
                error('Uknown solver type: ''%s''.', solverType)
                
     end % switch (solver)



    switch (temporalFilterSynthesisMethod)
        case 'delayLeadLagFilter'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.delayLeadLagFilter(...
                modelParams.finalValues, temporalFrequencySupportHz);

        case 'delayLeadLagFilter2'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.delayLeadLagFilter2(...
                modelParams.finalValues, temporalFrequencySupportHz);

        case 'delayHighPassFilter'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.delayHighPassFilter(...
                modelParams.finalValues, temporalFrequencySupportHz);

        case 'asymmetricBandPassFilter'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.asymmetricBandPassFilter(...
                modelParams.finalValues, temporalFrequencySupportHz);

        case 'dampedOscillationFilter'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter(...
                modelParams.finalValues, temporalFrequencySupportHz);
         
        case 'dampedOscillationLowPassCascadeFilter'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.dampedOscillationLowPassCascadeFilter(...
                modelParams.finalValues, temporalFrequencySupportHz);

        case 'differenceOfLowPassFilters'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters(...
                modelParams.finalValues, temporalFrequencySupportHz);

        case 'differenceOfLowPassFilters2'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters2(...
                modelParams.finalValues, temporalFrequencySupportHz);


        case 'sumOfLowPassFilters'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.sumOfLowPassFilters(...
                modelParams.finalValues, temporalFrequencySupportHz);

        otherwise
            error('Unknown filterTye: ''%s''.', temporalFilterSynthesisMethod);
    end


    theInnerRetinaTTF  = theInnerRetinaTTF  * waveformPolarity;


    % Nested function
    function theResidual = theObjectiveFunctionToMinimize(theCurrentParams, ...
        temporalFrequencySupportHz, ...
        theTargetTTF, thePhotocurrentsBasedTTF,  ...
        frequencyWeights, temporalFilterSynthesisMethod, amplitudeSpectrumVsComplexSpectrumBias, ...
        fittedModelParamsAxes, progressFigureHandle, ...
        modelParams, showCurrentParamValuesPlot, showFitResults)

        switch (temporalFilterSynthesisMethod)
            case 'delayLeadLagFilter'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    RGCMosaicConstructor.temporalFilterEngine.delayLeadLagFilter(theCurrentParams, temporalFrequencySupportHz);
    
            case 'delayLeadLagFilter2'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    RGCMosaicConstructor.temporalFilterEngine.delayLeadLagFilter2(theCurrentParams, temporalFrequencySupportHz);
    
            case 'delayHighPassFilter'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    RGCMosaicConstructor.temporalFilterEngine.delayHighPassFilter(theCurrentParams, temporalFrequencySupportHz);
    
            case 'asymmetricBandPassFilter'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    RGCMosaicConstructor.temporalFilterEngine.asymmetricBandPassFilter(theCurrentParams, temporalFrequencySupportHz);

            case 'dampedOscillationFilter'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter(theCurrentParams, temporalFrequencySupportHz);

            case 'dampedOscillationLowPassCascadeFilter'
                 [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    RGCMosaicConstructor.temporalFilterEngine.dampedOscillationLowPassCascadeFilter(theCurrentParams, temporalFrequencySupportHz);

            case 'differenceOfLowPassFilters'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters(theCurrentParams, temporalFrequencySupportHz);
       
            case 'differenceOfLowPassFilters2'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters2(theCurrentParams, temporalFrequencySupportHz);
       
            case 'sumOfLowPassFilters'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    RGCMosaicConstructor.temporalFilterEngine.sumOfLowPassFilters(theCurrentParams, temporalFrequencySupportHz);
    
        end


        theCurrentInnerRetinaFilterTTF  = theCurrentInnerRetinaFilterTTF * waveformPolarity;


        % Compute the time domain-residual
        theCurrentInnerRetinaFilterResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                theCurrentInnerRetinaFilterTTF, temporalFrequencySupportHz, ...
                'causal', false);
        theCurrentInnerRetinaFilterResponseData.amplitude = theCurrentInnerRetinaFilterResponseData.amplitude/max(abs(theCurrentInnerRetinaFilterResponseData.amplitude));

        theTimeBinsOfInterest = find((theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds>=temporalWeightingLimitsSeconds(1)) & (theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds<=temporalWeightingLimitsSeconds(2)));
        theTimeDomainResidual = norm(theIdealInnerRetinaFilterResponseData.amplitude(theTimeBinsOfInterest) - theCurrentInnerRetinaFilterResponseData.amplitude(theTimeBinsOfInterest)) / max(theIdealInnerRetinaFilterResponseData.amplitude);


        % Compute the spectral-domain redidual
        if (residualIsBasedOnTTFofCascadedPhotocurrentInnerRetinaFilter)
            desiredTTF = theTargetTTF;
            achievedTTF = theCurrentInnerRetinaFilterTTF .* thePhotocurrentsBasedTTF;

        else
            desiredTTF = theIdealInnerRetinaTTF;
            achievedTTF = theCurrentInnerRetinaFilterTTF;
        end

        
        diffAmplitude = frequencyWeights(:) .* (abs(desiredTTF(:)) - abs(achievedTTF(:)));
        diffComplex = frequencyWeights(:) .* (desiredTTF(:) - achievedTTF(:));
        theSpectralDomainResidual = amplitudeSpectrumVsComplexSpectrumBias * norm(diffAmplitude) / max(norm(desiredTTF(:))) + ...
                                   (1-amplitudeSpectrumVsComplexSpectrumBias) * norm(diffComplex) / max(norm(desiredTTF(:)));


        % Combine the two residuals
        %
        theResidual = timeDomainResidualWeighting * theTimeDomainResidual + (1-timeDomainResidualWeighting)* theSpectralDomainResidual; 
        residualsSequence(numel(residualsSequence)+1) = theResidual;


        if (theResidual == min(residualsSequence(:))) && (saveIntermediateResults)

            save(sprintf('OptimizationAtIteration_%d.mat', numel(residualsSequence)), ...
                'theCurrentParams', ...
                'temporalFrequencySupportHz', ...
                'theCurrentInnerRetinaFilterTTF', ...
                'thePhotocurrentsBasedTTF', ...
                'theTargetTTF')
        end

    
        if (showCurrentParamValuesPlot)
            % Show the current param values plot
            modelParams.finalValues = theCurrentParams;
            RGCMosaicConstructor.visualize.fittedModelParams(fittedModelParamsAxes, modelParams, temporalFilterSynthesisMethod);
        end
    
    
        if (showFitResults) && (theResidual == min(residualsSequence(:)))   


            figure(progressFigureHandle);
            clf;


            ax = subplot('Position', [0.02 0.08 0.5 0.7]);


            if (residualIsBasedOnTTFofCascadedPhotocurrentInnerRetinaFilter)
 
                p1 = plot(ax,theTargetImpulseResponseData.temporalSupportSeconds*1e3, ...
                          theTargetImpulseResponseData.amplitude, ...
                          'k-', 'LineWidth', 1.5);
                hold(ax, 'on')
                theCurrentCascadeImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                    achievedTTF, temporalFrequencySupportHz, ...
                    'causal', false);

                p2 = plot(ax,theCurrentCascadeImpulseResponseData.temporalSupportSeconds*1e3, ...
                          theCurrentCascadeImpulseResponseData.amplitude, ...
                          'r-', 'LineWidth', 1.5);

            else

                p1 = plot(ax,theIdealInnerRetinaFilterResponseData.temporalSupportSeconds*1e3, ...
                          theIdealInnerRetinaFilterResponseData.amplitude, ...
                          'k-', 'LineWidth', 1.5);
                hold(ax, 'on')
                p2 = plot(ax,theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds*1e3, ...
                          theCurrentInnerRetinaFilterResponseData.amplitude, ...
                          'r-', 'LineWidth', 1.5);
            
            end



            plot(ax,temporalWeightingLimitsSeconds(1)*1e3*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);
            plot(ax,temporalWeightingLimitsSeconds(2)*1e3*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);

            m = numel(theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds);
            set(ax, 'XLim', theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds(m)+[0 200]);
            set(ax, 'FontSize', 16);
            legend(ax, [p1 p2], {'ideal (direct deconvolution)', sprintf('fitted ''%s'' model', temporalFilterSynthesisMethod)});
            theParametersString = '';
            for iParam = 1:numel(modelParams.names)
                theParametersString = sprintf('%s%s: %g\n',theParametersString, modelParams.names{iParam}, theCurrentParams(iParam));
            end
            title(ax, theParametersString, 'FontWeight', 'normal', 'FontName', 'SourceCodePro');


            % The sequene of residuals (inset)
            ax = axes('Position', [0.25 0.13 0.25 0.3]);
            plot(ax,1:numel(residualsSequence),residualsSequence, 'b-', 'LineWidth', 1.5);
            set(ax, 'YLim', [0 residualsSequence(1)*1.2], 'Color', 'none', 'XColor', [0.4 0.4 1], 'YColor', [0.4 0.4 1]);
            box(ax, 'off');
            grid(ax, 'on')
            title(ax, sprintf('residuals: %2.2f (frequency domain, time windowed), %2.2f (time domain)', theSpectralDomainResidual, theTimeDomainResidual));

            
            ax = subplot('Position', [0.59 0.58 0.4 0.4]);
            RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
                ax, temporalFrequencySupportHz, desiredTTF, 'o', ...
                true, false, [0 0 0], 1.0, ...
                '');
            hold(ax, 'on');
            RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
                ax, temporalFrequencySupportHz, achievedTTF, '-', ...
                true, false, [1 0 0], 1.0, ...
                '');
            plot(ax, minFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);
            plot(ax, maxFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);


            ax = subplot('Position', [0.59 0.08 0.4 0.4]);
            RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
                ax, temporalFrequencySupportHz, desiredTTF, 'o', ...
                false, false, [0 0 0], 1.0, ...
                '');
            hold(ax, 'on');
            RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
                ax, temporalFrequencySupportHz, achievedTTF, '-', ...
                false, false, [1 0 0], 1.0, ...
                '');
            plot(ax, minFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);
            plot(ax, maxFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);

            drawnow;
        end
    end % Nested function

end



function generateMagnitudePhaseSpectraPlots(temporalFrequencySupportHz, ...
    frequencyWeightingLimits, ...
    theTargetTTF, thePhotocurrentsBasedTTF, theIdealInnerRetinaTTF)

    hFig = figure(999);
    set(hFig, 'Position', [10 10 1500 800]);
    ax = subplot(2,3,1);
    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
        ax, temporalFrequencySupportHz, theTargetTTF, 'o', ...
        true, false, [0 0 0], 1.0,  ...
        'target (macaque) TTF');
    % Add minimum & maximum weights line
    hold(ax, 'on');
    plot(ax, frequencyWeightingLimits(1)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);
    plot(ax, frequencyWeightingLimits(2)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);

    ax = subplot(2,3,4);
    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
        ax, temporalFrequencySupportHz, theTargetTTF, 'o', ...
        false, false, [0 0 0], 1.0, ...
        '');
    % Add minimum & maximum weights line
    hold(ax, 'on');
    plot(ax, frequencyWeightingLimits(1)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);
    plot(ax, frequencyWeightingLimits(2)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);


    ax = subplot(2,3,2);
    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
        ax, temporalFrequencySupportHz, thePhotocurrentsBasedTTF, 'o', ...
        true, true, [0 0 1], 1.0, ...
        'photocurrents-based TTF');

    % Add minimum & maximum weights line
    hold(ax, 'on');
    plot(ax, frequencyWeightingLimits(1)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);
    plot(ax, frequencyWeightingLimits(2)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);

    ax = subplot(2,3,5);
    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
        ax, temporalFrequencySupportHz, thePhotocurrentsBasedTTF, 'o', ...
        false, true, [0 0 1], 1.0, ...
        '');
    % Add minimum & maximum weights line
    hold(ax, 'on');
    plot(ax, frequencyWeightingLimits(1)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);
    plot(ax, frequencyWeightingLimits(2)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);


    ax = subplot(2,3,3);
    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
        ax, temporalFrequencySupportHz, theIdealInnerRetinaTTF, 'o', ...
        true, true, [1 0 0], 1.0, ...
        'desired inner retina TTF');

    % Add minimum & maximum weights line
    hold(ax, 'on');
    plot(ax, frequencyWeightingLimits(1)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);
    plot(ax, frequencyWeightingLimits(2)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);

    ax = subplot(2,3,6);
    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
        ax, temporalFrequencySupportHz, theIdealInnerRetinaTTF, 'o', ...
        false, true, [1 0 0], 1.0, ...
        '');

    % Add minimum & maximum weights line
    hold(ax, 'on');
    plot(ax, frequencyWeightingLimits(1)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);
    plot(ax, frequencyWeightingLimits(2)*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);
end
