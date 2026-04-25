%
% RGCMosaicConstructor.temporalFilterEngine.innerRetinaFilterBetweenPhotocurrentBasedAndTargetTTF
%

function [theInnerRetinaTTF, modelParams] = innerRetinaFilterBetweenPhotocurrentBasedAndTargetTTF(...
                temporalFrequencySupportHz, theTargetTTF, thePhotocurrentsBasedTTF, ...
                frequencyWeights, filterType, solverType, multiStartsNum, useParallel)


    theDesiredInnerRetinaTTF = theTargetTTF./thePhotocurrentsBasedTTF;
    theDesiredInnerRetinaFilterResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                theDesiredInnerRetinaTTF, temporalFrequencySupportHz);

    idx = find(frequencyWeights==1.0);
    minFrequencyToIncludeWithUnitWeight = temporalFrequencySupportHz(idx(1));
    generateMagnitudePhaseSpectraPlots(temporalFrequencySupportHz, minFrequencyToIncludeWithUnitWeight, theTargetTTF, thePhotocurrentsBasedTTF, theDesiredInnerRetinaTTF);


    hFig = figure(1000); clf;
    set(hFig, 'Position', [10 10 1100 1000]);
    ax = subplot(1,1,1);

    switch (filterType)
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

        case 'dampedOscillationFiter'
               [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter([], temporalFrequencySupportHz);

        case 'differenceOfLowPassFilters'
             [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters([], temporalFrequencySupportHz);

        otherwise
            error('Unknown filterTye: ''%s''.', filterType);
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

    objectiveFunctionToMinimize = @(x)theObjectiveFunctionToMinimize(x, temporalFrequencySupportHz, ...
        theTargetTTF, thePhotocurrentsBasedTTF, frequencyWeights, ...
        filterType, fittedModelParamsAxes, progressFigureHandle, modelParams, ...
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
            'Algorithm', 'sqp',... % 'sqp', ... % 'interior-point',...
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

  

    % For lead-lag filter order fixed to 1.0, modelParams.finalValues were found as follows:
    % 0.1306 -7.6329 1.1930 150 (MAX) 1 0.9796 20.9183 0.9796 13.9997

    % For lead-lag filter order fixed to 2.0, modelParams.finalValues were found as follows:
    % 0.0627 -0.8099 1.2076 71.8375   2 1.2076 14.4942 1.2076 19.8082

    % For lead-lag filter order fixed to 3.0, modelParams.finalValues were found as follows:
    %0.1335  -9.1854 1.6700 29.5640   3  1.6700 4.1152 1.6700 28.7999

    % For lead-lag filter order fixed to 4.0, modelParams.finalValues were found as follows:
    % 0.1853 -29.4469 1.6509 20.0019  4 1 .6509 18.0088 1.6508 31.3295


    % For lead-lag filter order fixed to 5.0, modelParams.finalValues were found as follows:
    % 0.1942 -5.9163 2.8167 15.7018 5 2.8167 8.9332 2.8166 11.6619

    % For lead-lag filter order fixed to 6.0, modelParams.finalValues were found as follows:
    % 0.2153 -13.8488 2.6219 13.5265 6 2.6218 24.4346 5 1.6532

    % For lead-lag filter order fixed to 4.0 (256 multi-starts), modelParams.finalValues were found as follows:
    % 0.1738 -9.5041 2.2091 20.0430 4 2.2091 23.6671 2.2091 3.2103

    switch (filterType)
        case 'delayLeadLagFilter'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.delayLeadLagFilter(modelParams.finalValues,temporalFrequencySupportHz);

        case 'delayLeadLagFilter2'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.delayLeadLagFilter2(modelParams.finalValues,temporalFrequencySupportHz);

        case 'delayHighPassFilter'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.delayHighPassFilter(modelParams.finalValues,temporalFrequencySupportHz);

        case 'asymmetricBandPassFilter'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.asymmetricBandPassFilter(modelParams.finalValues,temporalFrequencySupportHz);

        case 'dampedOscillationFiter'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter(modelParams.finalValues, temporalFrequencySupportHz);
         
        case 'differenceOfLowPassFilters'
            theInnerRetinaTTF = RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters(modelParams.finalValues, temporalFrequencySupportHz);

        otherwise
            error('Unknown filterTye: ''%s''.', filterType);
    end


    % Nested function
    function theResidual = theObjectiveFunctionToMinimize(theCurrentParams, ...
        temporalFrequencySupportHz, theTargetTTF, thePhotocurrentsBasedTTF,  ...
        frequencyWeights, filterType, fittedModelParamsAxes, progressFigureHandle, modelParams, showCurrentParamValuesPlot, showFitResults)

        switch (filterType)
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

            case 'dampedOscillationFiter'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    RGCMosaicConstructor.temporalFilterEngine.dampedOscillationFilter(theCurrentParams, temporalFrequencySupportHz);

            case 'differenceOfLowPassFilters'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    RGCMosaicConstructor.temporalFilterEngine.differenceOfLowPassFilters(theCurrentParams, temporalFrequencySupportHz);
        end


        theResidual = norm(frequencyWeights .* (theTargetTTF - theCurrentInnerRetinaFilterTTF.*thePhotocurrentsBasedTTF));
    
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
            RGCMosaicConstructor.visualize.fittedModelParams(fittedModelParamsAxes, modelParams, 'TTF fit');
        end
    
    
        if (showFitResults) && (theResidual == min(residualsSequence(:)))

            theCurrentInnerRetinaFilterResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                theCurrentInnerRetinaFilterTTF, temporalFrequencySupportHz);
   

            figure(progressFigureHandle);
            clf;

            ax = subplot('Position', [0.02 0.08 0.5 0.7]);
            p1 = plot(ax,theDesiredInnerRetinaFilterResponseData.temporalSupportSeconds*1e3, theDesiredInnerRetinaFilterResponseData.amplitude, ...
                'k-', 'LineWidth', 1.5);
            hold(ax, 'on')
            p2 = plot(ax,theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds*1e3, theCurrentInnerRetinaFilterResponseData.amplitude, ...
                'r-', 'LineWidth', 1.5);
            
            m = numel(theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds);
            set(ax, 'XLim', theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds(m)+[0 200]);
            set(ax, 'FontSize', 16);
            legend(ax, [p1 p2], {'direct deconvolution', sprintf('fitted ''%s'' model', filterType)});
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


            ax = subplot('Position', [0.59 0.58 0.4 0.4]);
            RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
                ax, temporalFrequencySupportHz, theDesiredInnerRetinaTTF, 'o', ...
                true, false, [0 0 0], ...
                '');
            hold(ax, 'on');
            RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
                ax, temporalFrequencySupportHz, theCurrentInnerRetinaFilterTTF, '-', ...
                true, false, [1 0 0], ...
                '');
            plot(ax, minFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);

            ax = subplot('Position', [0.59 0.08 0.4 0.4]);
            RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
                ax, temporalFrequencySupportHz, theDesiredInnerRetinaTTF, 'o', ...
                false, false, [0 0 0], ...
                '');
            hold(ax, 'on');
            RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
                ax, temporalFrequencySupportHz, theCurrentInnerRetinaFilterTTF, '-', ...
                false, false, [1 0 0], ...
                '');
            plot(ax, minFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.5);

            drawnow;
        end
    end % Nested function

end



function generateMagnitudePhaseSpectraPlots(temporalFrequencySupportHz, minFrequencyToIncludeWithUnitWeight, ...
    theTargetTTF, thePhotocurrentsBasedTTF, theDesiredInnerRetinaTTF)

    hFig = figure(999);
    set(hFig, 'Position', [10 10 1500 800]);
    ax = subplot(2,3,1);
    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
        ax, temporalFrequencySupportHz, theTargetTTF, 'o', ...
        true, false, [0 0 0], ...
        'target (macaque) TTF');

    % Add minimum weights line
    hold(ax, 'on');
    plot(ax, minFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.0);

    ax = subplot(2,3,4);
    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
        ax, temporalFrequencySupportHz, theTargetTTF, 'o', ...
        false, false, [0 0 0], ...
        '');

    % Add minimum weights line
    hold(ax, 'on');
    plot(ax, minFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'k--', 'LineWidth', 1.0);

    ax = subplot(2,3,2);
    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
        ax, temporalFrequencySupportHz, thePhotocurrentsBasedTTF, 'o', ...
        true, true, [0 0 1], ...
        'photocurrents-based TTF');

    % Add minimum weights line
    hold(ax, 'on');
    plot(ax, minFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'b--', 'LineWidth', 1.0);

    ax = subplot(2,3,5);
    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
        ax, temporalFrequencySupportHz, thePhotocurrentsBasedTTF, 'o', ...
        false, true, [0 0 1], ...
        '');

    % Add minimum weights line
    hold(ax, 'on');
    plot(ax, minFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'b--', 'LineWidth', 1.0);

    ax = subplot(2,3,3);
    RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot(...
        ax, temporalFrequencySupportHz, theDesiredInnerRetinaTTF, 'o', ...
        true, true, [1 0 0], ...
        'desired inner retina TTF');

     % Add minimum weights line
    hold(ax, 'on');
    plot(ax, minFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);

    ax = subplot(2,3,6);
    RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot(...
        ax, temporalFrequencySupportHz, theDesiredInnerRetinaTTF, 'o', ...
        false, true, [1 0 0], ...
        '');

     % Add minimum weights line
    hold(ax, 'on');
    plot(ax, minFrequencyToIncludeWithUnitWeight*[1 1], get(ax, 'YLim'), 'r--', 'LineWidth', 1.0);
end
