function theRFcomputeStruct = optimizeSurroundConePooling(obj, ...
    theRGCindex, targetVisualSTFparams,  mosaicParams, opticsParams, ...
    initialRetinalConePoolingParams, ...
    displayFittingProgress, exportedFittingProgressFolder, figNo, figTitle)
    
    % Compute optimization components
    [modelConstants, retinalConePoolingParams, visualRcDegs] = ...
        obj.computeOptimizationComponents(theRGCindex, displayFittingProgress, exportedFittingProgressFolder);

    % Override intial retinal cone pooling params
    if (~isempty(initialRetinalConePoolingParams))
        retinalConePoolingParams.initialValues = initialRetinalConePoolingParams;
    end

    rmseSequence = [];
    
    % Compute initial visual STF to make sure everything is running ok
    [theInitialRMSE, theInitialSTFdata] = ...
        theObjectiveFunction(retinalConePoolingParams.initialValues);

    % Create the optimization problem
    % Optimization optics
    options = optimoptions(...
                'fmincon',...
                'Display','iter',...
                'Algorithm','interior-point', ... %'interior-point', 'sqp'...
                'GradObj', 'off', ...
                'DerivativeCheck', 'off', ...
                'MaxFunEvals', 10^5, ...
                'MaxIter', 256, ...
                'UseParallel', true ...
    );

    if (obj.multiStartsNumRetinalPooling  > 1)
        options.Display = 'final-detailed';
    end

    problem = createOptimProblem('fmincon',...
              'objective', @theObjectiveFunction, ...
              'x0', retinalConePoolingParams.initialValues, ...
              'lb', retinalConePoolingParams.lowerBounds, ...
              'ub', retinalConePoolingParams.upperBounds, ...
              'options', options...
              );

    % Start timer 
    tic

    % Fit model
    if (obj.multiStartsNumRetinalPooling == 1)
        % Run the solver
        retinalConePoolingParams.finalValues = fmincon(problem);

        if (isempty(retinalConePoolingParams.finalValues))
            % Re-run the  solver
            fprintf(2, 'Solver failed. Trying once again\n');
            retinalConePoolingParams.finalValues = fmincon(problem);
        end

    else
        ms = MultiStart(...
              'Display', 'final', ...
              'StartPointsToRun','bounds', ...  % run only initial points that are feasible with respect to bounds
              'UseParallel', false);
    
        % Run the multi-start solver
        [retinalConePoolingParams.finalValues, ~, ~, ~, allMins] = ...
            run(ms, problem, obj.multiStartsNumRetinalPooling); 

        if (isempty(retinalConePoolingParams.finalValues))
            % Re-run the multi-start solver
            fprintf(2, 'Solver failed. Trying once again\n');
            [retinalConePoolingParams.finalValues, ~, ~, ~, allMins] = ...
                run(ms, problem, obj.multiStartsNumRetinalPooling);
        end

    end

    % Done with fitting. Report time to fit the RVFT model
    fprintf('\n===========================================\n');
    fprintf('Fitting RVFT model finished in %2.2f hours\n', toc/60/60);
    fprintf('===========================================\n');


    

    % Compute the fitted visual STF
    [theFinalRMSE, theFinalSTFdata, theFinalPooledConeIndicesAndWeights] = ...
        theObjectiveFunction(retinalConePoolingParams.finalValues);

    % Assemble theRFcomputeStruct
    theRFcomputeStruct = struct();
    theRFcomputeStruct.modelConstants = modelConstants;
    theRFcomputeStruct.retinalConePoolingParams = retinalConePoolingParams;
    theRFcomputeStruct.theTargetSTFparams = targetVisualSTFparams;
    theRFcomputeStruct.theAchievedSTFdata = theFinalSTFdata;
    theRFcomputeStruct.theFinalRMSE = theFinalRMSE;
    theRFcomputeStruct.theFinalPooledConeIndicesAndWeights = theFinalPooledConeIndicesAndWeights;
    theRFcomputeStruct.rmseSequence = rmseSequence;

    hFig = MosaicPoolingOptimizer.visualizeOptimizationProgress(figNo, figTitle, ...
                targetVisualSTFparams, opticsParams, ...
                theFinalSTFdata, retinalConePoolingParams, modelConstants.retinalConePoolingModel, ...
                theFinalPooledConeIndicesAndWeights, ...
                rmseSequence);

    [~,~,pdfsDirectory] = MosaicPoolingOptimizer.resourceFileNameAndPath('pdfsDirectory', ...
        'mosaicParams', mosaicParams);
    pdfFilename = fullfile(pdfsDirectory, sprintf('%s.pdf',figTitle));
   
    NicePlot.exportFigToPDF(pdfFilename, hFig, 300);
    close(hFig);
    
    %  ------- theObjectiveFunction --------
    function [theCurrentRMSE, theCurrentSTFdata, pooledConeIndicesAndWeights] = ...
            theObjectiveFunction(currentRetinalPoolingParamValues)

        pooledConeIndicesAndWeights = modelConstants.weightsComputeFunctionHandle(...
            modelConstants, currentRetinalPoolingParamValues);
           
        theCurrentSTFdata = obj.rgcSTFfromPooledConeMosaicSTFresponses(...
            pooledConeIndicesAndWeights, visualRcDegs);

        % Compute RMSE
        RsRcRatioResidual = theCurrentSTFdata.fittedDoGModelRsRcRatio/targetVisualSTFparams.surroundToCenterRcRatio - 1;
        SCintSensRatioResidual = theCurrentSTFdata.fittedDoGModelSCIntSensRatio/targetVisualSTFparams.surroundToCenterIntegratedSensitivityRatio - 1;
        rmseWeights = [modelConstants.rmseWeightForRsRcResidual modelConstants.rmseWeightForSCintSensResidual];    
        rmseWeights = rmseWeights / sum(rmseWeights);
        theCurrentRMSE = sqrt(rmseWeights(1) * RsRcRatioResidual^2 + rmseWeights(2)*SCintSensRatioResidual^2);
        rmseSequence(size(rmseSequence,1)+1,:) = [theCurrentRMSE RsRcRatioResidual SCintSensRatioResidual];

       % Display fitting progress
       if (displayFittingProgress)

            if (~isempty(exportedFittingProgressFolder))
                exportFigsForPaper(obj.theRGCMosaic.inputConeMosaic, ...
                    obj.theRGCMosaic.rgcRFpositionsDegs(theRGCindex,:), ...
                    theCurrentSTFdata, targetVisualSTFparams, ...
                    pooledConeIndicesAndWeights, exportedFittingProgressFolder);
                pause
            end

            retinalConePoolingParamsStruct = retinalConePoolingParams;
            retinalConePoolingParamsStruct.finalValues = currentRetinalPoolingParamValues;

            MosaicPoolingOptimizer.visualizeOptimizationProgress(figNo, figTitle, ...
                targetVisualSTFparams, opticsParams, ...
                theCurrentSTFdata, retinalConePoolingParamsStruct, modelConstants.retinalConePoolingModel, ...
                pooledConeIndicesAndWeights, ...
                rmseSequence);

       end % displayFittingProgress

    end  % ------- Nested objective function --------

end




function exportFigsForPaper(theInputConeMosaic, theCurrentRGCposition, theSTFdata, targetVisualSTFparams, pooledConeIndicesAndWeights, exportedFittingProgressFolder)
    % Generate demo figures

    % The center weights

    spatialSupportRangeArcMin = 40;
    tickSeparationArcMin = 10;

    hFig = figure(1990); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    centerLineWeightingFunctions =  mRGCMosaic.renderSubregionConePoolingPlot(theAxes{1,1}, ...
            theInputConeMosaic, ...
            theCurrentRGCposition, ...
            pooledConeIndicesAndWeights.centerConeIndices, ...
            pooledConeIndicesAndWeights.centerConeWeights, ...
            'withFigureFormat', ff, ...
            'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'plotTitle', '', ...
            'noXLabel', true, ...
            'noYLabel', true, ...
            'noYTicks', true, ...
            'xAxisTickAngleRotationDegs', 0);

    rawFiguresRoot = exportedFittingProgressFolder;
    pdfFileName = fullfile(rawFiguresRoot,'currentCenterdWeights.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

    % The current surround weights
    hFig = figure(1991); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

   
    surroundLineWeightingFunctions =  mRGCMosaic.renderSubregionConePoolingPlot(theAxes{1,1}, ...
            theInputConeMosaic, ...
            theCurrentRGCposition, ...
            pooledConeIndicesAndWeights.surroundConeIndices, ...
            pooledConeIndicesAndWeights.surroundConeWeights, ...
            'withFigureFormat', ff, ...
            'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'plotTitle', '', ...
            'noXLabel', true, ...
            'noYLabel', true, ...
            'noYTicks', true, ...
            'xAxisTickAngleRotationDegs', 0);

    rawFiguresRoot = exportedFittingProgressFolder;
    pdfFileName = fullfile(rawFiguresRoot,'currentSurroundWeights.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

    % The current line weighting functions
    hFig = figure(1992); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    normalizedPeakSurroundSensitivity = 0.4;
    sensitivityRange(2) =  max([max(centerLineWeightingFunctions.x.amplitude(:)) max(centerLineWeightingFunctions.y.amplitude(:))]);
    sensitivityRange(1) = -normalizedPeakSurroundSensitivity*sensitivityRange(2);

    centerLineWeightingFunctions.x.amplitude = centerLineWeightingFunctions.x.amplitude / max(sensitivityRange);
    centerLineWeightingFunctions.y.amplitude = centerLineWeightingFunctions.y.amplitude / max(sensitivityRange);
    surroundLineWeightingFunctions.x.amplitude = surroundLineWeightingFunctions.x.amplitude / max(sensitivityRange);
    surroundLineWeightingFunctions.y.amplitude = surroundLineWeightingFunctions.y.amplitude / max(sensitivityRange);
    sensitivityRange = sensitivityRange / max(sensitivityRange);

    mRGCMosaic.renderSubregionConePoolingLineWeightingFunctions(theAxes{1,1}, ...
            centerLineWeightingFunctions.x, surroundLineWeightingFunctions.x, ...
            sensitivityRange, 'x', ...
            'withFigureFormat', ff, ...
            'spatialSupportRangeArcMin', spatialSupportRangeArcMin, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'plotTitle', '', ...
            'noYTicks', ~true, ...
            'xAxisTickAngleRotationDegs', 0);

    rawFiguresRoot = exportedFittingProgressFolder;
    pdfFileName = fullfile(rawFiguresRoot,'currentLineWeightingFunctions.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

    % The currently achieved visualSTF
    hFig = figure(1997); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    normFactor = max(theSTFdata.fittedDoGModelToVisualSTF.centerSTFHiRes);
    
    plot(theAxes{1,1}, theSTFdata.fittedDoGModelToVisualSTF.sfHiRes, ...
                       theSTFdata.fittedDoGModelToVisualSTF.compositeSTFHiRes/normFactor, 'k-', ...
                       'LineWidth', 2*ff.lineWidth);
    hold(theAxes{1,1}, 'on');
    
    plot(theAxes{1,1}, theSTFdata.fittedDoGModelToVisualSTF.sfHiRes, ...
                       theSTFdata.fittedDoGModelToVisualSTF.centerSTFHiRes/normFactor, 'r--', ...
                       'LineWidth', ff.lineWidth);
    plot(theAxes{1,1}, theSTFdata.fittedDoGModelToVisualSTF.sfHiRes, ...
                       theSTFdata.fittedDoGModelToVisualSTF.surroundSTFHiRes/normFactor, 'b--', ...
                       'LineWidth', ff.lineWidth);

    plot(theAxes{1,1}, theSTFdata.spatialFrequencySupport, theSTFdata.visualSTF/normFactor, 'ko', ...
        'LineWidth', 1.5*ff.lineWidth, 'MarkerSize', ff.markerSize, 'MarkerFaceColor', [0.7 0.7 0.7]);        

    legend(theAxes{1,1}, {'composite (DoG fit)', 'center (DoG fit)', 'surround (DoG fit)', 'computed'}, ...
        'Location', 'SouthWest', 'FontSize', ff.legendFontSize);
    legend('boxoff')

    % Axes and limits
    set(theAxes{1,1}, 'XLim', [0.09 100], ...
        'XTick', [0.1 0.3 1 3 10 30 100], 'XTickLabel', {'.1', '.3', '1', '3', '10', '30', '100'}, 'XScale', 'log');
    set(theAxes{1,1}, 'YLim', [ff.axisOffsetFactor 1.0], 'YTick', 0:0.2:1);

    xlabel(theAxes{1,1},'spatial frequency (c/deg)', 'FontAngle', ff.axisFontAngle);
    ylabel(theAxes{1,1}, 'visual STF', 'FontAngle', ff.axisFontAngle);

    % Grid
    grid(theAxes{1,1}, 'on'); box(theAxes{1,1}, 'off');

    % Ticks
    set(theAxes{1,1}, 'TickDir', 'both');

    % Font size
    set(theAxes{1,1}, 'FontSize', ff.fontSize);

    % axis color and width
    set(theAxes{1,1}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);

    rawFiguresRoot = exportedFittingProgressFolder;
    pdfFileName = fullfile(rawFiguresRoot,'currentVisualSTF.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);


    % The currently achieved correspondence to Croner&Kaplan mean values
    hFig = figure(1997); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    achievedRsRcRatio = theSTFdata.fittedDoGModelRsRcRatio;
    achievedSCintSensRatio = theSTFdata.fittedDoGModelSCIntSensRatio;
    targetRsRcRatio = targetVisualSTFparams.surroundToCenterRcRatio;
    targetSCintSensRatio = targetVisualSTFparams.surroundToCenterIntegratedSensitivityRatio;

    bar(theAxes{1,1}, 1, 100*((achievedRsRcRatio/targetRsRcRatio)-1), 0.5, 'BaseValue', 0.0);
    hold(theAxes{1,1}, 'on');
    bar(theAxes{1,1}, 2, 100*((achievedSCintSensRatio/targetSCintSensRatio)-1), 0.5, 'BaseValue', 0.0);
    grid(theAxes{1,1}, 'on');
    box(theAxes{1,1}, 'off');

    ylabel(theAxes{1,1},'100 x (achieved/target - 1) (%)', 'FontAngle', ff.axisFontAngle);
    xlabel(theAxes{1,1}, 'STF shape metric', 'FontAngle', ff.axisFontAngle);

    % Font size
    set(theAxes{1,1}, 'FontSize', ff.fontSize);

    % axis color and width
    set(theAxes{1,1}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);

    % 0-deg roated x-ticks
    xtickangle(theAxes{1,1}, 0);

    set(theAxes{1,1}, 'XLim', [0.5+3*ff.axisOffsetFactor 2.5],  'XTick', [1 2], ...
    'XTickLabel', {sprintf('radius (srnd:cntr)'), sprintf('int. sens. (srnd:cntr)')}, ...
    'YLim', [-180+200*ff.axisOffsetFactor  180], ...
    'YTick', -150:50:150);


    rawFiguresRoot = exportedFittingProgressFolder;
    pdfFileName = fullfile(rawFiguresRoot,'currentDoGModelParamsCorrespondenceToCronnerKaplan.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
    
    % Generate paper-ready figures (scaled versions of the figures in
    % the rawFiguresRoot directory) which are stored in the PaperReady folder
    dd = strrep(rawFiguresRoot, 'Raw', 'cpdf');

    PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
    commandString = sprintf('%s -args %s/generatePLOSOnePaperReadyFigures.txt', dd, PLOSdirectory);
    system(commandString);


       

end
