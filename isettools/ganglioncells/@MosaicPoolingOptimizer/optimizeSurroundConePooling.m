function theRFcomputeStruct = optimizeSurroundConePooling(obj, ...
    theRGCindex, targetVisualSTFparams,  mosaicParams, opticsParams, ...
    initialRetinalConePoolingParams, ...
    displayFittingProgress, figNo, figTitle)
    
    % Compute optimization components
    [modelConstants, retinalConePoolingParams, visualRcDegs] = ...
        obj.computeOptimizationComponents(theRGCindex);

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
