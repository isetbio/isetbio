function theRFcomputeStruct = retinalConePoolingParamsForTargetVisualRF(obj, ...
                centerConeMajorityType, initialRetinalConePoolingParamsStruct, progressFigureName)

    displayFittedModel = true;
    displayFittingProgress = true;
    recordProgressVideo = false;

    indicesOfConesPooledByTheRFcenter = obj.targetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter;
    weightsOfConesPooledByTheRFcenter = obj.targetVisualRFDoGparams.weightsOfConesPooledByTheRFcenter;

    switch (centerConeMajorityType)
        case cMosaic.LCONE_ID
            theRFCenterConeMajorityPSF  = obj.spectrallyWeightedPSFData.LconeWeighted;

            spatialPositionDegs = mean(obj.coneMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);
            if (isempty(progressFigureName))
                figureName = sprintf('RF center with %d L-cone(s) at (%2.2f,%2.2f) degs', ...
                    numel(indicesOfConesPooledByTheRFcenter), spatialPositionDegs(1), spatialPositionDegs(2));
            else
                figureName = sprintf('%s (L)', progressFigureName);
            end
            figNo = 1000 + numel(weightsOfConesPooledByTheRFcenter)*10+1;
            
            
        case cMosaic.MCONE_ID
            theRFCenterConeMajorityPSF  = obj.spectrallyWeightedPSFData.MconeWeighted;

            spatialPositionDegs = mean(obj.coneMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);
            if (isempty(progressFigureName))
                figureName = sprintf('RF center with %d M-cone(s) at (%2.2f,%2.2f) degs', ...
                    numel(indicesOfConesPooledByTheRFcenter), spatialPositionDegs(1), spatialPositionDegs(2));
            else
                figureName = sprintf('%s (M)', progressFigureName);
            end

            figNo = 1000 + numel(weightsOfConesPooledByTheRFcenter)*10+2;
     

        otherwise
            error('Not L or M cone: %d', centerConeType);
    end

    % Spatial support
    spatialSupportDegs = [...
        obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs(:) ...
        obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs(:)];

    [Xdegs,Ydegs] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));
    Rdegs2 = Xdegs.^2+Ydegs.^2;

    % Retinal cone pooling model constants
    modelConstants = struct();

    % Spatial support
    modelConstants.spatialSupportDegs = spatialSupportDegs;
    modelConstants.Rdegs2 = Rdegs2;

    % The cone mosaic and the spectrally-weighted PSFs
    modelConstants.theConeMosaic = obj.coneMosaic;
    modelConstants.theRFCenterConeMajorityPSF = theRFCenterConeMajorityPSF;
    modelConstants.theSurroundLconePlusMconePSF = obj.spectrallyWeightedPSFData.LMconeWeighted;

    % The connectable cone types to the center and surroud
    modelConstants.surroundConnectableConeTypes = obj.targetVisualRFDoGparams.surroundConnectableConeTypes;
    modelConstants.centerConnectableConeTypes = obj.targetVisualRFDoGparams.centerConnectableConeTypes;

    
   if (strcmp(obj.targetVisualRFDoGparams.retinalConePoolingModel, 'arbitrary center cone weights, double exponential surround from H1 cell with index 1'))
       fprintf(2,'OLD H1 MODEL NAMING USED. Translating')
         obj.targetVisualRFDoGparams.retinalConePoolingModel = 'arbitraryCenterConeWeights_doubleExpH1cellIndex1SurroundWeights';
   end
    
    switch (obj.targetVisualRFDoGparams.retinalConePoolingModel)
        case 'arbitraryCenterConeWeights_doubleExpH1cellIndex1SurroundWeights'

            modelConstants.weightsComputeFunctionHandle = @RTVF.conePoolingCoefficientsForArbitraryCenterDoubleExpSurround;
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;

            % From the 4 cells in Figure 6 of Packer & Dacey (2002)
            RnarrowToRwideRatios  = [152/515 170/718 115/902 221/1035];
            NWvolumeRatios        = [1.0     0.8     0.3     0.2];

            % Range for RwDegs, based to characteristic radius of cones and # of cones in RF center
            RwDegsInitial    = 4.00 * 1/mean(RnarrowToRwideRatios) * obj.anatomicalRFcenterCharacteristicRadiusDegs * sqrt(2.3) * sqrt(numel(modelConstants.indicesOfCenterCones));
            RwDegsLowerBound = max([0.02 0.01*RwDegsInitial]);
            RwDegsUpperBound = min([max(spatialSupportDegs(:))*2 4*RwDegsInitial]);

            %                                        Kc      Ks/KcRatio    narrowToWideFieldVolumeRatio  RwideDegs            RnarrowToRwideRatio
            retinalConePoolingParams.names =         {'Kc',  'KsKcRatio',  'VnVwRatio',                  'RwDegs',             'RnRwRatio'};
            retinalConePoolingParams.scaling =       {'log', 'log',        'log',                        'linear',                'log'};
            retinalConePoolingParams.initialValues = [1.       0.06        mean(NWvolumeRatios)           RwDegsInitial         mean(RnarrowToRwideRatios)];
            retinalConePoolingParams.lowerBounds   = [0.5      0.005       min(NWvolumeRatios)            RwDegsLowerBound      min(RnarrowToRwideRatios)];
            retinalConePoolingParams.upperBounds   = [2        1e0         max(NWvolumeRatios)            RwDegsUpperBound      max(RnarrowToRwideRatios)];

            H1cellIndex = 1;
            parameterTolerance = 0.1;
                
            idx = find(ismember(retinalConePoolingParams.names, 'VnVwRatio'));
            measuredValue = NWvolumeRatios(H1cellIndex);
            retinalConePoolingParams.initialValues(idx) = measuredValue;
            retinalConePoolingParams.lowerBounds(idx) = measuredValue - parameterTolerance;
            retinalConePoolingParams.upperBounds(idx) = measuredValue + parameterTolerance;

            idx = find(ismember(retinalConePoolingParams.names, 'RnRwRatio'));
            measuredValue = RnarrowToRwideRatios(H1cellIndex);
            retinalConePoolingParams.initialValues(idx) = measuredValue;
            retinalConePoolingParams.lowerBounds(idx) = measuredValue - parameterTolerance;
            retinalConePoolingParams.upperBounds(idx) = measuredValue + parameterTolerance;

        otherwise
            error('Retinal cone pooling model ''%s'' is not implemented\n', obj.targetVisualRFDoGparams.retinalConePoolingModel)
    end % switch (obj.targetVisualRFDoGparams.retinalConePoolingModel)


    % Override initial params if we were passed some
    if (isempty(initialRetinalConePoolingParamsStruct))
        fprintf('\n*** RTVF will be fitted using default initial params. *** \n');
    else
        % Ensure that the imported params are for the current model
        for iParam = 1:numel(retinalConePoolingParams.names)
            assert(strcmp(retinalConePoolingParams.names{iParam}, initialRetinalConePoolingParamsStruct.names{iParam}), ...
                'model parameter name mismatch');
        end

        % All good, so replace the default initial params,
        %   retinalConePoolingParams.initialValues with
        % with the imported final values,
        %   initialRetinalConePoolingParamsStruct.finalValues

        fprintf('\n*** RTVF will be fitted with initial params set to the final values of the previously fitted RTVFobj. ***\n');
        retinalConePoolingParams.initialValues = initialRetinalConePoolingParamsStruct.finalValues;
    end

    
    videoOBJ = [];
    rmseSequence = [];

    if (displayFittingProgress)
        hFigProgress = figure(figNo); clf;
        set(hFigProgress, 'Position', [10 10 1500 950], 'Name', figureName);
        
        if (recordProgressVideo)
            videoFileName = 'ProgressVideo';
            videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
            videoOBJ.FrameRate = 10;
            videoOBJ.Quality = 100;
            videoOBJ.open();           
        end

    end


    % Compute initial visual RF map
    [theInitialRMSE, theInitiaVisualRFmap , theInitialSTFdata] = theObjectiveFunction(retinalConePoolingParams.initialValues);


    if (displayFittedModel)
        targetParams = struct(...
            'RcDegs', obj.visualRFcenterRcDegs, ...
            'RsRcRatio', obj.targetVisualRFDoGparams.surroundToCenterRcRatio, ...
            'SCintSensRatio', obj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio);

        RTVF.displayFittedModel(100, modelConstants, theInitiaVisualRFmap, theInitialSTFdata, ...
            targetParams);
    end

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

    if (obj.multiStartsNumRetinalPooling == 1)
        retinalConePoolingParams.finalValues = fmincon(problem);
    else
        ms = MultiStart(...
              'Display', 'final', ...
              'StartPointsToRun','bounds', ...  % run only initial points that are feasible with respect to bounds
              'UseParallel', false);
    
        % Run the multi-start solver
        [retinalConePoolingParams.finalValues, ~, ~, ~, allMins] = run(ms, problem, obj.multiStartsNumRetinalPooling); 
    end

    % Done with fitting. Report time to fit the RVFT model
    fprintf('\n===========================================\n');
    fprintf('Fitting RVFT model finished in %2.2f hours\n', toc/60/60);
    fprintf('===========================================\n');

     % Compute the fitted visual RF map
    [theFinalRMSE, theFinalVisualRFmap, theFinalSTFdata, theFinalPooledConeIndicesAndWeights] = ...
        theObjectiveFunction(retinalConePoolingParams.finalValues);


    if (displayFittedModel)
        RTVF.displayFittedModel(101, modelConstants, theFinalVisualRFmap, theFinalSTFdata, ...
            targetParams);
    end

    if (~isempty(videoOBJ))
        videoOBJ.close();
    end


    % Assemble theRFcomputeStruct
    theRFcomputeStruct = struct();
    theRFcomputeStruct.modelConstants = modelConstants;
    theRFcomputeStruct.retinalConePoolingParams = retinalConePoolingParams;
    theRFcomputeStruct.rmseSequence = rmseSequence;

    %  ------- Nested objective function --------
    function [theCurrentRMSE, theCurrentVisualRFmap, theCurrentSTFdata, pooledConeIndicesAndWeights] = ...
            theObjectiveFunction(currentRetinalPoolingParamValues)

        if (strcmp(obj.stfComputeMethod, RTVF.modeledSTFcomputeMethod))
            % Compute the STF data for the  current visual RF map
            % First, compute the visual RF map
            [theCurrentVisualRFmap, theCurrentRetinalRFcenterConeMap, ...
             theCurrentRetinalRFsurroundConeMap, pooledConeIndicesAndWeights] = ...
                obj.visualRFfromRetinalConePooling(...
                    modelConstants, currentRetinalPoolingParamValues);
  
            % Then compute the visual STF from the visualRF map
            theCurrentSTFdata = obj.visualSTFfromCronerKaplanAnalysisOfVisualRF(theCurrentVisualRFmap);
        else
            % Compute the indices and retinal weights of cones pooled by the center and the
            % surround mechanism
            pooledConeIndicesAndWeights = modelConstants.weightsComputeFunctionHandle(modelConstants, currentRetinalPoolingParamValues);
            theCurrentSTFdata = obj.visualSTFfromCronerKaplanAnalysisOfconeMosaicSTFresponses(pooledConeIndicesAndWeights);
            
            % No computation of visual RF maps
            theCurrentVisualRFmap = [];
            theCurrentRetinalRFcenterConeMap = [];
            theCurrentRetinalRFsurroundConeMap = [];
        end

        % Compute RMSE
        RsRcRatioResidual = theCurrentSTFdata.fittedDoGModelRsRcRatio/obj.targetVisualRFDoGparams.surroundToCenterRcRatio - 1;
        SCintSensRatioResidual = theCurrentSTFdata.fittedDoGModelSCIntSensRatio/obj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio - 1;
        rmseWeights = abs([1 1]);    
        rmseWeights = rmseWeights / sum(rmseWeights);
        theCurrentRMSE = sqrt(rmseWeights(1) * RsRcRatioResidual^2 + rmseWeights(2)*SCintSensRatioResidual^2);

        % Display fitting progress
        if (displayFittingProgress)
            rmseSequence = RTVF.displayFittingProgress(hFigProgress, videoOBJ, rmseSequence, ...
                RsRcRatioResidual, SCintSensRatioResidual, theCurrentRMSE, ...
                retinalConePoolingParams, currentRetinalPoolingParamValues, ...
                theCurrentSTFdata, ...
                obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs, ...
                obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs, ...
                theCurrentRetinalRFcenterConeMap, ...
                theCurrentRetinalRFsurroundConeMap);
        end  % displayFittingProgress

    end
end
