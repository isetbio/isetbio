function [retinalRFparamsStruct, weightsComputeFunctionHandle, ...
          theFittedVisualRF, theFittedRetinalRFcenter, ...
          theFittedRetinalRFsurround] = fitVisualRFByAdjustingRetinalPoolingParameters(modelConstants, theVisualRF)

    
    switch (modelConstants.retinalConePoolingModel)
        case 'GaussianCenterGaussianSurroundBased'
            %                                          Kc   RcDegs   Rs/Rc integratedS/C
            retinalConePoolingParams.initialValues = [  1    0.05    5     0.5];
            retinalConePoolingParams.lowerBounds   = [ 1e-3  0.1/60  1     0.0];
            retinalConePoolingParams.upperBounds   = [ 1e5   2       100   100];

            if (modelConstants.conesNumPooledByTheRFcenter == 1)
                % Fix RcDegs
                retinalConePoolingParams.initialValues(2) = modelConstants.centerConeRcDegs;
                retinalConePoolingParams.lowerBounds(2) = modelConstants.centerConeRcDegs;
                retinalConePoolingParams.upperBounds(2) = modelConstants.centerConeRcDegs;
            end

        case 'GaussianCenterDoubleExponentSurroundBased'
            %                                          Kc   RcDegs    KsToKcPeakRatio  narrowToWideVolumeRatio  RwideDegs  RnarrowToRwideRatio
            retinalConePoolingParams.initialValues = [  1    0.05    0.3              0.3                      0.2         0.1];
            retinalConePoolingParams.lowerBounds   = [ 1e-3  0.1/60  0.00             0.01                     0.1/60      0.01];
            retinalConePoolingParams.upperBounds   = [ 1e2   2       1                1.0                      40          10.0];

            if (modelConstants.conesNumPooledByTheRFcenter == 1)
                % Fix RcDegs
                retinalConePoolingParams.initialValues(2) = modelConstants.centerConeRcDegs;
                retinalConePoolingParams.lowerBounds(2) = modelConstants.centerConeRcDegs;
                retinalConePoolingParams.upperBounds(2) = modelConstants.centerConeRcDegs;
            end

        otherwise
           error('Unknown retinalConePoolingModel: ''%s''.', retinalConePoolingModel);
    end


    % The optimization objective
    objective = @(retinalPoolingParams) sum((visualRFfromRetinalPooling(retinalPoolingParams, modelConstants) - theVisualRF).^2, 'all');

    % Ready to fit
    options = optimset(...
            'Display', 'off', ...
            'Algorithm', 'interior-point',... % 'sqp', ... % 'interior-point',...
            'GradObj', 'off', ...
            'DerivativeCheck', 'off', ...
            'MaxFunEvals', 10^5, ...
            'MaxIter', 10^3);
        
     % Multi-start
     problem = createOptimProblem('fmincon',...
          'objective', objective, ...
          'x0', retinalConePoolingParams.initialValues, ...
          'lb', retinalConePoolingParams.lowerBounds, ...
          'ub', retinalConePoolingParams.upperBounds, ...
          'options', options...
          );

     ms = MultiStart(...
          'Display', 'off', ...
          'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
          'UseParallel', true);
      
     % Run the multi-start
     multiStartsNum = 8;
     theFittedParamsVector = run(ms, problem, multiStartsNum);

     fprintf('*** Fit results****\n');
     initialParamValues = retinalConePoolingParams.initialValues 
     lowerBounds = retinalConePoolingParams.lowerBounds
     upperBounds = retinalConePoolingParams.upperBounds
     fittedParamValues = theFittedParamsVector
     
     fprintf('*****************\n');

     % Compute the fitted visual RF
     [theFittedVisualRF, weightsComputeFunctionHandle, ...
      retinalRFparamsStruct, theFittedRetinalRFcenter, ...
      theFittedRetinalRFsurround] = visualRFfromRetinalPooling(theFittedParamsVector, modelConstants);
end


function [theFittedVisualRF, weightsComputeFunctionHandle, ...
          retinalRFparamsStruct, theRetinalRFcenter, theRetinalRFsurround] = visualRFfromRetinalPooling(retinalPoolingParamsVector, modelConstants)
    
    retinalRFparamsStruct = RetinaToVisualFieldTransformer.paramsVectorToParamsStruct(...
        retinalPoolingParamsVector, modelConstants.retinalConePoolingModel);


    switch (modelConstants.retinalConePoolingModel)
        case 'GaussianCenterGaussianSurroundBased'

            % Return a handle to the function that will take the
            % retinal RF params and produce cone weights for a
            % given cone mosaic
            weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.retinalConeWeightsFromDoGmodelParameters;

            % Determine the indices and weights of cones feeding into the
            % retinal RF center and surround mechanisms. These are specific
            % to the modelConstants.theConeMosaic by applying the DoG model
            % with the above retinalRFparamsStruct params
            pooledConeIndicesAndWeightsStruct = weightsComputeFunctionHandle(...
                                        retinalRFparamsStruct, ...
                                        modelConstants.conesNumPooledByTheRFcenter, ...
                                        [], [], ...
                                        modelConstants.theConeMosaic);

        case 'GaussianCenterDoubleExponentSurroundBased'

            % Return a handle to the function that will take the
            % retinal RF params and produce cone weights for a
            % given cone mosaic
            weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.retinalConeWeightsFromDoGDEmodelParameters;

            % Determine the indices and weights of cones feeding into the
            % retinal RF center and surround mechanisms. These are specific
            % to the modelConstants.theConeMosaic by applying the DoG model
            % with the above retinalRFparamsStruct params
            pooledConeIndicesAndWeightsStruct = weightsComputeFunctionHandle(...
                                        retinalRFparamsStruct, ...
                                        modelConstants.conesNumPooledByTheRFcenter, ...
                                        [], [], ...
                                        modelConstants.theConeMosaic);
        otherwise
           error('Unknown retinalConePoolingModel: ''%s''.', retinalConePoolingModel);
    end

    % Compute the retinalRF by summing the weighted cone apertures in the
    % center and surround as specified in the computed pooledConeIndicesAndWeightsStruct
    [theRetinalRFcenter, theRetinalRFsurround] = RetinaToVisualFieldTransformer.generateRFsubregionMapsFromPooledCones(...
                       modelConstants.rfSupportX, modelConstants.rfSupportY, modelConstants.theConeMosaic, pooledConeIndicesAndWeightsStruct);
        
    % And the full cone-pooling based retinal RF
    theRetinalRF = theRetinalRFcenter - theRetinalRFsurround;

    % Convolve the cone-pooling based retinal RF to get the corresponding visual RF
    theFittedVisualRF = conv2(theRetinalRF, modelConstants.thePSF.data, 'same');
end
