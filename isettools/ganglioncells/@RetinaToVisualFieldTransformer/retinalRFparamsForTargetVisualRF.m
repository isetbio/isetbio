function retinalRFparamsForTargetVisualRF(obj, indicesOfConesPooledByTheRFcenter, ...
    weightsOfConesPooledByTheRFcenter, targetVisualRFDoGparams)
    
    % Spatial support
    spatialSupportDegs = [obj.thePSFData.supportXdegs(:) obj.thePSFData.supportYdegs(:)];
    [Xdegs,Ydegs] = meshgrid(obj.thePSFData.supportXdegs(:), obj.thePSFData.supportYdegs(:));
    Rdegs2 = Xdegs.^2+Ydegs.^2;

    % Compute the visual RF center and its characteristic radius
    [visualRFcenterConeMap, visualRFcenterCharacteristicRadiusDegs, ...
     visualRFcenterCharacteristicRadiiDegs, visualRFcenterFlatTopExponents, ...
     visualRFcenterXYpos, visualRFcenterOrientationDegs] = obj.analyzeRFcenter(...
               indicesOfConesPooledByTheRFcenter, ...
               weightsOfConesPooledByTheRFcenter, ...
               spatialSupportDegs);

    
    modelConstants = struct();
    modelConstants.spatialSupportDegs = spatialSupportDegs;
    modelConstants.Rdegs2 = Rdegs2;

    % The weight given to regions of the targetVisualRFmap to which only
    % the surround mechanism is active
    modelConstants.surroundWeightBias = 0.0;

    switch (targetVisualRFDoGparams.visualRFmodel)
        case 'arbitrary center, gaussian surround'
            % The method handle
            targetRFfunctionHandle = @RetinaToVisualFieldTransformer.differenceOfArbitraryCenterAndGaussianSurroundRF;

            % The method free params
            visualRFparamsVector(1) = 1;
            visualRFparamsVector(2) = targetVisualRFDoGparams.surroundToCenterRcRatio;
            visualRFparamsVector(3) = targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;

            % Additional modelConstants for this model
            modelConstants.visualRFcenterCharacteristicRadiusDegs = visualRFcenterCharacteristicRadiusDegs;
            modelConstants.visualRFcenterConeMap = visualRFcenterConeMap;
            

        case 'gaussian center, gaussian surround'
            % The method handle
            targetRFfunctionHandle = @RetinaToVisualFieldTransformer.differenceOfGaussianCenterAndGaussianSurroundRF;

            % The method free params
            visualRFparamsVector(1) = 1;
            visualRFparamsVector(2) = visualRFcenterCharacteristicRadiusDegs;
            visualRFparamsVector(3) = visualRFcenterXYpos(1) + 1j*visualRFcenterXYpos(2);
            visualRFparamsVector(4) = targetVisualRFDoGparams.surroundToCenterRcRatio;
            visualRFparamsVector(5) = targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;

        case 'ellipsoidal gaussian center, gaussian surround'
            % The method handle
            targetRFfunctionHandle = @RetinaToVisualFieldTransformer.differenceOfEllipsoidalGaussianCenterAndGaussianSurroundRF;

            % The method free params
            visualRFparamsVector(1) = 1;
            visualRFparamsVector(2) = visualRFcenterCharacteristicRadiiDegs(1);
            visualRFparamsVector(3) = visualRFcenterCharacteristicRadiiDegs(2);
            visualRFparamsVector(4) = visualRFcenterXYpos(1) + 1j*visualRFcenterXYpos(2);
            visualRFparamsVector(5) = visualRFcenterFlatTopExponents(1);
            visualRFparamsVector(6) = visualRFcenterFlatTopExponents(2);
            visualRFparamsVector(7) = visualRFcenterOrientationDegs;
            visualRFparamsVector(8) = targetVisualRFDoGparams.surroundToCenterRcRatio;
            visualRFparamsVector(9) = targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;

        otherwise
            error('Unknown targetVisualMapScheme: ''%s''.', targetVisualMapScheme);
            
    end % switch


    % Retinal cone pooling model constants
    modelConstants.theConeMosaic = obj.theConeMosaic;
    modelConstants.thePSF = obj.theCircularPSFData.data;
    
    switch (targetVisualRFDoGparams.retinalConePoolingModel)
        case 'arbitrary center cone weights, double exponential surround weights'
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterDoubleExpSurround;

            %                                         Kc   Ks/KcRatio   narrowToWideFieldVolumeRatio     RwideDegs                           RnarrowToRwideRatio
            retinalConePoolingParams.names =       {'Kc', 'Ks/KcRatio', 'narrowToWideFieldVolumeRatio', 'RwideDegs',                        'RnarrowToRwideRatio'};
            retinalConePoolingParams.scaling =     {'log', 'log', 'linear', 'linear', 'linear'};
            retinalConePoolingParams.initialValues = [1      1e-1       0.6                          visualRFcenterCharacteristicRadiusDegs*6      0.1];
            retinalConePoolingParams.lowerBounds   = [1e-3   1e-3       0.2                          visualRFcenterCharacteristicRadiusDegs/5     0.01];
            retinalConePoolingParams.upperBounds   = [1e3    1e1         1.0                         visualRFcenterCharacteristicRadiusDegs*15     1.0];

        case 'arbitrary center cone weights, gaussian surround weights'
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterGaussianSurround;

            %                                         Kc   Ks    RsDegs
            retinalConePoolingParams.names         = {'Kc', 'Ks', 'RsDegs'};
            retinalConePoolingParams.scaling       = {'log', 'log', 'linear'};
            retinalConePoolingParams.initialValues = [1    1e-2  visualRFcenterCharacteristicRadiusDegs*6];
            retinalConePoolingParams.lowerBounds   = [1e-3 1e-4  visualRFcenterCharacteristicRadiusDegs/15];
            retinalConePoolingParams.upperBounds   = [1e3  1     visualRFcenterCharacteristicRadiusDegs*15];

        case 'arbitrary center cone weights, double gaussian surround weights'
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterDoubleGaussianSurround;

            %                                         Kc     Ks      RsDegs                                        Kwide      RsWideToRsRatio
            retinalConePoolingParams.names         = {'Kc',  'Ks',   'RsDegs',                                    'Kwide',    'RsWideToRsRatio' };
            retinalConePoolingParams.scaling       = {'log', 'log',  'linear',                                    'linear',   'linear'};
            retinalConePoolingParams.initialValues = [1      1e-2    visualRFcenterCharacteristicRadiusDegs*6      0.2        2.0];
            retinalConePoolingParams.lowerBounds   = [1e-3   1e-4    visualRFcenterCharacteristicRadiusDegs/15     0.0        1.0];
            retinalConePoolingParams.upperBounds   = [1e3    1       visualRFcenterCharacteristicRadiusDegs*15     1.0        10.0];

        case 'arbitrary center cone weights, gaussian surround weights with adjustments'
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterGaussianAdjustSurround;
            
            % Radially adjust surround cone weights at distances up to
            % 3*max(visualRFcenterCharacteristicRadiiDegs) from the RF
            % center using 10 sample points
            adjustmentSamplesNum = 16;
            adjustmentRangeDegs = 7*max(visualRFcenterCharacteristicRadiiDegs);
            modelConstants.arbitrarySurroundCorrectionRadialSupportDegs = linspace(0, adjustmentRangeDegs, adjustmentSamplesNum);

            %                                         Kc   Ks    RsDegs
            retinalConePoolingParams.names = {'Kc', 'Ks', 'RsDegs'};
            retinalConePoolingParams.scaling = {'log', 'log', 'linear'};
            retinalConePoolingParams.initialValues = [1    1e-2  visualRFcenterCharacteristicRadiusDegs*6];
            retinalConePoolingParams.lowerBounds   = [1e-3 1e-4  visualRFcenterCharacteristicRadiusDegs/15];
            retinalConePoolingParams.upperBounds   = [1e3  1e-1  visualRFcenterCharacteristicRadiusDegs*15];

            % Surround adjustment factors
            nn = 1:numel(modelConstants.arbitrarySurroundCorrectionRadialSupportDegs);
            for i = 1:numel(modelConstants.arbitrarySurroundCorrectionRadialSupportDegs)
                retinalConePoolingParams.names{numel(retinalConePoolingParams.names)+1} = sprintf('f-%d', i);
                retinalConePoolingParams.scaling{numel(retinalConePoolingParams.scaling)+1} = 'linear';
            end

            extraIndices = 3+nn;
            retinalConePoolingParams.initialValues(extraIndices) = 0.00*ones(1,numel(nn));
            retinalConePoolingParams.lowerBounds(extraIndices) = -0.5*ones(1,numel(nn));
            retinalConePoolingParams.upperBounds(extraIndices) =  0.5*ones(1,numel(nn));

        otherwise
            error('Unknown retinalConePoolingModel: ''%s''.', retinalConePoolingModel);
    end


    % Compute the target visual RF map
    [targetVisualRFmap, targetVisualRFcenterMap, targetVisualRFsurroundMap] = ...
        targetRFfunctionHandle(modelConstants, visualRFparamsVector);
    
    if (modelConstants.surroundWeightBias > 0)
        % Isolate target visual RF map points that are < 0 (surround only)
        decrementsIndices = find(targetVisualRFmap(:)<0);
        targetVisualRFmapDecrements = targetVisualRFmap(decrementsIndices);
        maxDecrementsRFmap = max(abs(targetVisualRFmapDecrements(:)));
        maxFullRFmap = max(abs(targetVisualRFmap(:)));
    end

    % Do a dry run to make sure all components are running first
    if (obj.doDryRunFirst)
        theInitialFittedVisualRFmap = RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(...
            modelConstants, retinalConePoolingParams.initialValues);
        visualizeRF(55, spatialSupportDegs, targetVisualRFmap, 'target RF');
        visualizeRF(56, spatialSupportDegs, theInitialFittedVisualRFmap, 'initial fitted RF');
        pause();
    end

    if (obj.multiStartsNum == 1)
        % Local search
        retinalConePoolingParams.finalValues = lsqcurvefit(RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling, ...
            retinalConePoolingParams.initialValues, ...
            modelConstants, targetVisualRFmap, ...
            pretinalConePoolingParams.lowerBounds,...
            pretinalConePoolingParams.upperBounds);
    else
        % Get ready to fit: options
        options = optimset(...
            'Display', 'off', ...
            'Algorithm', 'interior-point',... % 'sqp', ... % 'interior-point',...
            'GradObj', 'off', ...
            'DerivativeCheck', 'off', ...
            'MaxFunEvals', 10^5, ...
            'MaxIter', 10^3);

        % Multi-start
        problem = createOptimProblem('fmincon',...
          'objective', @matchedVisualRFobjective, ...
          'x0', retinalConePoolingParams.initialValues, ...
          'lb', retinalConePoolingParams.lowerBounds, ...
          'ub', retinalConePoolingParams.upperBounds, ...
          'options', options...
          );

        % Generate multi-start problem
        ms = MultiStart(...
          'Display', 'off', ...
          'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
          'UseParallel', true);

        % Run the multi-start to obtain the optimal retinalRFparamsVector
        retinalConePoolingParams.finalValues = run(ms, problem, obj.multiStartsNum);
    end


    % Visualize the fitted params
    RetinaToVisualFieldTransformer.visualizeFittedParamValues(retinalConePoolingParams);

    % Compute the fitted visual RF
    [theFittedVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap] = ...
        RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(modelConstants, retinalConePoolingParams.finalValues);

    % Form the rfComputeStruct
    obj.rfComputeStruct = struct();
    obj.rfComputeStruct.theFittedVisualRFMap = theFittedVisualRF;

    % Convolve the center retinal cone-pooling based retinal RF to get the corresponding visual center RF
    obj.rfComputeStruct.theFittedVisualRFcenterConeMap = conv2(theRetinalRFcenterConeMap, modelConstants.thePSF, 'same');
    obj.rfComputeStruct.theFittedVisualRFsurroundConeMap = conv2(theRetinalRFsurroundConeMap, modelConstants.thePSF, 'same');
    
    obj.rfComputeStruct.theRetinalRFcenterConeMap = theRetinalRFcenterConeMap;
    obj.rfComputeStruct.theRetinalRFsurroundConeMap = theRetinalRFsurroundConeMap;
    obj.rfComputeStruct.retinalConePoolingParams = retinalConePoolingParams;
    obj.rfComputeStruct.targetVisualRFMap = targetVisualRFmap;
    obj.rfComputeStruct.targetVisualRFcenterMap = targetVisualRFcenterMap;
    obj.rfComputeStruct.targetVisualRFsurroundMap = targetVisualRFsurroundMap;
    obj.rfComputeStruct.modelConstants = modelConstants;

    % Nested function rfObjective
    function rmse = matchedVisualRFobjective(currentRetinalPoolingParams)
        % Compute the fitted visual RF given the current cone pooling params
        fittedVisualRF = RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(modelConstants, currentRetinalPoolingParams);
        
        % Compute RMSE
        if (modelConstants.surroundWeightBias > 0)
            fullRMSE = (1.0/maxFullRFmap * (fittedVisualRF(:) - targetVisualRFmap(:))).^2;
            decrRMSE = (1.0/maxDecrementsRFmap * (fittedVisualRF(decrementsIndices) - targetVisualRFmapDecrements)).^2;
            rmse =  sqrt(mean(fullRMSE,1)) + modelConstants.surroundWeightBias * sqrt(mean(decrRMSE,1));
        else
            fullRMSE = ((fittedVisualRF(:) - targetVisualRFmap(:))).^2;
            rmse = sqrt(mean(fullRMSE,1));
        end
    end
end



function visualizeRF(figNo, spatialSupportDegs, RFmap, theTitle)
    figure(figNo); clf;
    profileRF = sum(RFmap,1);

    subplot(1,3,1)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RFmap);
    set(gca, 'CLim', [-1 1])
    axis 'image'
    title(theTitle)
    colormap(gray)

    subplot(1,3,2)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RFmap);
    set(gca, 'CLim', 0.01*[-1 1])
    axis 'image'

    subplot(1,3,3)
    plot(spatialSupportDegs(:,1),  profileRF, 'r-');
    axis 'square'
end



