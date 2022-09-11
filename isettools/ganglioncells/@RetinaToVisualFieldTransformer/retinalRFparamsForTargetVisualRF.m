function retinalRFparamsForTargetVisualRF(obj, indicesOfConesPooledByTheRFcenter, ...
    weightsOfConesPooledByTheRFcenter, targetVisualRFDoGparams)
    
    % Spatial support
    spatialSupportDegs = [obj.thePSFData.supportXdegs(:) obj.thePSFData.supportYdegs(:)];
    [Xdegs,Ydegs] = meshgrid(obj.thePSFData.supportXdegs(:), obj.thePSFData.supportYdegs(:));
    Rdegs2 = Xdegs.^2+Ydegs.^2;

    % Compute the visual RF center and its characteristic radius
    [visualRFcenterConeMap, visualRFcenterCharacteristicRadiusDegs, ...
     visualRFcenterCharacteristicRadiiDegs, visualRFcenterFlatTopExponents, ...
     visualRFcenterXYpos, visualRFcenterOrientationDegs] = ...
       obj.analyzeRFcenter(indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter, spatialSupportDegs);

    
    modelConstants = struct();
    modelConstants.spatialSupportDegs = spatialSupportDegs;
    modelConstants.Rdegs2 = Rdegs2;
    
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
        case 'arbitrary center cone weights, gaussian surround weights'
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterGaussianSurround;

            %                                         Kc   Ks    RsDegs
            retinalConePoolingParams.initialValues = [1    1e-2  visualRFcenterCharacteristicRadiusDegs*6];
            retinalConePoolingParams.lowerBounds   = [1e-3 1e-4  visualRFcenterCharacteristicRadiusDegs/15];
            retinalConePoolingParams.upperBounds   = [1e3  1e-1  visualRFcenterCharacteristicRadiusDegs*15];

        case 'arbitrary center cone weights, gaussian surround weights with adjustments'
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterGaussianAdjustSurround;
            
            % Radially adjust surround cone weights at distances up to
            % 3*max(visualRFcenterCharacteristicRadiiDegs) from the RF
            % center using 10 sample points
            adjustmentSamplesNum = 10;
            adjustmentRangeDegs = 3*max(visualRFcenterCharacteristicRadiiDegs);
            modelConstants.arbitrarySurroundCorrectionRadialSupportDegs = linspace(0, adjustmentRangeDegs, adjustmentSamplesNum);

            %                                         Kc   Ks    RsDegs
            retinalConePoolingParams.initialValues = [1    1e-2  visualRFcenterCharacteristicRadiusDegs*6];
            retinalConePoolingParams.lowerBounds   = [1e-3 1e-4  visualRFcenterCharacteristicRadiusDegs/15];
            retinalConePoolingParams.upperBounds   = [1e3  1e-1  visualRFcenterCharacteristicRadiusDegs*15];

            % Surround adjustment factors
            nn = 1:numel(modelConstants.arbitrarySurroundCorrectionRadialSupportDegs);
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
    
    % Do a dry run to make sure all components are running first
    if (obj.doDryRunFirst)
        theInitialFittedVisualRFmap = RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(...\
            modelConstants, retinalConePoolingParams.initialValues);
        visualizeRF(55, spatialSupportDegs, targetVisualRFmap, 'target RF');
        visualizeRF(56, spatialSupportDegs, theInitialFittedVisualRFmap, 'target RF');
        pause();
    end


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
    retinalConePoolingParams.final = run(ms, problem, obj.multiStartsNum);

    % Compute the fitted visual RF
    [theFittedVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap] = ...
        RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(modelConstants, retinalConePoolingParams.final);

    % Form the rfComputeStruct
    obj.rfComputeStruct = struct();
    obj.rfComputeStruct.theFittedVisualRFMap = theFittedVisualRF;

    % Convolve the center retinal cone-pooling based retinal RF to get the corresponding visual center RF
    obj.rfComputeStruct.theFittedVisualRFcenterConeMap = conv2(theRetinalRFcenterConeMap, modelConstants.thePSF, 'same');
    obj.rfComputeStruct.theFittedVisualRFsurroundConeMap = conv2(theRetinalRFsurroundConeMap, modelConstants.thePSF, 'same');
    
    obj.rfComputeStruct.theRetinalRFcenterConeMap = theRetinalRFcenterConeMap;
    obj.rfComputeStruct.theRetinalRFsurroundConeMap = theRetinalRFsurroundConeMap;
    obj.rfComputeStruct.retinalConePoolingParams = retinalConePoolingParams.final;
    obj.rfComputeStruct.targetVisualRFMap = targetVisualRFmap;
    obj.rfComputeStruct.targetVisualRFcenterMap = targetVisualRFcenterMap;
    obj.rfComputeStruct.targetVisualRFsurroundMap = targetVisualRFsurroundMap;
    obj.rfComputeStruct.modelConstants = modelConstants;

    % Nested function rfObjective
    function rmse = matchedVisualRFobjective(currentRetinalPoolingParams)
        fittedVisualRF = RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(modelConstants, currentRetinalPoolingParams);
        fullRMSE = ((fittedVisualRF(:) - targetVisualRFmap(:))).^2;
        rmse =  sqrt(mean(fullRMSE,1));
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



