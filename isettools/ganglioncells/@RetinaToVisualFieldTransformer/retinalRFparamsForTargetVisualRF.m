function dataOut = retinalRFparamsForTargetVisualRF(obj,indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter, targetVisualRFDoGparams)
    
    % Spatial support
    spatialSupportDegs = [obj.thePSFData.supportXdegs(:) obj.thePSFData.supportYdegs(:)];
    [Xdegs,Ydegs] = meshgrid(obj.thePSFData.supportXdegs(:), obj.thePSFData.supportYdegs(:));
    Rdegs2 = Xdegs.^2+Ydegs.^2;

    % Compute the visual RF center and its characteristic radius
    [visualRFcenterConeMap, visualRFcenterCharacteristicRadiusDegs] = ...
       obj.analyzeRFcenter(indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter, spatialSupportDegs);

    % Target visual RF map model
    % Choose between {'gaussian center, gaussian surround', and 'arbitrary center, gaussian surround'}
    targetVisualMapScheme = 'gaussian center, gaussian surround';

    % Retinal cone pooling model
    retinalConePoolingModel = 'arbitrary center cone weights, gaussian surround weights';

    modelConstants = struct();
    modelConstants.spatialSupportDegs = spatialSupportDegs;
    modelConstants.Rdegs2 = Rdegs2;
    modelConstants.targetVisualMapScheme = targetVisualMapScheme;

    switch (targetVisualMapScheme)
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
            visualRFparamsVector(3) = targetVisualRFDoGparams.surroundToCenterRcRatio;
            visualRFparamsVector(4) = targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;

        otherwise
            error('Unknown targetVisualMapScheme: ''%s''.', targetVisualMapScheme);
            
    end % switch


    % Retinal cone pooling model constants
    modelConstants.theConeMosaic = obj.theConeMosaic;
    modelConstants.thePSF = obj.theCircularPSFData.data;
    
    switch (retinalConePoolingModel)
        case 'arbitrary center cone weights, gaussian surround weights'
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterGaussianSurround;

            %                                          Ks RsDegs
            retinalConePoolingParams.initialValues = [0.01 visualRFcenterCharacteristicRadiusDegs*6];
            retinalConePoolingParams.lowerBounds   = [0   visualRFcenterCharacteristicRadiusDegs/15];
            retinalConePoolingParams.upperBounds   = [1   visualRFcenterCharacteristicRadiusDegs*15];

        otherwise
            error('Unknown retinalConePoolingModel: ''%s''.', retinalConePoolingModel);
    end


    % Compute the target visual RF map
    targetVisualRFmap = targetRFfunctionHandle(modelConstants, visualRFparamsVector);
    
    
    
    debugMode = ~true;
    if (debugMode)
        % Dry run:
        theInitialFittedVisualRFmap = RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(modelConstants, retinalConePoolingParams.initialValues);

        visualizeRF(55, spatialSupportDegs, targetVisualRFmap, 'target RF');
        visualizeRF(56, spatialSupportDegs, theInitialFittedVisualRFmap, 'target RF');
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
    multiStartsNum = 10;
    retinalConePoolingParams.final = run(ms, problem, multiStartsNum);


    % Compute the fitted visual RF
    theFittedVisualRF = RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(modelConstants, retinalConePoolingParams.final);

    visualizeRF(55, spatialSupportDegs, targetVisualRFmap, 'target RF');
    visualizeRF(56, spatialSupportDegs, theFittedVisualRF, 'fitted RF');
    visualizeRF(57, spatialSupportDegs, targetVisualRFmap-theFittedVisualRF, 'target -fitted RF');

    % Form output data struct
    dataOut = struct();
    dataOut.theFittedVisualRF = theFittedVisualRF;
    dataOut.retinalConePoolingParams = retinalConePoolingParams.final;
    dataOut.targetVisualRFmap = targetVisualRFmap;
    dataOut.modelConstants = modelConstants;

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

        subplot(1,2,1)
        imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RFmap);
        set(gca, 'CLim', [-1 1])
        axis 'image'
        title(theTitle)
        colormap(gray)
    
        subplot(1,2,2)
        plot(spatialSupportDegs(:,1),  profileRF, 'r-');
        axis 'square'
       
end



