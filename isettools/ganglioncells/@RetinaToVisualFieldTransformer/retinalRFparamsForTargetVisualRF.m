function retinalRFparamsForTargetVisualRF(obj, indicesOfConesPooledByTheRFcenter, ...
    weightsOfConesPooledByTheRFcenter, targetVisualRFDoGparams)
    
    % Spatial support
    spatialSupportDegs = [obj.thePSFData.supportXdegs(:) obj.thePSFData.supportYdegs(:)];
    [Xdegs,Ydegs] = meshgrid(obj.thePSFData.supportXdegs(:), obj.thePSFData.supportYdegs(:));
    Rdegs2 = Xdegs.^2+Ydegs.^2;

    % Compute the visual RF center and its characteristic radius
    [visualRFcenterCharacteristicRadiusDegs, visualRFcenterConeMap, ...
     visualRFcenterCharacteristicRadiiDegs, visualRFcenterFlatTopExponents, ...
     visualRFcenterXYpos, visualRFcenterOrientationDegs, ...
     anatomicalRFcenterCharacteristicRadiusDegs] = obj.analyzeRFcenter(...
               indicesOfConesPooledByTheRFcenter, ...
               weightsOfConesPooledByTheRFcenter, ...
               spatialSupportDegs);

    
    if (obj.simulateCronerKaplanEstimation)
        if (~(strcmp(targetVisualRFDoGparams.visualRFmodel, 'gaussian center, gaussian surround')))
            error('When simulating the Croner & Kaplan estimation, the targetVisualRFDoGparams.visualRFmodel must be set to ''gaussian center, gaussian surround''.');
        end
    end


    % Model constants
    modelConstants = struct();
    
    % Spatial support
    modelConstants.spatialSupportDegs = spatialSupportDegs;
    modelConstants.Rdegs2 = Rdegs2;

    

    switch (targetVisualRFDoGparams.visualRFmodel)

        case 'gaussian center, gaussian surround'
            % The target visual RF compute handle
            targetVisualRFfunctionHandle = @RetinaToVisualFieldTransformer.differenceOfGaussianCenterAndGaussianSurroundRF;

            % The target visual RF params
            targetVisualRFparamsVector(1) = max(visualRFcenterConeMap(:));
            targetVisualRFparamsVector(2) = visualRFcenterCharacteristicRadiusDegs;
            targetVisualRFparamsVector(3) = visualRFcenterXYpos(1) + 1j*visualRFcenterXYpos(2);
            targetVisualRFparamsVector(4) = targetVisualRFDoGparams.surroundToCenterRcRatio;
            targetVisualRFparamsVector(5) = targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;
        
        case 'ellipsoidal gaussian center, gaussian surround'
            % The target visual RF compute handle
            targetVisualRFfunctionHandle = @RetinaToVisualFieldTransformer.differenceOfEllipsoidalGaussianCenterAndGaussianSurroundRF;

            % The target visual RF params
            targetVisualRFparamsVector(1) = max(visualRFcenterConeMap(:));
            targetVisualRFparamsVector(2) = visualRFcenterCharacteristicRadiiDegs(1);
            targetVisualRFparamsVector(3) = visualRFcenterCharacteristicRadiiDegs(2);
            targetVisualRFparamsVector(4) = visualRFcenterXYpos(1) + 1j*visualRFcenterXYpos(2);
            targetVisualRFparamsVector(5) = visualRFcenterFlatTopExponents(1);
            targetVisualRFparamsVector(6) = visualRFcenterFlatTopExponents(2);
            targetVisualRFparamsVector(7) = visualRFcenterOrientationDegs;
            targetVisualRFparamsVector(8) = targetVisualRFDoGparams.surroundToCenterRcRatio;
            targetVisualRFparamsVector(9) = targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;

        case 'arbitrary center, gaussian surround'
            % The target visual RF compute handle
            targetVisualRFfunctionHandle = @RetinaToVisualFieldTransformer.differenceOfArbitraryCenterAndGaussianSurroundRF;

            % The target visual RF params
            targetVisualRFparamsVector(1) = max(visualRFcenterConeMap(:));
            targetVisualRFparamsVector(2) = targetVisualRFDoGparams.surroundToCenterRcRatio;
            targetVisualRFparamsVector(3) = targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;

            % Additional modelConstants for this model
            modelConstants.visualRFcenterCharacteristicRadiusDegs = visualRFcenterCharacteristicRadiusDegs;
            modelConstants.visualRFcenterConeMap = visualRFcenterConeMap;

        otherwise
            error('Unknown targetVisualMapScheme: ''%s''.', targetVisualMapScheme);
            
    end % switch


    % Retinal cone pooling model constants
    modelConstants.theConeMosaic = obj.theConeMosaic;
    modelConstants.thePSF = obj.thePSFData.data;
    

    switch (targetVisualRFDoGparams.retinalConePoolingModel)
        case {'arbitrary center cone weights, double exponential surround weights-free', ...
              'arbitrary center cone weights, double exponential surround weights-meanVnVwRatio', ...
              'arbitrary center cone weights, double exponential surround weights-meanRnRwRatio'}
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterDoubleExpSurround;

            % From the 4 cells in Figure 6 of Packer & Dacey (2002)
            RnarrowToRwideRatios  = [152/515 170/718 115/902 221/1035];
            NWvolumeRatios        = [1.0     0.8     0.3     0.2];

            %                                        Kc   Ks/KcRatio   narrowToWideFieldVolumeRatio  RwideDegs                         RnarrowToRwideRatio
            retinalConePoolingParams.names =         {'Kc', 'KsKcRatio', 'VnVwRatio',               'RwDegs',                        'RnRwRatio'};
            retinalConePoolingParams.scaling =       {'log', 'log',      'log',                        'linear',                           'log'};
            retinalConePoolingParams.initialValues = [1      0.05       mean(NWvolumeRatios)         anatomicalRFcenterCharacteristicRadiusDegs*40       mean(RnarrowToRwideRatios)];
            retinalConePoolingParams.lowerBounds   = [1e-1   1e-6       min(NWvolumeRatios)          anatomicalRFcenterCharacteristicRadiusDegs          min(RnarrowToRwideRatios)];
            retinalConePoolingParams.upperBounds   = [10     1e-1         max(NWvolumeRatios)          anatomicalRFcenterCharacteristicRadiusDegs*100      max(RnarrowToRwideRatios)];

            if (~isempty(strfind(targetVisualRFDoGparams.retinalConePoolingModel,'meanVnVwRatio')))
                idx = find(ismember(retinalConePoolingParams.names, 'VnVwRatio'));
                retinalConePoolingParams.initialValues(idx) = mean(NWvolumeRatios);
                retinalConePoolingParams.lowerBounds(idx) = mean(NWvolumeRatios);
                retinalConePoolingParams.upperBounds(idx) = mean(NWvolumeRatios);
            end
            if (~isempty(strfind(targetVisualRFDoGparams.retinalConePoolingModel,'meanRnRwRatio')))
                idx = find(ismember(retinalConePoolingParams.names, 'RnarrowToRwideRatio'));
                retinalConePoolingParams.initialValues(idx) = mean(RnarrowToRwideRatios);
                retinalConePoolingParams.lowerBounds(idx) = mean(RnarrowToRwideRatios);
                retinalConePoolingParams.upperBounds(idx) = mean(RnarrowToRwideRatios);
            end

        case 'arbitrary center cone weights, gaussian surround weights'
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterGaussianSurround;

            %                                         Kc   Ks    RsDegs
            retinalConePoolingParams.names         = {'Kc', 'Ks', 'RsDegs'};
            retinalConePoolingParams.scaling       = {'log', 'log', 'linear'};
            retinalConePoolingParams.initialValues = [1    1e-2  anatomicalRFcenterCharacteristicRadiusDegs*6];
            retinalConePoolingParams.lowerBounds   = [1e-1 1e-6  anatomicalRFcenterCharacteristicRadiusDegs];
            retinalConePoolingParams.upperBounds   = [10    10     anatomicalRFcenterCharacteristicRadiusDegs*100];

        case 'arbitrary center cone weights, variable exponential surround weights'
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterGaussianSurround;

            %                                         Kc   Ks    RsDegs  exponent
            retinalConePoolingParams.names         = {'Kc', 'Ks', 'RsDegs', 'exponent'};
            retinalConePoolingParams.scaling       = {'log', 'log', 'linear', 'linear'};
            retinalConePoolingParams.initialValues = [1    1e-2  anatomicalRFcenterCharacteristicRadiusDegs*6    1.5];
            retinalConePoolingParams.lowerBounds   = [1e-1 1e-4  anatomicalRFcenterCharacteristicRadiusDegs*3    1];
            retinalConePoolingParams.upperBounds   = [10    10     anatomicalRFcenterCharacteristicRadiusDegs*100  2];


        case 'arbitrary center cone weights, double gaussian surround weights'
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterDoubleGaussianSurround;

            %                                         Kc     Ks      Kwide      RwideDegs                                  RnarrowToRwideRatio
            retinalConePoolingParams.names         = {'Kc',  'Ks',   'Kw',      'RwDegs',                                 'RnRwRatio' };
            retinalConePoolingParams.scaling       = {'log', 'log',  'linear',  'linear',                                      'linear'};
            retinalConePoolingParams.initialValues = [1      1e-2    0.2         anatomicalRFcenterCharacteristicRadiusDegs*6   0.3];
            retinalConePoolingParams.lowerBounds   = [1e-1   1e-4    0.0         anatomicalRFcenterCharacteristicRadiusDegs*2   0.1];
            retinalConePoolingParams.upperBounds   = [10     10       1.0         anatomicalRFcenterCharacteristicRadiusDegs*100  0.5];

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
            retinalConePoolingParams.names         = {'Kc', 'Ks', 'RsDegs'};
            retinalConePoolingParams.scaling       = {'log', 'log', 'linear'};
            retinalConePoolingParams.initialValues = [1    1e-2  anatomicalRFcenterCharacteristicRadiusDegs*6];
            retinalConePoolingParams.lowerBounds   = [1e-1 1e-4  anatomicalRFcenterCharacteristicRadiusDegs*2];
            retinalConePoolingParams.upperBounds   = [10    10  anatomicalRFcenterCharacteristicRadiusDegs*100];

            % Surround adjustment factors
            nn = 1:numel(modelConstants.arbitrarySurroundCorrectionRadialSupportDegs);
            for i = 1:numel(modelConstants.arbitrarySurroundCorrectionRadialSupportDegs)
                retinalConePoolingParams.names{numel(retinalConePoolingParams.names)+1} = sprintf('F-%d', i);
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
        targetVisualRFfunctionHandle(modelConstants, targetVisualRFparamsVector);
   

    % RF assessment type
    modelConstants.simulateCronerKaplanEstimation = obj.simulateCronerKaplanEstimation;
    if (modelConstants.simulateCronerKaplanEstimation)
        theTargetRFmeasurement = sum(targetVisualRFmap,1);
        [stfSupport, theTargetSTF] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            modelConstants.spatialSupportDegs(:,1),theTargetRFmeasurement);
    else
        theTargetRFmeasurement = targetVisualRFmap(:);
    end

    

    % Do a dry run to make sure all components are running first
    if (obj.doDryRunFirst)
        theInitialFittedVisualRFmap = RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(...
            modelConstants, retinalConePoolingParams.initialValues);
        visualizeRFs(55, spatialSupportDegs, targetVisualRFmap, theInitialFittedVisualRFmap, {'target RF',  'initial fitted RF'});
        
        if (modelConstants.simulateCronerKaplanEstimation)

            theInitialFittedRFmeasurement = sum(theInitialFittedVisualRFmap,1);
            [~, theInitialFittedSTF] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
                    modelConstants.spatialSupportDegs(:,1),theInitialFittedRFmeasurement);

            hFig = figure(57); clf;
            set(hFig, 'Name', 'DryRun');
            subplot(1,2,1);
            plot(spatialSupportDegs(:,1), theTargetRFmeasurement, 'k-');
            hold on;
             plot(spatialSupportDegs(:,1), theInitialFittedRFmeasurement, 'r-');
            xlabel('space (degs)');
            
            xlabel('space (degs)');
            legend({'target RF', 'initial RF'})

            subplot(1,2,2);
            plot(stfSupport, theTargetSTF, 'k-');
            hold on;
            plot(stfSupport, theInitialFittedSTF, 'r-');

            xlabel('spatial frequency (c/deg)');
            set(gca, 'XScale', 'log');
             legend({'target RF', 'initial RF'})
        end
        pause(1)

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
        [retinalConePoolingParams.finalValues, ~, ~, ~, allMins] = run(ms, problem, obj.multiStartsNum);
        allMins
        allMins(1).X
        allMins(1).X0
    end

    % Gain of retinal cone pooling params must be 1
    
    
    % Visualize the fitted params
    RetinaToVisualFieldTransformer.visualizeFittedParamValues(retinalConePoolingParams);

    % Compute the fitted visual RF
    [theFittedVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap] = ...
        RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(modelConstants, retinalConePoolingParams.finalValues);


    if (modelConstants.simulateCronerKaplanEstimation)
        theFittedRFmeasurement = sum(theFittedVisualRF,1);
        [stfSupport, theFittedSTF] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
                 modelConstants.spatialSupportDegs(:,1),theFittedRFmeasurement);

        theTargetRFmeasurement = sum(targetVisualRFmap,1);
        [~, theTargetSTF] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
            modelConstants.spatialSupportDegs(:,1),theTargetRFmeasurement);
    end

    % Form the rfComputeStruct
    obj.rfComputeStruct = struct();
    obj.rfComputeStruct.modelConstants = modelConstants;

    % All solutions
    obj.rfComputeStruct.allSolutions = allMins;
    obj.rfComputeStruct.retinalConePoolingParams = retinalConePoolingParams;

    
    if (modelConstants.simulateCronerKaplanEstimation)
        obj.rfComputeStruct.theSTF = struct('support', stfSupport, 'fitted', theFittedSTF, 'target', theTargetSTF);
    end

    % Convolve the center retinal cone-pooling based retinal RF to get the corresponding visual center RF
    obj.rfComputeStruct.theFittedVisualRFMap = theFittedVisualRF;
    obj.rfComputeStruct.theFittedVisualRFcenterConeMap = conv2(theRetinalRFcenterConeMap, modelConstants.thePSF, 'same');
    obj.rfComputeStruct.theFittedVisualRFsurroundConeMap = conv2(theRetinalRFsurroundConeMap, modelConstants.thePSF, 'same');
    
    obj.rfComputeStruct.theRetinalRFcenterConeMap = theRetinalRFcenterConeMap;
    obj.rfComputeStruct.theRetinalRFsurroundConeMap = theRetinalRFsurroundConeMap;
    
    obj.rfComputeStruct.targetVisualRFfunctionHandle = targetVisualRFfunctionHandle;
    obj.rfComputeStruct.targetVisualRFparamsVector = targetVisualRFparamsVector;
    obj.rfComputeStruct.targetVisualRFMap = targetVisualRFmap;
    obj.rfComputeStruct.targetVisualRFcenterMap = targetVisualRFcenterMap;
    obj.rfComputeStruct.targetVisualRFsurroundMap = targetVisualRFsurroundMap;

    

    % Nested function rfObjective
    function rmse = matchedVisualRFobjective(currentRetinalPoolingParams)
        % Compute the fitted visual RF given the current cone pooling params
        fittedVisualRF = RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(modelConstants, currentRetinalPoolingParams);

        % RF assessment
        if (modelConstants.simulateCronerKaplanEstimation)
            theFittedRFmeasurement = sum(fittedVisualRF,1);
            [~, theFittedSTF] = RetinaToVisualFieldTransformer.spatialTransferFunction(...
                 modelConstants.spatialSupportDegs(:,1),theFittedRFmeasurement);
            % RMSE based on STFs
            fullRMSE = ((theFittedSTF(:) - theTargetSTF(:))).^2;
        else
            theFittedRFmeasurement = fittedVisualRF(:);
            % RMSE based on RFs
            fullRMSE = ((theFittedRFmeasurement(:) - theTargetRFmeasurement(:))).^2;
        end

        % Compute RMSE
        rmse = sqrt(mean(fullRMSE,1));
    end
end



function visualizeRFs(figNo, spatialSupportDegs, RFmap, RFmap2, theTitle)
    figure(figNo); clf;
    profileRF = sum(RFmap,1);
    profileRF2 = sum(RFmap2,1);

    subplot(2,3,1)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RFmap);
    set(gca, 'CLim', [-1 1])
    axis 'image'
    title(theTitle)
    colormap(gray)

    subplot(2,3,2)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RFmap);
    set(gca, 'CLim', 0.01*[-1 1])
    axis 'image'

    subplot(2,3,3)
    plot(spatialSupportDegs(:,1),  profileRF, 'k-');
    axis 'square'

    subplot(2,3,4)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RFmap2);
    set(gca, 'CLim', [-1 1])
    axis 'image'
    title(theTitle)
    colormap(gray)

    subplot(2,3,5)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RFmap2);
    set(gca, 'CLim', 0.01*[-1 1])
    axis 'image'

    subplot(2,3,6)
    plot(spatialSupportDegs(:,1),  profileRF2, 'r-');
    axis 'square'

end


