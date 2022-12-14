function retinalRFparamsForTargetVisualRF(obj, indicesOfConesPooledByTheRFcenter, ...
    weightsOfConesPooledByTheRFcenter, targetVisualRFDoGparams)
    
    % Spatial support
    spatialSupportDegs = [obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs(:) obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs(:)];
    [Xdegs,Ydegs] = meshgrid(obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs(:), obj.theVlambdaWeightedPSFData.spatialSupportForRFmapXdegs(:));
    Rdegs2 = Xdegs.^2+Ydegs.^2;

    % Compute the visual RF center and its characteristic radius
    [visualRFcenterCharacteristicRadiusDegs, visualRFcenterConeMap, ...
     visualRFcenterCharacteristicRadiiDegs, visualRFcenterFlatTopExponents, ...
     visualRFcenterXYpos, bestHorizontalResolutionRotationDegs, ...
     anatomicalRFcenterCharacteristicRadiusDegs] = obj.analyzeRFcenter(...
               indicesOfConesPooledByTheRFcenter, ...
               weightsOfConesPooledByTheRFcenter, ...
               spatialSupportDegs);

    if (obj.simulateCronerKaplanEstimation)
        if (~(strcmp(targetVisualRFDoGparams.visualRFmodel, 'gaussian center, gaussian surround')))
            fprintf(2,'WARNINGL When simulating the Croner & Kaplan estimation, the targetVisualRFDoGparams.visualRFmodel is supposed to be ''gaussian center, gaussian surround''.');
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
    modelConstants.thePSF = obj.theVlambdaWeightedPSFData.vLambdaWeightedData;
    modelConstants.surroundConnectableConeTypes = targetVisualRFDoGparams.surroundConnectableConeTypes;
    modelConstants.centerConnectableConeTypes = targetVisualRFDoGparams.centerConnectableConeTypes;
    modelConstants.coneWeightsCompensateForVariationsInConeEfficiency = targetVisualRFDoGparams.coneWeightsCompensateForVariationsInConeEfficiency;
     
    switch (targetVisualRFDoGparams.retinalConePoolingModel)
        case {'arbitrary center cone weights, double exponential surround weights-free', ...
              'arbitrary center cone weights, double exponential surround weights-meanVnVwRatio', ...
              'arbitrary center cone weights, double exponential surround weights-meanRnRwRatio', ...
              'arbitrary center cone weights, double exponential surround from H1 cell with index 1', ...
              'arbitrary center cone weights, double exponential surround from H1 cell with index 2', ...
              'arbitrary center cone weights, double exponential surround from H1 cell with index 3', ...
              'arbitrary center cone weights, double exponential surround from H1 cell with index 4'}
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
            modelConstants.coneCharacteristicRadiusConversionFactor = obj.coneCharacteristicRadiusConversionFactor;
            modelConstants.weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.conePoolingCoefficientsForArbitraryCenterDoubleExpSurround;

            % From the 4 cells in Figure 6 of Packer & Dacey (2002)
            RnarrowToRwideRatios  = [152/515 170/718 115/902 221/1035];
            NWvolumeRatios        = [1.0     0.8     0.3     0.2];

            % Range for RwDegs, based to characteristic radius of cones and # of cones in RF center
            RwDegsInitial    = 4.00 * 1/mean(RnarrowToRwideRatios) * anatomicalRFcenterCharacteristicRadiusDegs * sqrt(2.3) * sqrt(numel(modelConstants.indicesOfCenterCones));
            RwDegsLowerBound = 0.25 * RwDegsInitial;
            RwDegsUpperBound = min([max(spatialSupportDegs(:))*2 15*RwDegsInitial]);

            %                                        Kc   Ks/KcRatio   narrowToWideFieldVolumeRatio  RwideDegs            RnarrowToRwideRatio
            retinalConePoolingParams.names =         {'Kc', 'KsKcRatio', 'VnVwRatio',                'RwDegs',             'RnRwRatio'};
            retinalConePoolingParams.scaling =       {'log', 'log',      'log',                      'linear',                'log'};
            retinalConePoolingParams.initialValues = [1       0.05        mean(NWvolumeRatios)       RwDegsInitial         mean(RnarrowToRwideRatios)];
            retinalConePoolingParams.lowerBounds   = [0.5    1e-3         min(NWvolumeRatios)        RwDegsLowerBound      min(RnarrowToRwideRatios)];
            retinalConePoolingParams.upperBounds   = [100      1e0         max(NWvolumeRatios)        RwDegsUpperBound      max(RnarrowToRwideRatios)];

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

            if (~isempty(strfind(targetVisualRFDoGparams.retinalConePoolingModel,'H1 cell with index')))

                if (~isempty(strfind(targetVisualRFDoGparams.retinalConePoolingModel,'H1 cell with index 1')))
                    % Get values from first H1 cell
                    H1cellIndex = 1;
                elseif (~isempty(strfind(targetVisualRFDoGparams.retinalConePoolingModel,'H1 cell with index 2')))
                    % Get values from second H1 cell
                    H1cellIndex = 2;
                elseif (~isempty(strfind(targetVisualRFDoGparams.retinalConePoolingModel,'H1 cell with index 3')))
                    % Get values from third H1 cell
                    H1cellIndex = 3;
                elseif (~isempty(strfind(targetVisualRFDoGparams.retinalConePoolingModel,'H1 cell with index 4')))
                    % Get values from fourth H1 cell
                    H1cellIndex = 4;
                else
                    error('Unknown model: ''%s''.', targetVisualRFDoGparams.retinalConePoolingModel);
                end

                tolerance = 20/100;
                
                idx = find(ismember(retinalConePoolingParams.names, 'VnVwRatio'));
                measuredValue = NWvolumeRatios(H1cellIndex);
                retinalConePoolingParams.initialValues(idx) = measuredValue;
                retinalConePoolingParams.lowerBounds(idx) = measuredValue / (1.0+tolerance);
                retinalConePoolingParams.upperBounds(idx) = measuredValue * (1.0+tolerance);

                idx = find(ismember(retinalConePoolingParams.names, 'RnRwRatio'));
                measuredValue = RnarrowToRwideRatios(H1cellIndex);
                retinalConePoolingParams.initialValues(idx) = measuredValue;
                retinalConePoolingParams.lowerBounds(idx) = measuredValue / (1.0+tolerance);
                retinalConePoolingParams.upperBounds(idx) = measuredValue * (1.0+tolerance);
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
            error('Unknown retinalConePoolingModel: ''%s''.', targetVisualRFDoGparams.retinalConePoolingModel);
    end


    % Compute the target visual RF map
    [theTargetVisualRFmap, theTargetVisualRFcenterMap, theTargetVisualRFsurroundMap] = ...
        targetVisualRFfunctionHandle(modelConstants, targetVisualRFparamsVector);
   

    % RF assessment type
    modelConstants.simulateCronerKaplanEstimation = obj.simulateCronerKaplanEstimation;
    
    debugRFrotation = true;

    if (modelConstants.simulateCronerKaplanEstimation)

        % TargetSTF match mode: either the params of the DoG model fit to
        % the visual STF or the STF itself
        modelConstants.targetSTFmatchMode = obj.targetSTFmatchMode;

        if (strcmp(modelConstants.targetSTFmatchMode, 'STFDoGparams'))
            multiStartsNumDoGFit = obj.multiStartsNumDoGFit;
        else
            multiStartsNumDoGFit = [];
        end

        % Compute theTargetSTF from targetVisualRFmap
        [~, rotatedTargetVisualRFmap, theTargetRFprofile, theTargetSTF, stfSupportCPD] = ...
            performCronerKaplanSimulation(...
                theTargetVisualRFmap, ...
                visualRFcenterCharacteristicRadiusDegs, ...
                bestHorizontalResolutionRotationDegs, ...
                modelConstants.spatialSupportDegs(:,1), ...
                modelConstants.targetSTFmatchMode, ...
                [], ...
                targetVisualRFDoGparams, ...
                anatomicalRFcenterCharacteristicRadiusDegs, ...
                multiStartsNumDoGFit);

        
        if (debugRFrotation)
            hFig = figure(333);
            set(hFig, 'Position', [10 10 1000 1000]);
            subplot(2,2,1)
            imagesc(theTargetVisualRFmap)
            set(gca, 'CLim', 0.1*max(abs(theTargetVisualRFmap(:)))*[-1 1]);
            axis 'image'
            title('target')
            
            subplot(2,2,2)
            imagesc(rotatedTargetVisualRFmap)
            set(gca, 'CLim', 0.1*max(abs(theTargetVisualRFmap(:)))*[-1 1]);
            axis 'image'
            title('target (rotated)')
            colormap(gray(1024));
        end
    end

    
    % Compute initial visual RF map
    [theInitialFittedVisualRFmap, theInitialRetinalRFcenterConeMap, theInitialRetinalRFsurroundConeMap] = ...
        RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(modelConstants, retinalConePoolingParams.initialValues);
    
    visualizeRFs(53, spatialSupportDegs, theInitialRetinalRFcenterConeMap, theInitialRetinalRFsurroundConeMap, {'initial retinal RF center',  'initial retinal RF surround'});
    visualizeRFs(54, spatialSupportDegs, theTargetVisualRFmap, theInitialFittedVisualRFmap, {'target RF',  'initial fitted RF'});
    

    if (modelConstants.simulateCronerKaplanEstimation)  

        [~, theRotatedInitialFittedVisualRFmap, theInitialRFprofile, theInitialFittedSTF, stfSupportCPD] = performCronerKaplanSimulation(...
                theInitialFittedVisualRFmap, ...
                visualRFcenterCharacteristicRadiusDegs, ...
                bestHorizontalResolutionRotationDegs, ...
                modelConstants.spatialSupportDegs(:,1), ...
                modelConstants.targetSTFmatchMode, ...
                theTargetSTF, ...
                targetVisualRFDoGparams, ...
                anatomicalRFcenterCharacteristicRadiusDegs, ...
                multiStartsNumDoGFit);

        displaySpatialProfileAndSTF(57, spatialSupportDegs, stfSupportCPD, ...
            theTargetRFprofile, theInitialRFprofile, ...
            theTargetSTF, theInitialFittedSTF, ...
            [], [], [], [], 'initial');


        if (debugRFrotation)
            % Add to figure 333
            figure(333)
            subplot(2,2,3)
            imagesc(theInitialFittedVisualRFmap, 0.1*max(abs(theRotatedInitialFittedVisualRFmap(:)))*[-1 1])
            axis 'image'
            title('initial fitted visual map')
            
            subplot(2,2,4)
            imagesc(theRotatedInitialFittedVisualRFmap, 0.1*max(abs(theRotatedInitialFittedVisualRFmap(:)))*[-1 1])
            title('initial fitted visual map (rotated)')
            axis 'image'
        end

    end

    % Optimization optics
    options = optimoptions(...
                'fmincon',...
                'Display','iter',...
                'Algorithm','sqp', ... %'interior-point', ...
                'GradObj', 'off', ...
                'DerivativeCheck', 'off', ...
                'MaxFunEvals', 10^5, ...
                'MaxIter', 256 ...
     );

    
    % Zero the rmseSequence
    rmseSequence = [];

    % Start timer 
    tic

    if (obj.multiStartsNum == 1)
         problem = createOptimProblem('fmincon',...
              'objective', @matchedVisualRFobjective, ...
              'x0', retinalConePoolingParams.initialValues, ...
              'lb', retinalConePoolingParams.lowerBounds, ...
              'ub', retinalConePoolingParams.upperBounds, ...
              'options', options...
              );
        
        retinalConePoolingParams.finalValues = fmincon(problem);

    elseif (obj.multiStartsNum ~= 1)

        options.Display = 'final-detailed';
        % Create problem structure
        problem = createOptimProblem('fmincon',...
              'objective', @matchedVisualRFobjective, ...
              'x0', retinalConePoolingParams.initialValues, ...
              'lb', retinalConePoolingParams.lowerBounds, ...
              'ub', retinalConePoolingParams.upperBounds, ...
              'options', options...
              );

        if (obj.multiStartsNum > 1) && (~isinf(obj.multiStartsNum))
            % Generate multi-start problem
            ms = MultiStart(...
              'Display', 'final', ...
              'StartPointsToRun','bounds', ...  % run only initial points that are feasible with respect to bounds
              'UseParallel', false);
    
            % Run the multi-start solver
            [retinalConePoolingParams.finalValues, ~, ~, ~, allMins] = run(ms, problem, obj.multiStartsNum);
        
        elseif (isinf(obj.multiStartsNum))
            % Global solver
            gs = GlobalSearch(...
                'Display', 'final', ...
                'NumStageOnePoints', 500, ...
                'NumTrialPoints', 2500, ...
                'StartPointsToRun','bounds' ...  % run only initial points that are feasible with respect to bounds
                );

            % Run the global search solver
            retinalConePoolingParams.finalValues = run(gs,problem);
        end
    end

    % Done with fitting. Report time to fit the RVFT model
    fprintf('Fitting RVFT model finished in %2.2f hours\n', toc/60/60);

    % Compute the final fitted visual RF
    [theFittedVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, pooledConeIndicesAndWeights] = ...
            RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(modelConstants, retinalConePoolingParams.finalValues);


    if (modelConstants.simulateCronerKaplanEstimation)

        % Compute final result
        [theFinalRMSEvector, theRotatedFittedVisualRF, theFittedRFprofile, ...
         theFinalFittedSTF, stfSupportCPD, ...
         theFinalFittedSTFsurroundToCenterRcRatio, ...
         theFinalFittedSTFsurroundToCenterIntegratedSensitivityRatio] = performCronerKaplanSimulation(...
                theFittedVisualRF, ...
                visualRFcenterCharacteristicRadiusDegs, ...
                bestHorizontalResolutionRotationDegs, ...
                modelConstants.spatialSupportDegs(:,1), ...
                modelConstants.targetSTFmatchMode, ...
                theTargetSTF, ...
                targetVisualRFDoGparams, ...
                anatomicalRFcenterCharacteristicRadiusDegs, ...
                multiStartsNumDoGFit);

        
        if (obj.multiStartsNum ~= 1)
            plotTitle = 'fitted (single-run)';
        else
            plotTitle = 'fitted (multi-run)';
        end

        displaySpatialProfileAndSTF(158, spatialSupportDegs, stfSupportCPD, ...
                theTargetRFprofile, theFittedRFprofile, ...
                theTargetSTF, theFinalFittedSTF, ...
                targetVisualRFDoGparams.surroundToCenterRcRatio, ...
                targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio, ...
                theFinalFittedSTFsurroundToCenterRcRatio, ...
                theFinalFittedSTFsurroundToCenterIntegratedSensitivityRatio, ...
                plotTitle);
    else
        % Compute final RMSEvector
        theFinalRMSEvector  = ((theFittedVisualRF(:) - theTargetVisualRFmap(:))).^2;
    end

    % Final RMSE
    theFinalRMSE = sqrt(mean(theFinalRMSEvector ,1));


    % Visualize the fitted params
    hFig = figure(1999); clf;
    set(hFig, 'Position', [1000 300 600 500], 'Color', [1 1 1 ]);
    ax = subplot('Position', [0.1 0.1 0.9 0.9]);
    RetinaToVisualFieldTransformer.visualizeRetinalSurroundModelParametersAndRanges(ax, retinalConePoolingParams);


    % Form the rfComputeStruct
    obj.rfComputeStruct = struct();
    obj.rfComputeStruct.modelConstants = modelConstants;

    if (obj.multiStartsNum > 1) && (~isinf(obj.multiStartsNum))
        % All solutions
        obj.rfComputeStruct.allSolutions = allMins;
    end


    obj.rfComputeStruct.retinalConePoolingParams = retinalConePoolingParams;
    obj.rfComputeStruct.pooledConeIndicesAndWeights = pooledConeIndicesAndWeights;
    
    if (modelConstants.simulateCronerKaplanEstimation)
        obj.rfComputeStruct.theSTF = struct('support', stfSupportCPD, 'fitted', theFinalFittedSTF, 'target', theTargetSTF);
        obj.rfComputeStruct.theRotationOfTheFittedVisualRF = bestHorizontalResolutionRotationDegs;
        obj.rfComputeStruct.theRotatedFittedVisualRF = theRotatedFittedVisualRF;
    end

    % Convolve the center retinal cone-pooling based retinal RF to get the corresponding visual center RF
    obj.rfComputeStruct.theFittedVisualRFMap = theFittedVisualRF;
    obj.rfComputeStruct.theFittedVisualRFcenterConeMap = conv2(theRetinalRFcenterConeMap, modelConstants.thePSF, 'same');
    obj.rfComputeStruct.theFittedVisualRFsurroundConeMap = conv2(theRetinalRFsurroundConeMap, modelConstants.thePSF, 'same');
    
    obj.rfComputeStruct.theRetinalRFcenterConeMap = theRetinalRFcenterConeMap;
    obj.rfComputeStruct.theRetinalRFsurroundConeMap = theRetinalRFsurroundConeMap;
    

    obj.rfComputeStruct.targetVisualRFfunctionHandle = targetVisualRFfunctionHandle;
    obj.rfComputeStruct.targetVisualRFparamsVector = targetVisualRFparamsVector;
    obj.rfComputeStruct.targetVisualRFMap = theTargetVisualRFmap;
    obj.rfComputeStruct.targetVisualRFcenterMap = theTargetVisualRFcenterMap;
    obj.rfComputeStruct.targetVisualRFsurroundMap = theTargetVisualRFsurroundMap;

    

    % Nested function rfObjective
    function rmse = matchedVisualRFobjective(currentRetinalPoolingParams)

        % Compute the current visual RF given the current cone pooling params
        currentVisualRF = RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(...
            modelConstants, currentRetinalPoolingParams);

        % RF assessment
        if (modelConstants.simulateCronerKaplanEstimation)

            % Compute current RMSEvector
            theCurrentRMSEvector = performCronerKaplanSimulation(...
                currentVisualRF, ...
                visualRFcenterCharacteristicRadiusDegs, ...
                bestHorizontalResolutionRotationDegs, ...
                modelConstants.spatialSupportDegs(:,1), ...
                modelConstants.targetSTFmatchMode, ...
                theTargetSTF, ...
                targetVisualRFDoGparams, ...
                anatomicalRFcenterCharacteristicRadiusDegs, ...
                multiStartsNumDoGFit);

        else
            % RMSE based on RFs
            theCurrentRMSEvector = ((currentVisualRF(:) - theTargetVisualRFmap(:))).^2;
        end

        % Compute RMSE
        rmse = sqrt(mean(theCurrentRMSEvector,1));

        rmseSequence(numel(rmseSequence)+1) = rmse;
        figure(1555);
        ax = subplot(1,3,3);
        plot(ax,1:numel(rmseSequence), rmseSequence, 'r.', 'MarkerSize', 12, 'LineWidth', 1.0); 
        set(ax, 'YLim', [min([0.01 min(rmseSequence)]), max(rmseSequence)], 'YScale', 'log', 'FontSize', 14);
        xlabel(ax, 'iteration');
        ylabel(ax, 'RMSE');
        title('RMSE sequence');
        drawnow

    end
end


function [theRMSEvector, theRotatedRF, theRFprofile, theVisualSTF, theSpatialFrequencySupport, ...
    theFittedSTFsurroundToCenterRcRatio, theFittedSTFsurroundToCenterIntegratedSensitivityRatio] = ...
        performCronerKaplanSimulation(...
            theVisualRF, visualRFcenterCharacteristicRadiusDegs, bestHorizontalResolutionRotationDegs, sfSupport, ...
            targetSTFmatchMode, theTargetSTF, targetVisualRFDoGparams, ...
            anatomicalRFcenterCharacteristicRadiusDegs, multiStartsNumDoGFit)

    % Rotate theVisualRFmap according to the rotation
    % that maximizes horizontal resolution of the targetVisualRFmap
    theRotatedRF = RetinaToVisualFieldTransformer.bestHorizontalResolutionRFmap(...
        theVisualRF, bestHorizontalResolutionRotationDegs);

    % Integrate along Y to generate the visual RF profile
    theRFprofile = sum(theRotatedRF,1);

    % Compute the visual STF corresponding to the visual RF profile
    [theSpatialFrequencySupport, theVisualSTF] = ...
        RetinaToVisualFieldTransformer.spatialTransferFunction(...
                sfSupport, theRFprofile);

    if (isempty(theTargetSTF))
        theRMSEvector = [];
        theFittedSTFsurroundToCenterRcRatio = [];
        theFittedSTFsurroundToCenterIntegratedSensitivityRatio = [];
    else

        % Compute RMSEvector
        switch  (targetSTFmatchMode) 
            case 'STFcurve'
               % RMSE based on the difference of the STF curves
               theRMSEvector = ((theVisualSTF(:) - theTargetSTF(:))).^2;
    
               theFittedSTFsurroundToCenterRcRatio = [];
               theFittedSTFsurroundToCenterIntegratedSensitivityRatio = [];
    
               figure(1555); clf;
               ax = subplot(1,3,1);
               plot(ax,theSpatialFrequencySupport, theVisualSTF, 'ko-', 'MarkerSize', 12, 'LineWidth', 1.0); hold on;
               plot(ax,theSpatialFrequencySupport, theTargetSTF, 'b-', 'LineWidth', 1.5);
               set(ax, 'XScale', 'log', 'XLim', [0.1 100]);
               grid(ax,'on')
               

            case 'STFDoGparams'

                theVisualSTF = theVisualSTF/max(theVisualSTF(:)) * max(theTargetSTF(:));

                % Fit the final STF with a DoG model
                [theFittedSTFDoGparams, theFittedVisualSTF] = RetinaToVisualFieldTransformer.fitDoGmodelToMeasuredSTF(...
                      theSpatialFrequencySupport, ...
                      theVisualSTF, ...
                      anatomicalRFcenterCharacteristicRadiusDegs, ...
                      multiStartsNumDoGFit);
                      %'rangeForRc', visualRFcenterCharacteristicRadiusDegs*[1/sqrt(2.0) 1.0 sqrt(2.0)]);

                %theRMSEvector = ((theFittedVisualSTF.compositeSTF(:) - theTargetSTF(:))).^2;

                
                theRMSEvector = ((theVisualSTF(:) - theTargetSTF(:))).^2;

                
                theFittedSTFsurroundToCenterRcRatio = theFittedSTFDoGparams.finalValues(3);
                theFittedSTFsurroundToCenterIntegratedSensitivityRatio = theFittedSTFDoGparams.finalValues(2) * (theFittedSTFDoGparams.finalValues(3))^2;
    
                
                figure(1555); clf;
                ax = subplot(1,3,1);
                pMeasured = plot(ax,theSpatialFrequencySupport, theVisualSTF, 'ro', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', 12, 'LineWidth', 1.0); hold on;
                pMeasuredDoGfit = plot(ax,theFittedVisualSTF.sfHiRes, theFittedVisualSTF.compositeSTFHiRes, 'b-', 'LineWidth', 4.0);
                plot(ax,theFittedVisualSTF.sfHiRes, theFittedVisualSTF.compositeSTFHiRes, 'c-', 'LineWidth', 2.0);
                pTarget = plot(ax,theSpatialFrequencySupport, theTargetSTF, 'r-', 'LineWidth', 1.5);
                legend(ax, [pMeasured pMeasuredDoGfit pTarget], {'measured', 'measured (DoG fit)', 'target'}, 'Location', 'SouthWest');

                set(ax, 'XScale', 'log', 'XLim', [0.1 100]);
                grid(ax, 'on')
                title(ax,sprintf('Rs/Rc= (current/target: %2.1f/%2.1f)\nintS/C= (current/target: %2.2f/%2.1f)', ...
                    theFittedSTFsurroundToCenterRcRatio, targetVisualRFDoGparams.surroundToCenterRcRatio, ...
                    theFittedSTFsurroundToCenterIntegratedSensitivityRatio, ...
                    targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio), 'FontSize', 16);
                
                ax = subplot(1,3,2);
                RetinaToVisualFieldTransformer.visualizeRetinalSurroundModelParametersAndRanges(ax, theFittedSTFDoGparams);
                title('fitted DoG model parameters')
                drawnow;

                % Target STF DoG model fit params
                res1 = ((targetVisualRFDoGparams.surroundToCenterRcRatio - ...
                          theFittedSTFsurroundToCenterRcRatio)/targetVisualRFDoGparams.surroundToCenterRcRatio)^2;
    
                res2 = ((targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio -  ...
                          theFittedSTFsurroundToCenterIntegratedSensitivityRatio)/targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio)^2;
    

                %errors = [mean(theRMSEvector) 0.1*res1 0.1*res2]
                theRMSEvector = [mean(theRMSEvector); 0.1*res1; 0.1*res2];
    
            otherwise
                    error('modelConstants.targetSTFmatchMode must be set to either ''STFcurve'' or ''STFDoGparams''.');
        end % switch
    end

end


function visualizeRFs(figNo, spatialSupportDegs, RFmap, RFmap2, theTitles)
    figure(figNo); clf;
    profileRF = sum(RFmap,1);
    profileRF2 = sum(RFmap2,1);

    subplot(2,3,1)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RFmap);
    set(gca, 'CLim', max(abs(RFmap(:)))*[-1 1])
    axis 'image'
    title(theTitles{1})
    colormap(gray(1024))

    subplot(2,3,2)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RFmap);
    set(gca, 'CLim', max(abs(RFmap(:)))*0.1*[-1 1])
    axis 'image'

    subplot(2,3,3)
    plot(spatialSupportDegs(:,1),  profileRF, 'k-');
    axis 'square'

    subplot(2,3,4)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RFmap2);
    set(gca, 'CLim', max(abs(RFmap2(:)))*[-1 1])
    axis 'image'
    title(theTitles{2})
    colormap(gray(1024))

    subplot(2,3,5)
    imagesc(spatialSupportDegs(:,1), spatialSupportDegs(:,2), RFmap2);
    set(gca, 'CLim', max(abs(RFmap2(:)))*0.1*[-1 1])
    axis 'image'

    subplot(2,3,6)
    plot(spatialSupportDegs(:,1),  profileRF2, 'r-');
    axis 'square'
    drawnow;
end

function displaySpatialProfileAndSTF(figNo, spatialSupportDegs, stfSupportCPD, ...
                theTargetRFprofile, theFittedRFprofile, ...
                theTargetSTF,  theFittedSTF, ...
                targetSurroundToCenterRcRatio, ...
                targetSurroundToCenterIntegratedSensitivityRatio, ...
                achievedSurroundToCenterRcRatio, ...
                achievedSurroundToCenterIntegratedSensitivityRatio, ...
                theFittedLegend)

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [100 100 1200 900], 'Color', [1 1 1]);
    legends = {};
    ax = subplot(2,2,1);
    plot(ax, spatialSupportDegs(:,1), theTargetRFprofile, 'k-', 'LineWidth', 1.5);
    legends{1} = 'target';
    hold(ax, 'on');
    if (~isempty(theFittedRFprofile))
        plot(ax, spatialSupportDegs(:,1), theFittedRFprofile, 'b-', 'LineWidth', 1.5);
        legends{numel(legends)+1} = theFittedLegend;
    end
    set(ax, 'YLim', max(abs(theTargetRFprofile(:)))*[-0.5 1.1]);
    xlabel(ax, 'space (degs)');
    ylabel(ax, 'sensitivity');
    set(ax, 'FontSize', 16);
    grid(ax, 'on');
    legend(ax,legends);
    
    ax = subplot(2,2,2);
    plot(ax, stfSupportCPD, theTargetSTF, 'k-', 'LineWidth', 1.5);
    hold(ax, 'on');
    if (~isempty(theFittedSTF))
        plot(ax, stfSupportCPD, theFittedSTF, 'b-', 'LineWidth', 1.5);
    end

    ylabel(ax, 'STF');
    xlabel(ax,'spatial frequency (c/deg)');
    set(ax, 'XScale', 'log', 'YLim', [0 max(theTargetSTF(:))*1.1]);
    set(ax, 'FontSize', 16);
    grid(ax, 'on');

    legend(ax,legends);

    if (~isempty(targetSurroundToCenterRcRatio))
        ax = subplot(2,2,3);
        
        plot([1 10], [1 10], 'k-', 'LineWidth', 1.0); hold(ax, 'on');
        plot(ax, targetSurroundToCenterRcRatio, achievedSurroundToCenterRcRatio, 'ro', ...
            'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5); 
        set(ax, 'XLim', [1 10], 'YLim', [1 10], 'XTick', 1:1:10, 'YTick', 1:1:10);
        axis(ax, 'square');
        grid(ax, 'on'); box(ax, 'off');
        set(ax, 'FontSize', 16);
        xlabel('ax,target');
        ylabel(ax,'achieved');
        title(ax,'Rs/Rc ratio');
    
        ax = subplot(2,2,4);
        plot([0 1], [0 1], 'k-', 'LineWidth', 1.0); hold(ax, 'on');
        plot(ax, targetSurroundToCenterIntegratedSensitivityRatio, achievedSurroundToCenterIntegratedSensitivityRatio, 'ro', ...
            'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5);
        set(ax, 'XLim', [0 1], 'YLim', [0 1],  'XTick', 0:0.1:1.0,'YTick', 0:0.1:1.0);
        axis(ax, 'square');
        grid(ax, 'on'); box(ax, 'off');
        set(ax, 'FontSize', 16);
        xlabel('ax,target');
        title(ax,'S/C int. sensitivity ratio');
    end

    drawnow;
end