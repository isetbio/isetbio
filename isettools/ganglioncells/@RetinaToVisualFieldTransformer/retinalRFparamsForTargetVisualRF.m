function theRFcomputeStruct = retinalRFparamsForTargetVisualRF(obj, indicesOfConesPooledByTheRFcenter, ...
    weightsOfConesPooledByTheRFcenter, targetVisualRFDoGparams, ...
    centerConeType, initialRetinalConePoolingParamsStruct)

    switch (centerConeType)
        case cMosaic.LCONE_ID
            theEmployedPSFData = obj.theSpectrallyWeightedPSFData.LconeWeighted;

            spatialPositionDegs = mean(obj.theConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);
            figureName = sprintf('RF center with %d L-cone(s) at (%2.2f,%2.2f) degs', ...
                numel(indicesOfConesPooledByTheRFcenter), spatialPositionDegs(1), spatialPositionDegs(2));
            summaryFigNo = 2000 + numel(weightsOfConesPooledByTheRFcenter)*10+1;
            
        case cMosaic.MCONE_ID
            theEmployedPSFData = obj.theSpectrallyWeightedPSFData.MconeWeighted;

            spatialPositionDegs = mean(obj.theConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);
            figureName = sprintf('RF center with %d M-cone(s) at (%2.2f,%2.2f) degs', ...
                numel(indicesOfConesPooledByTheRFcenter), spatialPositionDegs(1), spatialPositionDegs(2));
            summaryFigNo = 2000 + numel(weightsOfConesPooledByTheRFcenter)*10+2;

        otherwise
            error('Not L or M cone: %d', centerConeType);
    end


    % Spatial support
    spatialSupportDegs = [...
        obj.theSpectrallyWeightedPSFData.spatialSupportForRFmapXdegs(:) ...
        obj.theSpectrallyWeightedPSFData.spatialSupportForRFmapXdegs(:)];

    [Xdegs,Ydegs] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));
    Rdegs2 = Xdegs.^2+Ydegs.^2;

    % Compute the visual RF center and its characteristic radius
    [visualRFcenterCharacteristicRadiusDegs, visualRFcenterConeMap, ...
     visualRFcenterCharacteristicRadiiDegs, visualRFcenterFlatTopExponents, ...
     visualRFcenterXYpos, bestHorizontalResolutionRotationDegs, ...
     anatomicalRFcenterCharacteristicRadiusDegs] = obj.analyzeRFcenter(...
               indicesOfConesPooledByTheRFcenter, ...
               weightsOfConesPooledByTheRFcenter, ...
               spatialSupportDegs, ...
               theEmployedPSFData);

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
    modelConstants.thePSF = theEmployedPSFData;
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
            RwDegsLowerBound = max([0.02 0.01*RwDegsInitial]);
            RwDegsUpperBound = min([max(spatialSupportDegs(:))*2 4*RwDegsInitial]);

            %                                        Kc   Ks/KcRatio   narrowToWideFieldVolumeRatio  RwideDegs            RnarrowToRwideRatio
            retinalConePoolingParams.names =         {'Kc', 'KsKcRatio', 'VnVwRatio',                'RwDegs',             'RnRwRatio'};
            retinalConePoolingParams.scaling =       {'log', 'log',      'log',                      'linear',                'log'};
            retinalConePoolingParams.initialValues = [1.       0.06        mean(NWvolumeRatios)       RwDegsInitial         mean(RnarrowToRwideRatios)];
            retinalConePoolingParams.lowerBounds   = [0.5    0.005        min(NWvolumeRatios)        RwDegsLowerBound      min(RnarrowToRwideRatios)];
            retinalConePoolingParams.upperBounds   = [2      1e0         max(NWvolumeRatios)        RwDegsUpperBound      max(RnarrowToRwideRatios)];

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

                tolerance = 0.1;
                
                idx = find(ismember(retinalConePoolingParams.names, 'VnVwRatio'));
                measuredValue = NWvolumeRatios(H1cellIndex);
                retinalConePoolingParams.initialValues(idx) = measuredValue;
                retinalConePoolingParams.lowerBounds(idx) = measuredValue - tolerance;
                retinalConePoolingParams.upperBounds(idx) = measuredValue + tolerance;

                idx = find(ismember(retinalConePoolingParams.names, 'RnRwRatio'));
                measuredValue = RnarrowToRwideRatios(H1cellIndex);
                retinalConePoolingParams.initialValues(idx) = measuredValue;
                retinalConePoolingParams.lowerBounds(idx) = measuredValue - tolerance;
                retinalConePoolingParams.upperBounds(idx) = measuredValue + tolerance;
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

        targetRsRcRatio = targetVisualRFDoGparams.surroundToCenterRcRatio;
        targetIntSensSCRatio = targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;

        if (strcmp(modelConstants.targetSTFmatchMode, 'STFDoGparams'))
            multiStartsNumDoGFit = obj.multiStartsNumDoGFit;
        else
            multiStartsNumDoGFit = [];
        end

        % Compute theTargetSTF from targetVisualRFmap
        [~, rotatedTargetVisualRFmap, theTargetRFprofile, theTargetSTF, stfSupportCPD] = ...
             RetinaToVisualFieldTransformer.performCronerKaplanSimulation(...
                theTargetVisualRFmap, ...
                [], [], ...
                retinalConePoolingParams.initialValues, ...
                retinalConePoolingParams, ...
                bestHorizontalResolutionRotationDegs, ...
                modelConstants.spatialSupportDegs(:,1), ...
                modelConstants.targetSTFmatchMode, ...
                [], ...
                targetRsRcRatio, targetIntSensSCRatio, ...
                visualRFcenterCharacteristicRadiusDegs, ...
                multiStartsNumDoGFit, ...
                1555, '');

        
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

        [~, theRotatedInitialFittedVisualRFmap, theInitialRFprofile, theInitialFittedSTF, stfSupportCPD] =  RetinaToVisualFieldTransformer.performCronerKaplanSimulation(...
                theInitialFittedVisualRFmap, ...
                theInitialRetinalRFcenterConeMap, theInitialRetinalRFsurroundConeMap, ...
                retinalConePoolingParams.initialValues, ...
                retinalConePoolingParams, ...
                bestHorizontalResolutionRotationDegs, ...
                modelConstants.spatialSupportDegs(:,1), ...
                modelConstants.targetSTFmatchMode, ...
                theTargetSTF, ...
                targetRsRcRatio, targetIntSensSCRatio, ...
                visualRFcenterCharacteristicRadiusDegs, ...
                multiStartsNumDoGFit, ...
                1555, figureName);

        displaySpatialProfileAndSTF(57, '', [], spatialSupportDegs, stfSupportCPD, ...
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
                'Algorithm','interior-point', ... %'interior-point', 'sqp'...
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
         theFinalFittedSTFsurroundToCenterIntegratedSensitivityRatio, ...
         ratioWeights] =  RetinaToVisualFieldTransformer.performCronerKaplanSimulation(...
                theFittedVisualRF, ...
                theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, ...
                retinalConePoolingParams.finalValues, ...
                retinalConePoolingParams, ...
                bestHorizontalResolutionRotationDegs, ...
                modelConstants.spatialSupportDegs(:,1), ...
                modelConstants.targetSTFmatchMode, ...
                theTargetSTF, ...
                targetRsRcRatio, targetIntSensSCRatio, ...
                visualRFcenterCharacteristicRadiusDegs, ...
                multiStartsNumDoGFit, ...
                1556, figureName);

        
        if (obj.multiStartsNum ~= 1)
            plotTitle = 'fitted (multi-run)';
        else
            plotTitle = 'fitted (single-run)';
        end

        displaySpatialProfileAndSTF(summaryFigNo, figureName, [1300 400], spatialSupportDegs, stfSupportCPD, ...
                theTargetRFprofile, theFittedRFprofile, ...
                theTargetSTF, theFinalFittedSTF, ...
                targetRsRcRatio, targetIntSensSCRatio, ...
                theFinalFittedSTFsurroundToCenterRcRatio, ...
                theFinalFittedSTFsurroundToCenterIntegratedSensitivityRatio, ...
                plotTitle);
    else
        % Compute final RMSEvector
        theFinalRMSEvector  = ((theFittedVisualRF(:) - theTargetVisualRFmap(:))).^2;
    end

  
    % Visualize the fitted params
    hFig = figure(1999); clf;
    set(hFig, 'Position', [1000 300 600 500], 'Color', [1 1 1 ]);
    ax = subplot('Position', [0.1 0.1 0.9 0.9]);
    RetinaToVisualFieldTransformer.visualizeRetinalSurroundModelParametersAndRanges(ax, retinalConePoolingParams);


    % Form theRFcomputeStruct
    theRFcomputeStruct = struct();
    theRFcomputeStruct.modelConstants = modelConstants;

    if (obj.multiStartsNum > 1) && (~isinf(obj.multiStartsNum))
        % All solutions
        theRFcomputeStruct.allSolutions = allMins;
    end


    theRFcomputeStruct.retinalConePoolingParams = retinalConePoolingParams;
    theRFcomputeStruct.pooledConeIndicesAndWeights = pooledConeIndicesAndWeights;
    
    if (modelConstants.simulateCronerKaplanEstimation)
        theRFcomputeStruct.theSTF = struct(...
            'support', stfSupportCPD, ...
            'fitted', theFinalFittedSTF, ...
            'target', theTargetSTF, ...
            'fittedRsRcRatio', theFinalFittedSTFsurroundToCenterRcRatio, ...
            'fittedSCIntSensRatio', theFinalFittedSTFsurroundToCenterIntegratedSensitivityRatio);
        theRFcomputeStruct.theRotationOfTheFittedVisualRF = bestHorizontalResolutionRotationDegs;
        theRFcomputeStruct.theRotatedFittedVisualRF = theRotatedFittedVisualRF;
        theRFcomputeStruct.theFinalRMSEvector = theFinalRMSEvector;
        theRFcomputeStruct.theFinalRMSEvectorRatioWeights = ratioWeights;
    end

    % Convolve the center retinal cone-pooling based retinal RF to get the corresponding visual center RF
    theRFcomputeStruct.theFittedVisualRFMap = theFittedVisualRF;
    theRFcomputeStruct.theFittedVisualRFcenterConeMap = conv2(theRetinalRFcenterConeMap, modelConstants.thePSF, 'same');
    theRFcomputeStruct.theFittedVisualRFsurroundConeMap = conv2(theRetinalRFsurroundConeMap, modelConstants.thePSF, 'same');
    
    theRFcomputeStruct.theRetinalRFcenterConeMap = theRetinalRFcenterConeMap;
    theRFcomputeStruct.theRetinalRFsurroundConeMap = theRetinalRFsurroundConeMap;
    

    theRFcomputeStruct.targetVisualRFfunctionHandle = targetVisualRFfunctionHandle;
    theRFcomputeStruct.targetVisualRFparamsVector = targetVisualRFparamsVector;
    theRFcomputeStruct.targetVisualRFMap = theTargetVisualRFmap;
    theRFcomputeStruct.targetVisualRFcenterMap = theTargetVisualRFcenterMap;
    theRFcomputeStruct.targetVisualRFsurroundMap = theTargetVisualRFsurroundMap;

    

    % Nested function rfObjective
    function theCurrentRMSE = matchedVisualRFobjective(currentRetinalPoolingParams)

        % Compute the current visual RF given the current cone pooling params
        [currentVisualRF, currentRetinalRFcenterConeMap, currentRetinalRFsurroundConeMap] = RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(...
            modelConstants, currentRetinalPoolingParams);

        % RF assessment
        if (modelConstants.simulateCronerKaplanEstimation)

            % Compute current RMSEvector
            theCurrentRMSEvector =  RetinaToVisualFieldTransformer.performCronerKaplanSimulation(...
                currentVisualRF, ...
                currentRetinalRFcenterConeMap, currentRetinalRFsurroundConeMap, ...
                currentRetinalPoolingParams, ...
                retinalConePoolingParams, ...
                bestHorizontalResolutionRotationDegs, ...
                modelConstants.spatialSupportDegs(:,1), ...
                modelConstants.targetSTFmatchMode, ...
                theTargetSTF, ...
                targetRsRcRatio, targetIntSensSCRatio, ...
                visualRFcenterCharacteristicRadiusDegs, ...
                multiStartsNumDoGFit, ...
                1555, figureName);

        else
            % RMSE based on RFs
            theCurrentRMSEvector = ((currentVisualRF(:) - theTargetVisualRFmap(:))).^2;
        end

        % Compute RMSE
        theCurrentRMSE = sqrt(mean(theCurrentRMSEvector,1));
       
        % Update the sequence of RMSEs
        rmseSequence(:,size(rmseSequence,2)+1) = theCurrentRMSEvector;


        figure(1555);
        ax = subplot(2,4, 4);
        iterations = 1:size(rmseSequence,2);
        plot(ax,iterations, rmseSequence(1,:), 'r-', 'MarkerSize', 12, 'LineWidth', 1.0); hold on
        plot(ax,iterations, rmseSequence(2,:), 'b-', 'MarkerSize', 12, 'LineWidth', 1.0); 
        plot(ax,iterations, rmseSequence(3,:), 'k-', 'MarkerSize', 12, 'LineWidth', 1.0); 
        legend(ax, {'shape', 'Rs/Rc', 'int S/C'}, 'Location', 'NorthOutside', 'Orientation','horizontal');
        set(ax, 'YLim', [min([0.01 min(rmseSequence(:))]), max(rmseSequence(:))], 'YScale', 'log', 'FontSize', 14);
        xlabel(ax, 'iteration');
        ylabel(ax, 'RMSE');
        
        drawnow

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

function displaySpatialProfileAndSTF(figNo, figureName, figXYpos, spatialSupportDegs, stfSupportCPD, ...
                theTargetRFprofile, theFittedRFprofile, ...
                theTargetSTF,  theFittedSTF, ...
                targetSurroundToCenterRcRatio, ...
                targetSurroundToCenterIntegratedSensitivityRatio, ...
                achievedSurroundToCenterRcRatio, ...
                achievedSurroundToCenterIntegratedSensitivityRatio, ...
                theFittedLegend)

    [~, idx] = max(theTargetRFprofile(:));
    peakPositionTarget = spatialSupportDegs(idx,1);
                                                                                                                                                                                                        
    [~, idx] = max(theFittedRFprofile(:));
    peakPositionFitted = spatialSupportDegs(idx,1);

    hFig = figure(figNo); clf;
    if (isempty(figXYpos))
        figXYpos = [100 100];
    end

    set(hFig, 'Position', [figXYpos(1) figXYpos(2) 1200 900], 'Color', [1 1 1], 'Name', figureName);
    legends = {};
    ax = subplot(2,2,1);
    plot(ax, spatialSupportDegs(:,1)-peakPositionTarget, theTargetRFprofile, 'k-', 'LineWidth', 2);
    legends{1} = 'target';
    hold(ax, 'on');
    if (~isempty(theFittedRFprofile))
        plot(ax, spatialSupportDegs(:,1)-peakPositionFitted , theFittedRFprofile, 'r-', 'LineWidth', 2);
        legends{numel(legends)+1} = theFittedLegend;
    end
    set(ax, 'YLim', max(abs(theTargetRFprofile(:)))*[-0.5 1.1], 'XLim', 0.3*[-1 1], 'XTick', -1:0.1:1);
    xlabel(ax, 'space (degs)');
    ylabel(ax, 'sensitivity');
    set(ax, 'FontSize', 16);
    grid(ax, 'on');
    legend(ax,legends);
    
    ax = subplot(2,2,2);
    plot(ax, stfSupportCPD, theTargetSTF, 'k-', 'LineWidth', 1.5);
    hold(ax, 'on');
    if (~isempty(theFittedSTF))
        plot(ax, stfSupportCPD, theFittedSTF, 'r-', 'LineWidth', 1.5);
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