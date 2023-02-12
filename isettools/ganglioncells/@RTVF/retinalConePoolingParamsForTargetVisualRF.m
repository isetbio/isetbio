function RFcomputeStruct = retinalConePoolingParamsForTargetVisualRF(obj, ...
                centerConeMajorityType, initialRetinalConePoolingParamsStruct)

  
    RFcomputeStruct = [];

    indicesOfConesPooledByTheRFcenter = obj.targetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter;
    weightsOfConesPooledByTheRFcenter = obj.targetVisualRFDoGparams.weightsOfConesPooledByTheRFcenter;

    switch (centerConeMajorityType)
        case cMosaic.LCONE_ID
            theRFCenterConeMajorityPSF  = obj.spectrallyWeightedPSFData.LconeWeighted;

            spatialPositionDegs = mean(obj.coneMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);
            figureName = sprintf('RF center with %d L-cone(s) at (%2.2f,%2.2f) degs', ...
                numel(indicesOfConesPooledByTheRFcenter), spatialPositionDegs(1), spatialPositionDegs(2));
            summaryFigNo = 2000 + numel(weightsOfConesPooledByTheRFcenter)*10+1;
            
        case cMosaic.MCONE_ID
            theRFCenterConeMajorityPSF  = obj.spectrallyWeightedPSFData.MconeWeighted;

            spatialPositionDegs = mean(obj.coneMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);
            figureName = sprintf('RF center with %d M-cone(s) at (%2.2f,%2.2f) degs', ...
                numel(indicesOfConesPooledByTheRFcenter), spatialPositionDegs(1), spatialPositionDegs(2));
            summaryFigNo = 2000 + numel(weightsOfConesPooledByTheRFcenter)*10+2;

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

    % Zero the rmseSequence
    rmseSequence = [];

    hFigProgress = figure(1000); clf;
    set(hFigProgress, 'Position', [10 10 1580 430], 'Name', figureName);

    % Compute initial visual RF map
    [theInitialRMSE, theInitiaVisualRFmap , theInitialSTFdata] = theObjectiveFunction(retinalConePoolingParams.initialValues);

    
    debugCronerKaplanAnalysis = true;
    if (debugCronerKaplanAnalysis)
        displaySTFdata(100, modelConstants, theInitiaVisualRFmap, theInitialSTFdata, ...
            obj.visualRFcenterRcDegs, ...
            obj.targetVisualRFDoGparams.surroundToCenterRcRatio, ...
            obj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio);
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
    fprintf('Fitting RVFT model finished in %2.2f hours\n', toc/60/60);


     % Compute the fitted visual RF map
    [theFinalRMSE, theFinalVisualRFmap , theFinalSTFdata] = theObjectiveFunction(retinalConePoolingParams.initialValues);


    displaySTFdata(101, modelConstants, theFinalVisualRFmap, theFinalSTFdata, ...
            obj.visualRFcenterRcDegs, ...
            obj.targetVisualRFDoGparams.surroundToCenterRcRatio, ...
            obj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio);

    %  ------- Nested objective function --------
    function [theCurrentRMSE, theCurrentVisualRFmap, theCurrentSTFdata] = theObjectiveFunction(currentRetinalPoolingParamValues)

        % Compute initial visual RF map
        [theCurrentVisualRFmap , theCurrentRetinalRFcenterConeMap, theCurrentRetinalRFsurroundConeMap] = ...
            obj.visualRFfromRetinalConePooling(...
                modelConstants, currentRetinalPoolingParamValues);
  
        % Compute the STF data for the  initial visual RF map
        theCurrentSTFdata = obj.visualRFmapPropertiesFromCronerKaplanAnalysis(theCurrentVisualRFmap);
    
        % Compute RMSE
        RsRcRatioResidual = theCurrentSTFdata.fittedDoGModelRsRcRatio/obj.targetVisualRFDoGparams.surroundToCenterRcRatio - 1;
        SCintSensRatioResidual = theCurrentSTFdata.fittedDoGModelSCIntSensRatio/obj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio - 1;
        rmseWeights = abs([1 1]);    
        rmseWeights = rmseWeights / sum(rmseWeights);
        theCurrentRMSE = sqrt(rmseWeights(1) * RsRcRatioResidual^2 + rmseWeights(2)*SCintSensRatioResidual^2);

        % Update the sequence of RMSEs
        rmseSequence(:,size(rmseSequence,2)+1) = abs([RsRcRatioResidual SCintSensRatioResidual theCurrentRMSE]);


        currentRetinalConePoolingParams = retinalConePoolingParams;
        currentRetinalConePoolingParams.finalValues = currentRetinalPoolingParamValues;

        figure(hFigProgress);
        ax = subplot(2,3,1);
        RTVF.visualizeFittedModelParametersAndRanges(ax, currentRetinalConePoolingParams);

        ax = subplot(2,3,2);
        RTVF.visualizeFittedModelParametersAndRanges(ax, theCurrentSTFdata.fittedDoGModelParams);

        ax = subplot(2,3,3);
        cla(ax);
        iterations = 1:size(rmseSequence,2);
        
        plot(ax,iterations, rmseSequence(1,:), 'r-', 'MarkerSize', 12, 'LineWidth', 1.0); hold(ax, 'on')
        plot(ax,iterations, rmseSequence(2,:), 'b-', 'MarkerSize', 12, 'LineWidth', 1.0); 
        plot(ax,iterations, rmseSequence(3,:), 'k-', 'MarkerSize', 12, 'LineWidth', 1.0); 
        legend(ax, {'Rs/Rc', 'int S/C', 'total'}, 'Location', 'NorthOutside', 'Orientation','horizontal');
        set(ax, 'YLim', max(abs(rmseSequence(:)))*[1e-5 1.1], 'YScale', 'log', 'FontSize', 14);
        xlabel(ax, 'iteration');
        ylabel(ax, 'RMSE');

        if (theCurrentRMSE == min(squeeze(rmseSequence(3,:))))
            spatialSupportDegsX = obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs;
            spatialSupportDegsY = obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs;
            ax = subplot(2,3,4);
            imagesc(ax, spatialSupportDegsX, spatialSupportDegsY, theCurrentRetinalRFcenterConeMap);
            axis(ax, 'image');
            title(sprintf('RF at iteration %d', size(rmseSequence,2)));
            ax = subplot(2,3,5);
            imagesc(ax, spatialSupportDegsX, spatialSupportDegsY, theCurrentRetinalRFsurroundConeMap);
            axis(ax, 'image');
            colormap(gray(1024));
            drawnow;
        end


    end


end






function displaySTFdata(figNo, modelConstants, theVisualRFmap, theVisualSTFdata, ...
    visualRFcenterRcDegs, targetRsRcatio, targetSCIntSensRatio)

    
    figure(figNo); clf;
    ax = subplot(2,3,1);
    imagesc(ax,modelConstants.spatialSupportDegs(:,1), ...
            modelConstants.spatialSupportDegs(:,2), ...
            theVisualRFmap);
    axis(ax, 'image');  axis(ax, 'xy');

    title(sprintf('center max: %2.3f', max(theVisualRFmap(:))))



    ax = subplot(2,3,2);
    plot(ax,theVisualSTFdata.spatialFrequencySupport, theVisualSTFdata.visualSTF, 'ks', 'LineWidth', 1.5);
    hold(ax, 'on');
    plot(ax,theVisualSTFdata.spatialFrequencySupport, theVisualSTFdata.fittedDoGModelToVisualSTF.compositeSTF, 'r-', 'LineWidth', 1.5);
    plot(ax,theVisualSTFdata.spatialFrequencySupport, theVisualSTFdata.fittedDoGModelToVisualSTF.centerSTF, 'r:', 'LineWidth', 1.0);
    plot(ax,theVisualSTFdata.spatialFrequencySupport, theVisualSTFdata.fittedDoGModelToVisualSTF.surroundSTF, 'r--', 'LineWidth', 1.0);
    set(ax, 'XScale', 'log', 'XLim', [0.1 100], 'XTick', [0.1 0.3 1 3 10 30 100]);

    ax = subplot(2,3,4);
    plot(ax,visualRFcenterRcDegs*60, theVisualSTFdata.fittedRcDegs*60 , 'ro', 'LineWidth', 1.5);
    hold(ax, 'on')
    plot([0 20], [0 20], 'k-', 'LineWidth', 1.0)
    xlabel('target Rc (arcmin)');
    ylabel('fitted DoG model Rc (arcmin)')
    set(ax, 'XLim', [0 3], 'YLim', [0 3]);
    axis(ax, 'square');

    ax = subplot(2,3,5);
    plot(ax,targetRsRcatio, theVisualSTFdata.fittedDoGModelRsRcRatio, 'ro', 'LineWidth', 1.5);
    hold(ax, 'on')
    plot([1 20], [1 20], 'k-', 'LineWidth', 1.0)
    xlabel('target Rs/Rc');
    ylabel('fitted DoG model Rs/Rc')
    set(ax, 'XLim', [1 20], 'YLim', [1 20]);
    axis(ax, 'square');

    ax = subplot(2,3,6);
    plot(ax, targetSCIntSensRatio, theVisualSTFdata.fittedDoGModelSCIntSensRatio , 'ro', 'LineWidth', 1.5);
    set(ax, 'XLim', [0 1], 'YLim', [0 1]);
    hold(ax, 'on')
    plot([0 1], [0 1], 'k-', 'LineWidth', 1.0)
    xlabel('target S/C int sens');
    ylabel('fitted S/C int sens')
    axis(ax, 'square');

end