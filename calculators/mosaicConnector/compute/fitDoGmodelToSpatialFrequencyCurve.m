function [patchDogParams, spatialFrequenciesCPDHR, responseTuningHR, meanParams] = ...
    fitDoGmodelToSpatialFrequencyCurve(spatialFrequenciesCPD, responseTuning, responseTuningSE, initialParams, ...
    maxSpikeRateModulation, visualizeIndividualFits, LMScontrast, opticsPostFix, PolansSubjectID, exportFig, exportDir)

    % The DOF model function handle
    DoGFunction = @(params,sf)(...
        params(1)           * ( pi * (params(2))^2 * exp(-(pi*params(2)*sf).^2) ) - ...
        params(1)*params(3) * ( pi * (params(2)*params(4))^2 * exp(-(pi*params(2)*params(4)*sf).^2) ) );
    
    % Upper and lower values of DoG params
    %               Kc       Rc     kS/kC       Rs/Rc
    lowerBounds   = [1     0.001   1e-4   2];
    upperBounds   = [Inf   1.0     1e-1   20];
    
    % Fitting options
    oldoptions = optimoptions('lsqcurvefit');
    options = optimoptions(oldoptions,'MaxFunctionEvaluations',5000, 'FunctionTolerance', 1e-8, 'MaxIterations', 10000);
    
    % Preallocate memory
    rgcsNum = size(responseTuning,1);
    patchDogParams = cell(1, rgcsNum);
    fittedParams = zeros(rgcsNum, length(lowerBounds));
    
    % High-resolution spatial frequency axis
    spatialFrequenciesCPDHR = logspace(log10(spatialFrequenciesCPD(1)), log10(spatialFrequenciesCPD(end)), 100);
    responseTuningHR = zeros(rgcsNum, numel(spatialFrequenciesCPDHR));
    
    % Fit spatial frequency curves for each RGC
    for iRGC = 1:rgcsNum
        % Retrieve data to be fitted
        theSFtuning = responseTuning(iRGC,:);
        theSFtuningSE = responseTuningSE(iRGC,:);
        
        % Set initial params if not set
        if (isempty(initialParams))
            initialParams = [max(theSFtuning)   0.05    1/1000       0.5];
        end
        
        % Fit the model to the data
        fittedParams(iRGC,:)  = lsqcurvefit(DoGFunction, initialParams,spatialFrequenciesCPD, theSFtuning,lowerBounds,upperBounds, options);
        patchDogParams{iRGC} = struct(...
            'kC', fittedParams(iRGC,1), ...
            'rC', fittedParams(iRGC,2), ...
            'kS', fittedParams(iRGC,3)*fittedParams(iRGC,1), ...
            'rS', fittedParams(iRGC,4)*fittedParams(iRGC,2) ...
            );
        
        % Compute model fit on a high-res spatial frequency axis
        responseTuningHR(iRGC,:) = DoGFunction(squeeze(fittedParams(iRGC,:)), spatialFrequenciesCPDHR);
        
        % Visualize model fit and data
        if (visualizeIndividualFits)
            visualizeSpatialFrequencyTuning([], spatialFrequenciesCPD, theSFtuning, theSFtuningSE, maxSpikeRateModulation, ...
                spatialFrequenciesCPDHR, squeeze(responseTuningHR(iRGC,:)), patchDogParams{iRGC}, ...
                iRGC, LMScontrast, opticsPostFix, PolansSubjectID, exportFig, exportDir);
        end
    end
    
    % Mean model params for the RGC patch
    meanParams = mean(fittedParams,1);
end
