function [patchDogParams, spatialFrequenciesCPDHR, responseTuningHR, meanParams] = ...
    fitDoGmodelToSpatialFrequencyCurve(spatialFrequenciesCPD, responseTuning, plotFitResults, initialParams, maxSpikeRate)

    
    DoGFunction = @(params,sf)(...
        params(1)           * ( pi * (params(2))^2 * exp(-(pi*params(2)*sf).^2) ) - ...
        params(1)*params(3) * ( pi * (params(2)*params(4))^2 * exp(-(pi*params(2)*params(4)*sf).^2) ) );
    %               Kc       Rc     kS/kC       Rs/Rc
    lowerBounds   = [1     0.001   1e-4   2];
    upperBounds   = [Inf   1.0     1e-1   20];

    
    
    rgcsNum = size(responseTuning,1);
    patchDogParams = cell(1, rgcsNum);
    spatialFrequenciesCPDHR = logspace(log10(spatialFrequenciesCPD(1)), log10(spatialFrequenciesCPD(end)), 30);
    responseTuningHR = zeros(rgcsNum, numel(spatialFrequenciesCPDHR));
    
    for iRGC = 1:rgcsNum
        theSFtuning = responseTuning(iRGC,:);
        if (isempty(initialParams))
            initialParams = [max(theSFtuning)   0.05    1/1000       0.5];
        end
        fittedParams(iRGC,:)  = lsqcurvefit(DoGFunction, initialParams,spatialFrequenciesCPD, theSFtuning,lowerBounds,upperBounds);
        patchDogParams{iRGC} = struct(...
            'kC', fittedParams(iRGC,1), ...
            'rC', fittedParams(iRGC,2), ...
            'kS', fittedParams(iRGC,3)*fittedParams(iRGC,1), ...
            'rS', fittedParams(iRGC,4)*fittedParams(iRGC,2) ...
            );
        responseTuningHR(iRGC,:) = DoGFunction(squeeze(fittedParams(iRGC,:)), spatialFrequenciesCPDHR);
        if (plotFitResults)
            figure(334); clf;
            plot(spatialFrequenciesCPD, theSFtuning, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.8 0.8 0.8]); hold on;
            plot(spatialFrequenciesCPDHR, squeeze(responseTuningHR(iRGC,:)), 'r-', 'LineWidth', 1.5);
            set(gca, 'XScale', 'log', 'XLim', [0.1 60], 'YLim', [0 maxSpikeRate]);
            p = patchDogParams{iRGC};
            title(sprintf('kc=%2.0f, ks=%2.0f, rc=%2.4f, rs=%2.4f',p.kC, p.kS, p.rC, p.rS));
            drawnow;
            pause
        end
        
    end
    meanParams = mean(fittedParams,1);
end
