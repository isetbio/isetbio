function [patchDogParams, spatialFrequenciesCPDHR, responseTuningHR, meanParams] = ...
    fitDoGmodelToSpatialFrequencyCurve(spatialFrequenciesCPD, responseTuning, responseTuningSE, initialParams, ...
    maxSpikeRateModulation, visualizeIndividualFits, exportFig, LMScontrast)

    
    DoGFunction = @(params,sf)(...
        params(1)           * ( pi * (params(2))^2 * exp(-(pi*params(2)*sf).^2) ) - ...
        params(1)*params(3) * ( pi * (params(2)*params(4))^2 * exp(-(pi*params(2)*params(4)*sf).^2) ) );
    %               Kc       Rc     kS/kC       Rs/Rc
    lowerBounds   = [1     0.001   1e-4   2];
    upperBounds   = [Inf   1.0     1e-1   20];
    oldoptions = optimoptions('lsqcurvefit');
    options = optimoptions(oldoptions,'MaxFunctionEvaluations',5000, 'FunctionTolerance', 1e-8, 'MaxIterations', 10000);
    
    
    rgcsNum = size(responseTuning,1);
    patchDogParams = cell(1, rgcsNum);
    spatialFrequenciesCPDHR = logspace(log10(spatialFrequenciesCPD(1)), log10(spatialFrequenciesCPD(end)), 100);
    responseTuningHR = zeros(rgcsNum, numel(spatialFrequenciesCPDHR));
    
    if (visualizeIndividualFits)
        plotlabOBJ = setupPlotLab(0);
    end
    
    for iRGC = 1:rgcsNum
        theSFtuning = responseTuning(iRGC,:);
        theSFtuningSE = responseTuningSE(iRGC,:);
        if (isempty(initialParams))
            initialParams = [max(theSFtuning)   0.05    1/1000       0.5];
        end
        fittedParams(iRGC,:)  = lsqcurvefit(DoGFunction, initialParams,spatialFrequenciesCPD, theSFtuning,lowerBounds,upperBounds, options);
        patchDogParams{iRGC} = struct(...
            'kC', fittedParams(iRGC,1), ...
            'rC', fittedParams(iRGC,2), ...
            'kS', fittedParams(iRGC,3)*fittedParams(iRGC,1), ...
            'rS', fittedParams(iRGC,4)*fittedParams(iRGC,2) ...
            );
        responseTuningHR(iRGC,:) = DoGFunction(squeeze(fittedParams(iRGC,:)), spatialFrequenciesCPDHR);
        if (visualizeIndividualFits)
            hFig = figure(334); clf;
            for k = 1:numel(spatialFrequenciesCPD)
                line(spatialFrequenciesCPD(k)*[1 1], theSFtuning(k)+theSFtuningSE(k)*[-1 1], 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]); hold on;
            end
            scatter(spatialFrequenciesCPD, theSFtuning, 'MarkerEdgeColor', [0.3 0.3 0.3], 'MarkerFaceColor', [0.8 0.8 0.8]);
            line(spatialFrequenciesCPDHR, squeeze(responseTuningHR(iRGC,:)), 'Color', [1 0 0], 'LineWidth', 3);
            set(gca, 'XScale', 'log', 'XLim', [0.03 101], 'XTick', [0.03 0.1 0.3 1 3 10 30], 'YTick', [0:5:200], 'YScale', 'linear','YLim', [0 maxSpikeRateModulation]);
            p = patchDogParams{iRGC};
            title(sprintf('K_c: %2.0f, K_s: %2.0f, R_c: %2.1f arcmin, R_s: %2.1f arcmin',p.kC, p.kS, p.rC*60, p.rS*60));
            xlabel('spatial frequency (c/deg)');
            ylabel('response modulation');
            
            drawnow;
            if (exportFig)
                plotlabOBJ.exportFig(hFig, 'pdf', sprintf('SF_RGC_%d_LMS_%0.2f_%0.2f_%0.2f',iRGC,LMScontrast(1), LMScontrast(2), LMScontrast(3)), pwd());
            end
        end
    end
    
    meanParams = mean(fittedParams,1);
    
    if (visualizeIndividualFits)
        setupPlotLab(-1);
    end

end

function plotlabOBJ = setupPlotLab(mode)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [0 0 0; 1 0 0], ...
                'axesBox', 'off', ...
                'axesTickDir', 'both', ...
                'renderer', 'painters', ...
                'lineMarkerSize', 14, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'figureWidthInches', 7, ...
                'figureHeightInches', 7);
    else
        plotlab.resetAllDefaults();
    end
end 
