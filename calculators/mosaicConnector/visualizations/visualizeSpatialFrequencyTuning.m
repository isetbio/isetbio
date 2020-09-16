function visualizeSpatialFrequencyTuning(axesHandle, spatialFrequenciesCPD, theSFtuning, theSFtuningSE, maxSpikeRateModulation, ...
    spikeRateTicks, spatialFrequenciesCPDHR, responseTuningHR, modelParams, modelFitted, targetRGC, LMScontrast, opticsPostFix, PolansSubjectID, exportFig, exportDir, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('synthParams', [], @(x)(isempty(x)||(isstruct(x))));
    p.parse(varargin{:});
    synthParams = p.Results.synthParams;
    
    % Reset figure
    if (isempty(axesHandle))
        % Set plotLab
        figSize = [7 7];
        plotlabOBJ = setupPlotLab(0, figSize);
        hFig = figure(334); clf;
        axesHandle = axes('Position', [0.12*figSize(1) 0.12*figSize(2) 0.85*figSize(1) 0.8*figSize(2)]);
    else
        exportFig = false;
    end
    
    % Plot the error of the mean vertical bars
    for k = 1:numel(spatialFrequenciesCPD)
        line(axesHandle, spatialFrequenciesCPD(k)*[1 1], theSFtuning(k)+theSFtuningSE(k)*[-1 1], 'LineWidth', 2.0, 'Color', [1 .5 0.5]); hold on;
    end
            
    % Plot the model fit
    line(axesHandle, spatialFrequenciesCPDHR, responseTuningHR, 'Color', [1 0 0]);
    
    % Plot the mean points
    scatter(axesHandle, spatialFrequenciesCPD, theSFtuning, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]);
    
    if (~isempty(synthParams)) && (strcmp(modelFitted, 'DifferenceOfGaussians'))
        responseGain = max(theSFtuning);
        visualizeSFTuningOfUnderlyingModel(axesHandle, targetRGC, synthParams, responseGain, spatialFrequenciesCPDHR, max(responseTuningHR));
    end
            
    % Set the axes
    axis(axesHandle, 'square');
    set(axesHandle, 'XScale', 'log', 'XLim', [0.06 105], 'XTick', [0.03 0.1 0.3 1 3 10 30 100], ...
        'YTick', spikeRateTicks, 'YScale', 'linear','YLim', [0 maxSpikeRateModulation]);
    xlabel(axesHandle,'spatial frequency (c/deg)');
    ylabel(axesHandle, 'response modulation');
            
    title(axesHandle, sprintf('RGC #%d, C_{LMS} = <%0.1f, %0.1f, %0.1f>', targetRGC, LMScontrast(1), LMScontrast(2), LMScontrast(3)), ...
        'FontName', 'Source Code Pro', 'FontSize', 24);
    
    % Show the params of the fitted model
    if (~isempty(modelParams))
        if (strcmp(modelFitted, 'DifferenceOfGaussians'))
            theText = sprintf('K_c: %2.0f\nK_s: %2.1f\nR_c: %2.1f arcmin\nR_s: %2.1f arcmin', modelParams.kC, modelParams.kS, modelParams.rC*60, modelParams.rS*60);
        else
            theText = sprintf('K_c: %2.0f\nR_c: %2.1f arcmin', modelParams.kC,  modelParams.rC*60);
        end
        
        text(axesHandle, 0.1, maxSpikeRateModulation*0.85, theText, ...
             'FontSize', 12, 'FontName', 'Source Code Pro', 'BackgroundColor', [1 1 1]);
    end
    
    drawnow;
    
    if (exportFig)
        pdfFileName = sprintf('SF_RGC_%d_LMS_%0.2f_%0.2f_%0.2f_PolansSID_%d_%s',targetRGC, LMScontrast(1), LMScontrast(2), LMScontrast(3), PolansSubjectID, opticsPostFix);
        plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, exportDir);
        setupPlotLab(-1);
    end
end

function plotlabOBJ = setupPlotLab(mode, figSize)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [1 0 0; 1 0 0], ...
                'axesBox', 'off', ...
                'axesTickDir', 'both', ...
                'renderer', 'painters', ...
                'lineMarkerSize', 14, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'figureWidthInches', figSize(1), ...
                'figureHeightInches', figSize(2));
    else
        plotlab.resetAllDefaults();
    end
end 