function psychometricFunctionView(app, mode, varargin)
    p = inputParser;
    p.addParameter('withData',  []);
    p.parse(varargin{:});
    psychometricData = p.Results.withData;
            
    switch (mode)
        case 'initialize'
            initializePsychometricFunctionView(app);
        case 'update'
            updatePsychometricFunctionViewWithNewData(app, psychometricData);
    end
end

function updatePsychometricFunctionViewWithNewData(app, psychometricData)

    % Compute ranges of psychometric functions
    lapse = 0; guess = 0.5;
    log10Contrast = -6:0.01:0;
    log10ContrastMax = app.psychometricFunctionParams.log10ContrastMax;
    log10ContrastMin = app.psychometricFunctionParams.log10ContrastMin;
    minSlope = app.psychometricFunctionParams.slopeMin;
    maxSlope = app.psychometricFunctionParams.slopeMax;
    log10ContrastMean = 0.5*(log10ContrastMin+log10ContrastMax);
    pCorrectMaxSlope = 1 - (lapse - (guess + lapse - 1)*exp(-10.^(maxSlope*(log10Contrast - log10ContrastMean)/20)));
    pCorrectMinSlope = 1 - (lapse - (guess + lapse - 1)*exp(-10.^(minSlope*(log10Contrast - log10ContrastMean)/20)));
    pCorrectMaxSlopeMaxThreshold = 1 - (lapse - (guess + lapse - 1)*exp(-10.^(maxSlope*(log10Contrast - log10ContrastMax)/20)));
    pCorrectMinSlopeMaxThreshold = 1 - (lapse - (guess + lapse - 1)*exp(-10.^(minSlope*(log10Contrast - log10ContrastMax)/20)));
    pCorrectMaxSlopeMinThreshold = 1 - (lapse - (guess + lapse - 1)*exp(-10.^(maxSlope*(log10Contrast - log10ContrastMin)/20)));
    pCorrectMinSlopeMinThreshold = 1 - (lapse - (guess + lapse - 1)*exp(-10.^(minSlope*(log10Contrast - log10ContrastMin)/20)));
    log10CthresholdCoords = [log10Contrast, fliplr(log10Contrast)];
    inBetweenPCorrectMaxSlopeCoords = [pCorrectMaxSlopeMinThreshold, fliplr(pCorrectMaxSlopeMaxThreshold)];
    inBetweenPCorrectMinSlopeCoords = [pCorrectMinSlopeMinThreshold, fliplr(pCorrectMinSlopeMaxThreshold)];
            
    % Plot ranges of psychometric functions as filled plots
    hold(app.psychometricFunctionView, 'off');
    fill(app.psychometricFunctionView, log10CthresholdCoords, inBetweenPCorrectMaxSlopeCoords, ...
        [1 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeAlpha', 0.35, 'EdgeColor', [1 0 0], 'LineWidth', 0.5);
    hold(app.psychometricFunctionView, 'on');
    fill(app.psychometricFunctionView, log10CthresholdCoords, inBetweenPCorrectMinSlopeCoords, ...
        [0.5 0.5 1], 'FaceAlpha', 0.15, 'EdgeAlpha', 0.35, 'EdgeColor', [0 0 1], 'LineWidth', 0.5);
    if (isempty(psychometricData))
        plot(app.psychometricFunctionView, log10Contrast, pCorrectMaxSlope, 'r-', 'LineWidth', 1.5);
        plot(app.psychometricFunctionView, log10Contrast, pCorrectMinSlope, 'b-', 'LineWidth', 1.5);
    end
    
    % Also change the XLims
    %set(app.psychometricFunctionView, 'XLim', [app.psychometricFunctionParams.log10ContrastMin min([0 app.psychometricFunctionParams.log10ContrastMax])]);   
    set(app.psychometricFunctionView, 'XLim', [-6 0]);           
    if (~isempty(psychometricData))
        % Plot the psychometric functions up to this point
        for iSF = 1:numel(psychometricData)
            psychometricDataForThisSF = psychometricData{iSF};
            plot(app.psychometricFunctionView, ...
                 psychometricDataForThisSF.examinedContrastsFit,  psychometricDataForThisSF.pCorrectFit, ...
                 'k-', 'Color', [0.3 0.3 0.3], 'LineWidth', 4);
            plot(app.psychometricFunctionView, ...
                 psychometricDataForThisSF.examinedContrastsFit,  psychometricDataForThisSF.pCorrectFit, ...
                 'k-', 'Color', squeeze(app.csfLineColors(iSF,:)), 'LineWidth', 2);
            scatter(app.psychometricFunctionView, ...
                 psychometricDataForThisSF.examinedContrasts, psychometricDataForThisSF.pCorrect, 14*14, ...
                 'Marker', 'o', 'MarkerEdgeColor', [0.3 0.3 0.3], 'MarkerFaceColor', squeeze(app.csfLineColors(iSF,:)), ...
                 'MarkerFaceAlpha', 0.7, 'LineWidth', 2);
            drawnow;
        end % iSF
    end
    
    
end

function initializePsychometricFunctionView(app)
    set(app.psychometricFunctionView, 'XLim', [-6 0.5], 'YLim', [0.49 1.01], 'XScale', 'linear');
    set(app.psychometricFunctionView, 'XTick', [-6 -5 -4 -3 -2 -1 0], 'YTick', 0.5:0.1:1, 'FontSize', 14);
    xlabel(app.psychometricFunctionView, 'log10(contrast)');
    ylabel(app.psychometricFunctionView, 'pCorrect');
    grid(app.psychometricFunctionView, 'on');
    box(app.psychometricFunctionView, 'on');

    % Do not show the interactions toolbax
    app.psychometricFunctionView.Toolbar.Visible = 'off';

    % Only allow paning
    app.psychometricFunctionView.Interactions = [];
end

