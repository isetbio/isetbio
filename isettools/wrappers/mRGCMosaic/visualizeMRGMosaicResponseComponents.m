function visualizeMRGMosaicResponseComponents(theMRGCMosaic, visualizedMRGCindices, noiseFreeResponse, ...
    noisyResponseInstancesConeNoiseOnly, noisyResponseInstancesIntrinsicMRGCnoiseOnly, generatePDFs)

    %% Visualize the mRGCMosaic
    visualizedWidthDegs  = theMRGCMosaic.sizeDegs(1)*1.01;
    visualizedHeightDegs = theMRGCMosaic.sizeDegs(2)*1.01;
    domainVisualizationLimits(1:2) = theMRGCMosaic.eccentricityDegs(1) + visualizedWidthDegs  * 0.5*[-1 1];
    domainVisualizationLimits(3:4) = theMRGCMosaic.eccentricityDegs(2) + visualizedHeightDegs * 0.5*[-1 1];

    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [1 1 1]);
    ax = subplot('Position', [0.035 0.03 0.96 0.96]);
    theMRGCMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'identifyInputCones', ~true, ...
        'identifyPooledCones', ~true, ...
        'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
        'centerSubregionContourSamples', 24, ...
        'backgroundColor', [1 1 1], ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', struct('x', -10:10, 'y', -10:10), ...
        'plotTitle', '');

    
    if (isempty(visualizedMRGCindices))
        return;
    end
    
    %% Plot cone noise - only responses of visualized mRGCs
    [hFigNoiseFreeActivation, hFigSelectResponses] = visualizeResponsesOfSelectMRGCs(4, theMRGCMosaic, visualizedMRGCindices, ...
        noiseFreeResponse, noisyResponseInstancesConeNoiseOnly, ...
        domainVisualizationLimits, ...
        [1 0.5 0.5], [1 0 0], 'cone noise only');

    if (generatePDFs)
        % Export to PDF
        NicePlot.exportFigToPDF(...
            'noiseFreeMRGCactivation.pdf', ...
            hFigNoiseFreeActivation, 300);
    
        % Export to PDF
        NicePlot.exportFigToPDF(...
            sprintf('coneNoiseOnly_membraneSigma_%2.2f.pdf', theMRGCMosaic.vMembraneGaussianNoiseSigma), ...
            hFigSelectResponses, 300);
    end


    %% Plot intrinsic noise - only responses of visualized mRGCs
    [~, hFigSelectResponses] = visualizeResponsesOfSelectMRGCs(5, theMRGCMosaic, visualizedMRGCindices, ...
        noiseFreeResponse, noisyResponseInstancesIntrinsicMRGCnoiseOnly, ...
        domainVisualizationLimits, ...
        [1 0.5 0.5], [1 0 0], sprintf('intrinsic mRGC noise only (sigma:%2.2f)', theMRGCMosaic.vMembraneGaussianNoiseSigma));

    if (generatePDFs)
        % Export to PDF
        NicePlot.exportFigToPDF(...
            sprintf('intrinsicMRGCNoiseOnly_membraneSigma_%2.2f.pdf', theMRGCMosaic.vMembraneGaussianNoiseSigma), ...
            hFigSelectResponses, 300);
    end

end

function [hFigNoiseFreeActivation, hFigSelectResponses] = visualizeResponsesOfSelectMRGCs(figNo, theMRGCMosaic, visualizedMRGCindices, ...
    noiseFreeResponse, noisyResponseInstances, ...
    domainVisualizationLimits, faceColor, edgeColor, noiseLegend)

    % Visualize the mRGCMosaic noise-free response
    hFigNoiseFreeActivation = figure(3); clf;
    set(hFigNoiseFreeActivation, 'Position', [10 10 1900 1050], 'Color', [1 1 1]);
    ax = subplot('Position', [0.035 0.03 0.96 0.96]);
    theMRGCMosaic.visualize(...
        'figureHandle', hFigNoiseFreeActivation, ...
        'axesHandle', ax, ...
        'activation', noiseFreeResponse, ...
        'activationRange', max(abs(noiseFreeResponse(:)))*[-1 1], ...
        'verticalActivationColorBarInside', true, ...
        'centerSubregionContourSamples', 24, ...
        'backgroundColor', [0 0 0], ...
        'labelRGCsWithIndices', visualizedMRGCindices, ...
        'labeledRGCsColor', [1 0 0], ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', struct('x', -10:10, 'y', -10:10), ...
        'plotTitle', 'noise-free response');


    % Visualize the select mRGC responses
    hFigSelectResponses = figure(figNo); clf;
    set(hFigSelectResponses, 'Position', [10 10 1900 1050], 'Color', [1 1 1]);
    ax = subplot('Position', [0.035 0.05 0.96 0.4]);
    
    % Extract responses of mRGCs around the stimulus y-coord
    timeBin = 1;

    visualizedMRGCXcoords = squeeze(theMRGCMosaic.rgcRFpositionsDegs(visualizedMRGCindices,1));
    visualizedMRGCNoiseFreeResponse = noiseFreeResponse(1, timeBin, visualizedMRGCindices);
    visualizedMRGCNoisyResponseInstances = noisyResponseInstances(:, timeBin, visualizedMRGCindices);
   

    % The 5-95 range of the responses
    visualizedMRGCNoisyResponseRange = prctile(visualizedMRGCNoisyResponseInstances,[5 95], 1);

    % The standard deviations of the responses
    visualizedMRGCNoisyResponseStDev = std(visualizedMRGCNoisyResponseInstances, 0, 1);

    shadedAreaBetweenTwoLines(ax, visualizedMRGCXcoords, ...
        squeeze(visualizedMRGCNoisyResponseRange(1,:,:)), ...
        squeeze(visualizedMRGCNoisyResponseRange(2,:,:)), ...
        faceColor, edgeColor, 0.5, 1.5, '-');

    hold(ax, 'on');
    plot(ax,visualizedMRGCXcoords, squeeze(visualizedMRGCNoiseFreeResponse), 'k.', ...
         'MarkerSize', 18, 'LineWidth', 1.5);

    hold(ax, 'off')
    set(ax, 'XLim', domainVisualizationLimits(1:2), ...
            'YLim', [-1.0 1.0], 'FontSize', 16, ...
            'Color', 'none');
    grid(ax, 'on');
    legend(ax, {noiseLegend, 'noise-free'}, 'Location', 'NorthOutside', 'NumColumns', 2, 'Orientation','horizontal');
    xlabel('space (degs)')
    ylabel('mRGC response')


    ax = subplot('Position', [0.035 0.57 0.96 0.4]);
    plot(ax,visualizedMRGCXcoords, squeeze(visualizedMRGCNoisyResponseStDev), ...
        '-', 'Color', edgeColor, 'LineWidth', 1.5);
    stdRange = [0 0.3];
    set(ax, 'XLim', domainVisualizationLimits(1:2), ...
        'YLim', stdRange, 'FontSize', 16, ...
        'Color', 'none');

    grid(ax, 'on');
    legend(ax, noiseLegend, 'Location', 'NorthOutside', 'NumColumns', 2, 'Orientation','horizontal');
    xlabel('space (degs)')
    ylabel('mRGC response st. dev.')

end

function p = shadedAreaBetweenTwoLines(ax,x, y1, y2, ...
     faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)

    x = reshape(x, [1 numel(x)]);
    y1 = reshape(y1, [1 numel(x)]);
    y2 = reshape(y2, [1 numel(x)]);


    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    p = patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, 'edgeAlpha', 1.0,...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end