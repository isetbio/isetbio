function hFig = visualizeHorizontalConeActivationProfiles(obj, theConeMosaicResponse, coneTypesToVisualize, ...
    maxResponse, visualizedResponseScalingDegs, horizontalSliceYcoordDegs, horizontalSliceWidthDegs)

    % Visualized L-cone indices and x-coords
    if (ismember(cMosaic.LCONE_ID, coneTypesToVisualize))
        [theVisualizedLConeIndices, theROI] = ...
            obj.coneIndicesOfCertainTypeWithinHorizontalSlice(horizontalSliceYcoordDegs, horizontalSliceWidthDegs, cMosaic.LCONE_ID);
    else
        theVisualizedLConeIndices = [];
    end

    % Visualized M-cone indices and x-coords
    if (ismember(cMosaic.MCONE_ID, coneTypesToVisualize))
        theVisualizedMConeIndices = ...
            obj.coneIndicesOfCertainTypeWithinHorizontalSlice(horizontalSliceYcoordDegs, horizontalSliceWidthDegs, cMosaic.MCONE_ID);
    else
        theVisualizedMConeIndices = [];
    end


    % Visualized S-cone indices and x-coords
    if (ismember(cMosaic.SCONE_ID, coneTypesToVisualize))
        theVisualizedSConeIndices = ...
            obj.coneIndicesOfCertainTypeWithinHorizontalSlice(horizontalSliceYcoordDegs, horizontalSliceWidthDegs, cMosaic.SCONE_ID);
    else
        theVisualizedSConeIndices = [];
    end

    visualizationLimits = [...
        obj.eccentricityDegs(1) - obj.sizeDegs(1)*0.51 ...
        obj.eccentricityDegs(1) + obj.sizeDegs(1)*0.51 ...
        obj.eccentricityDegs(2) - obj.sizeDegs(2)*0.51 ...
        obj.eccentricityDegs(2) + obj.sizeDegs(2)*0.51 ...
        ];
    visualizationTicks = struct('x', -20:0.5:20, 'y', -20:0.5:20);

    hFig = figure(1);clf;
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [0 0 0]);
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);

    % Plot the cone mosaic response
    obj.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'activation', theConeMosaicResponse, ...
        'verticalActivationColorBarInside', true, ...
        'backgroundColor', [0.2 0.2 0.2], ...
        'domainVisualizationLimits', visualizationLimits, ...
        'domainVisualizationTicks', visualizationTicks);

    hold(ax, 'on');

    % ROI over which we visualize cone responses
    theROI.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'xLims', visualizationLimits(1:2), ...
        'yLims', visualizationLimits(3:4), ...
        'fillColor', [0.1 1.0 0.1 0.6]);

    % Y-offset
    yOffset = horizontalSliceWidthDegs*0.5+0.1;
    
    % Render visualization background
    renderVisualizationBackground(ax, obj.eccentricityDegs, obj.sizeDegs, ...
        horizontalSliceYcoordDegs, yOffset, visualizedResponseScalingDegs);

    % Render L-cone responses
    if (~isempty(theVisualizedLConeIndices))
        renderConeResponses(ax, squeeze(obj.coneRFpositionsDegs(theVisualizedLConeIndices,1)), theVisualizedLConeIndices, ...
            theConeMosaicResponse, horizontalSliceYcoordDegs+yOffset, visualizedResponseScalingDegs, maxResponse, [1.0 0.2 0.4]);
    end

    % Render M-cone responses
    if (~isempty(theVisualizedMConeIndices))
        renderConeResponses(ax, squeeze(obj.coneRFpositionsDegs(theVisualizedMConeIndices,1)), theVisualizedMConeIndices, ...
            theConeMosaicResponse, horizontalSliceYcoordDegs+yOffset, visualizedResponseScalingDegs, maxResponse, [0.1 1.0 0.6]);
    end


    % Render S-cone responses
    if (~isempty(theVisualizedSConeIndices))
        renderConeResponses(ax, squeeze(obj.coneRFpositionsDegs(theVisualizedSConeIndices,1)), theVisualizedSConeIndices, ...
            theConeMosaicResponse, horizontalSliceYcoordDegs+yOffset, visualizedResponseScalingDegs, maxResponse, [0.6 0.1 1.0]);
    end

    % Finalize plot
    set(hFig, 'Color', [0 0 0]);
    set(ax, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'LineWidth', 1.5, 'FontSize', 30);
    grid(ax, 'off')

end


function renderVisualizationBackground(ax, mosaicEccDegs, mosaicSizeDegs, horizontalSliceYcoordDegs, yOffsetDegs, visualizedResponseScalingDegs)
    xx = mosaicEccDegs(1)+0.5*mosaicSizeDegs(1)*[-1 1 1 -1 -1];
    ymin = -visualizedResponseScalingDegs*0.1;
    ymax =  visualizedResponseScalingDegs*1.1;
    yy = horizontalSliceYcoordDegs + yOffsetDegs + [ymin ymin ymax ymax ymin];
    yGrid(1) = horizontalSliceYcoordDegs + yOffsetDegs + 0;
    yGrid(2) = yGrid(1) + visualizedResponseScalingDegs*0.25;
    yGrid(3) = yGrid(1) + visualizedResponseScalingDegs*0.5;
    yGrid(4) = yGrid(1) + visualizedResponseScalingDegs*0.75;
    yGrid(5) = yGrid(1) + visualizedResponseScalingDegs*1.0;
    vertices = [xx(:) yy(:)];
    patch(ax, 'Faces',1:size(vertices,1),'Vertices',vertices, 'EdgeColor', [1 1 0], 'EdgeAlpha', 0.7, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha', 0.85);
    %plot(ax, xx, xx*0+yGrid(1), 'w:','LineWidth', 0.5);
    %plot(ax, xx, xx*0+yGrid(2), 'w:','LineWidth', 0.5);
    plot(ax, xx, xx*0+yGrid(3), 'y:','LineWidth', 0.5);
    %plot(ax, xx, xx*0+yGrid(4), 'w:','LineWidth', 0.5);
    %plot(ax, xx, xx*0+yGrid(5), 'w:','LineWidth', 0.5);
end


function renderConeResponses(ax, theVisualizedConeXcoords, theVisualizedConeIndices, ...
    theConeMosaicResponse, horizontalSliceYcoord, visualizedResponseScalingDegs, maxResponse, coneColor)

    theVisualizedConeResponses = squeeze(theConeMosaicResponse(1,1,theVisualizedConeIndices));
    scatter(ax, theVisualizedConeXcoords, horizontalSliceYcoord + theVisualizedConeResponses/maxResponse*visualizedResponseScalingDegs, 100, ...
        'MarkerEdgeColor', coneColor, 'MarkerFaceColor', coneColor, 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1.0, 'LineWidth', 1.5);
end