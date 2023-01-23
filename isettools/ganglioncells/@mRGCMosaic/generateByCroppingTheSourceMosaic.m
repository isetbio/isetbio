% Method to generate the mRGCMosaic by cropping the sourceMidgetRGCMosaic
function generateByCroppingTheSourceMosaic(obj, sourceMidgetRGCMosaic, varargin)
    p = inputParser;
    p.addParameter('visualizeSpatialRelationshipToSourceMosaic', [], @(x)(isempty(x))||(isscalar(x)));

    % Parse input
    p.parse(varargin{:});
    visualizeSpatialRelationshipToSourceMosaic = p.Results.visualizeSpatialRelationshipToSourceMosaic;

    % Get the input cone mosaic
    obj.inputConeMosaic = sourceMidgetRGCMosaic.inputConeMosaic;
 
    if ((isempty(obj.eccentricityDegs))&&(isempty(obj.sizeDegs))) || ...
       ((all(obj.eccentricityDegs==sourceMidgetRGCMosaic.eccentricityDegs))&&(all(obj.sizeDegs == sourceMidgetRGCMosaic.sizeDegs)))
        % When empty just set to the eccentricity and size of the 
        % sourceMidgetRGCMosaic
        obj.eccentricityDegs = sourceMidgetRGCMosaic.eccentricityDegs;
        obj.sizeDegs = sourceMidgetRGCMosaic.sizeDegs;

        % Get the rf positions and cone pooling matrices of the entire
        % sourceMidgetRGCMosaic
        obj.rgcRFpositionsDegs = sourceMidgetRGCMosaic.rgcRFpositionsDegs;
        obj.centerConePoolingMatrix = sourceMidgetRGCMosaic.rgcRFcenterConePoolingMatrix;
        obj.surroundConePoolingMatrix = sourceMidgetRGCMosaic.rgcRFsurroundConePoolingMatrix;
    else
        % Instantiate a rectangular ROI with the desired size and position
        theROI = regionOfInterest(...
            'geometryStruct', struct(...
                'units', 'degs', ...
                'shape', 'rect', ...
                'center', obj.eccentricityDegs, ...
                'width', obj.sizeDegs(1), ...
                'height', obj.sizeDegs(2), ...
                'rotation', 0 ...
                ));

        idx = theROI.indicesOfPointsInside(sourceMidgetRGCMosaic.rgcRFpositionsDegs);

        % Crop now
        obj.rgcRFpositionsDegs = sourceMidgetRGCMosaic.rgcRFpositionsDegs(idx,:);
        obj.centerConePoolingMatrix = sourceMidgetRGCMosaic.rgcRFcenterConePoolingMatrix(:,idx);
        obj.surroundConePoolingMatrix = sourceMidgetRGCMosaic.rgcRFsurroundConePoolingMatrix(:,idx);
    end

    obj.inputConesNum = size(obj.centerConePoolingMatrix,1);
    obj.rgcsNum = size(obj.centerConePoolingMatrix,2);

    if (visualizeSpatialRelationshipToSourceMosaic)
        visualizeMosaicBorders(obj,sourceMidgetRGCMosaic, theROI)
    end

end

function visualizeMosaicBorders(obj, sourceMidgetRGCMosaic, theROI)
    
    sourceConesXYmin = min(obj.inputConeMosaic.coneRFpositionsDegs, [], 1);
    sourceConesXYmax = max(obj.inputConeMosaic.coneRFpositionsDegs, [], 1);
    
    theRGCXYmin = min(obj.rgcRFpositionsDegs, [],1);
    theRGCXYmax = max(obj.rgcRFpositionsDegs, [],1);

    sourceRGCXYmin = min(sourceMidgetRGCMosaic.rgcRFpositionsDegs, [],1);
    sourceRGCXYmax = max(sourceMidgetRGCMosaic.rgcRFpositionsDegs, [],1);


    hFig = figure(); clf;
    set(hFig, 'Position', [10 10 1500 800], 'Color', [1 1 1]);
    ax = subplot('Position', [0.05 0.08 0.5 0.90]);
    theROI.visualize(...
       'figureHandle', hFig, ...
       'axesHandle', ax, ...
       'fillColor', [0.8 0.8 0.2]);
    hold(ax, 'on');
    pCones = plotRect(ax, sourceConesXYmin(1), sourceConesXYmax(1), sourceConesXYmin(2), sourceConesXYmax(2), [0.5 0.5 0.5]);
    pSourceRGC = plotRect(ax, sourceRGCXYmin(1), sourceRGCXYmax(1), sourceRGCXYmin(2), sourceRGCXYmax(2), [0 0 1]);
    pRGC = plotRect(ax, theRGCXYmin(1), theRGCXYmax(1), theRGCXYmin(2), theRGCXYmax(2), [1 0 0]);
    legend(ax, [pCones pSourceRGC pRGC], {'input cone mosaic', 'source RGC mosaic', 'RGC mosaic'}, ...
        'Location', 'NorthOutside', 'Orientation', 'horizontal');
    axis(ax, 'equal');
    set(ax, 'XLim', [sourceConesXYmin(1)-0.2 sourceConesXYmax(1)+0.2], 'YLim', [sourceConesXYmin(2)-0.2 sourceConesXYmax(2)+0.2]);
    set(ax, 'XTick', -20:0.5:20, 'YTick', -20:0.5:20);
    
    grid(ax, 'on');
    xtickangle(ax, 0);
    set(ax, 'FontSize', 16, 'TickDir', 'both');
    xlabel(ax, 'x (degs)');
    ylabel(ax, 'y (degs)');

    
    ax = subplot('Position', [0.58 0.08 0.4 0.90]);
    pRGC = plotRect(ax, theRGCXYmin(1), theRGCXYmax(1), theRGCXYmin(2), theRGCXYmax(2), [1 0 0]);
    hold(ax, 'on');
    plot(obj.rgcRFpositionsDegs(:,1), obj.rgcRFpositionsDegs(:,2), 'ko');
    axis(ax, 'equal');
    set(ax, 'XLim', [theRGCXYmin(1)-0.2 theRGCXYmax(1)+0.2], 'YLim', [theRGCXYmin(2)-0.2 theRGCXYmax(2)+0.2]);
    set(ax, 'XTick', -20:0.2:20, 'YTick', -20:0.2:20);
    grid(ax, 'on');
    xtickangle(ax, 0);
    set(ax, 'FontSize', 16, 'TickDir', 'both');
    xlabel(ax, 'x (degs)');
    ylabel(ax, 'y (degs)');
    title('cropped mRGCmosaic')
end

function p = plotRect(ax, xMin, xMax, yMin, yMax, color)
    x = [xMin xMin xMax xMax xMin];
    y = [yMin yMax yMax yMin yMin];
    plot(ax, x, y, '-', 'Color', [0 0 0], 'LineWidth', 2);
    hold(ax, 'on')
    p = plot(ax, x, y, '--', 'Color', color, 'LineWidth', 2.0);
end


