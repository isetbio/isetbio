function renderSubregionConePoolingLineWeightingFunctions(ax, ...
            centerLineWeightingFunction, surroundLineWeightingFunction, ...
            sensitivityRange, horizontalAxisDirection, varargin)

    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('noXTicks', false, @islogical);
    p.addParameter('noYTicks', false, @islogical);
    p.addParameter('gridless', false, @islogical);
    p.addParameter('plotTitle', '', @ischar);
    p.addParameter('tickSeparationArcMin', [],  @(x)(isempty(x)||isscalar(x)));
    p.addParameter('spatialSupportRangeArcMin', [],  @(x)(isempty(x)||isscalar(x)));
    p.addParameter('xAxisTickAngleRotationDegs', 90, @isscalar)
    p.addParameter('withFigureFormat', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('resetAxes', true, @islogical);

    p.parse(varargin{:});
    
    spatialSupportRangeArcMin = p.Results.spatialSupportRangeArcMin;
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    plotTitle = p.Results.plotTitle;
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    noXTicks = p.Results.noXTicks;
    noYTicks = p.Results.noYTicks;
    gridless = p.Results.gridless;
    xAxisTickAngleRotationDegs = p.Results.xAxisTickAngleRotationDegs;
    ff = p.Results.withFigureFormat;
    resetAxes = p.Results.resetAxes;

    if (isempty(spatialSupportRangeArcMin))
        spatialSupportRangeArcMin = 10;
    end

    if (isempty(tickSeparationArcMin))
        tickSeparationArcMin = 3.0;
    end

    centerLineWeightingFunction.spatialSupportDegs = (centerLineWeightingFunction.spatialSupportDegs(:))';
    centerLineWeightingFunction.amplitude = (centerLineWeightingFunction.amplitude(:))';
    surroundLineWeightingFunction.spatialSupportDegs = (surroundLineWeightingFunction.spatialSupportDegs(:))';
    surroundLineWeightingFunction.amplitude = (surroundLineWeightingFunction.amplitude(:))';
    
    xo = mean(centerLineWeightingFunction.spatialSupportDegs);
    centerLineWeightingFunction.spatialSupportArcMin = (centerLineWeightingFunction.spatialSupportDegs  - xo)*60;
    surroundLineWeightingFunction.spatialSupportArcMin = (surroundLineWeightingFunction.spatialSupportDegs - xo)*60;

    if (resetAxes)
        cla(ax, 'reset');
    end

    % Surround
    shadedAreaBetweenTwoLines(ax, surroundLineWeightingFunction.spatialSupportArcMin, ...
        -surroundLineWeightingFunction.amplitude, surroundLineWeightingFunction.amplitude*0, ...
        [0.35 0.35 0.35], [0.1 0.1 0.1], 0.3, 1.0, '-');
    hold(ax, 'on');

    % Center
    shadedAreaBetweenTwoLines(ax, centerLineWeightingFunction.spatialSupportArcMin, ...
        centerLineWeightingFunction.amplitude, centerLineWeightingFunction.amplitude*0, ...
        [0.85 0.85 0.85], [0.1 0.1 0.1], 0.3, 1.0, '-');
   
    % Composite in red line
    plot(ax,surroundLineWeightingFunction.spatialSupportArcMin, ...
        centerLineWeightingFunction.amplitude - surroundLineWeightingFunction.amplitude, ...
        'r-', 'LineWidth', 1.5);
            
    %plot(ax, surroundLineWeightingFunction.spatialSupportArcMin, surroundLineWeightingFunction.spatialSupportArcMin*0, 'k--', 'LineWidth', 0.5);
            
    XLims = spatialSupportRangeArcMin/2*[-1 1];
    xTicks = 0:(tickSeparationArcMin):60;
    xTicks = [-fliplr(xTicks(2:end)) xTicks];
    if (tickSeparationArcMin >= 6)
        xTickLabels = sprintf('%2.0f\n', xTicks);
    else
        xTickLabels = sprintf('%2.1f\n', xTicks);
    end

    sensitivityTicks = -1:0.2:1;
    set(ax, 'XLim', [XLims(1)+(XLims(2)-XLims(1))*ff.axisOffsetFactor XLims(2)], ...
            'YLim', [sensitivityRange(1)+(sensitivityRange(2)-sensitivityRange(1))*ff.axisOffsetFactor sensitivityRange(2)], ...
            'XTick', xTicks, 'YTick', sensitivityTicks , ...
            'XTickLabel', xTickLabels, ...
            'YTickLabel', sprintf('%2.1f\n', sensitivityTicks )); 

    if (noXTicks)
        set(ax, 'XTickLabel', {});
    end

    if (noYTicks)
        set(ax, 'YTickLabel', {});
    end

    if (isempty(ff))
        if (~noXLabel)
            xlabel(ax, sprintf('space, %s (arc min)', horizontalAxisDirection));
        end
    
        if (~noYLabel)
            if (strcmp(horizontalAxisDirection, 'x'))
                ylabel(ax, 'integrated (y) cone weights' ,'FontAngle', ff.axisFontAngle);
            else
                ylabel(ax, 'integrated (x) cone weights' ,'FontAngle', ff.axisFontAngle);
            end
        end

        if (~isempty(plotTitle))
            title(ax, plotTitle);
        end
    else
        % Font size
        set(ax, 'FontSize', ff.fontSize);
    
        % axis color and width
        set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);

        if (~noXLabel)
            xlabel(ax, sprintf('space, %s (arc min)', horizontalAxisDirection), 'FontAngle', ff.axisFontAngle);
        end
    
        if (~noYLabel)
            if (strcmp(horizontalAxisDirection, 'x'))
                ylabel(ax, 'integrated (y) cone weights' ,'FontAngle', ff.axisFontAngle);
            else
                ylabel(ax, 'integrated (x) cone weights' ,'FontAngle', ff.axisFontAngle);
            end
        end

        if (~isempty(plotTitle))
            title(ax, plotTitle, ...
                'Color', ff.titleColor, ...
                'FontSize', ff.titleFontSize, ...
                'FontWeight', ff.titleFontWeight);
        end

    end

    
    if (~gridless)
        grid(ax, 'on');
        box(ax, 'off');
    else
        grid(ax, 'off');
        box(ax, 'off');
    end

    
    axis(ax, 'square');
    xtickangle(ax, xAxisTickAngleRotationDegs);

end

function shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end
