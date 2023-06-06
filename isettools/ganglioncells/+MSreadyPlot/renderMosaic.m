function renderMosaic(hFig, theMosaicOBJ, plotTitle, varargin)

    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('plotRFoutlines', true, @islogical);
    p.addParameter('identifyPooledCones', true, @islogical);
    p.addParameter('identifyInputCones', true, @islogical);
    p.addParameter('backgroundColor', [], @(x)((ischar(x)&&(strcmp(x,'none')))||isempty(x)||((isvector(x))&&(numel(x) == 3))));
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    identifyPooledCones = p.Results.identifyPooledCones;
    identifyInputCones = p.Results.identifyInputCones;
    plotRFoutlines = p.Results.plotRFoutlines;
    backgroundColor = p.Results.backgroundColor;

    ff = MSreadyPlot.figureFormat('1x1 long');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    ax = theAxes{1,1};

    xMin =  theMosaicOBJ.eccentricityDegs(1) - 0.5*theMosaicOBJ.sizeDegs(1);
    xMax =  theMosaicOBJ.eccentricityDegs(1) + 0.5*theMosaicOBJ.sizeDegs(1);
    yMin =  theMosaicOBJ.eccentricityDegs(2) - 0.5*theMosaicOBJ.sizeDegs(2);
    yMax =  theMosaicOBJ.eccentricityDegs(2) + 0.5*theMosaicOBJ.sizeDegs(2);
    domainVisualizationLimits = [xMin xMax yMin yMax];
    domainVisualizationTicks = struct(...
        'x', [xMin theMosaicOBJ.eccentricityDegs(1) xMax], ...
        'y', [yMin theMosaicOBJ.eccentricityDegs(2) yMax]);


    theMosaicOBJ.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'identifyPooledCones', identifyPooledCones, ...
        'identifyInputCones', identifyInputCones, ...
        'plotRFoutlines', plotRFoutlines, ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'backgroundColor', backgroundColor);

    % xy axis labels
    if (~noXLabel)
        xlabel(ax, 'space (degs)', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'XTickLabel', {});
    end
    if (~noYLabel)
        ylabel(ax, 'space (degs)', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'YTickLabel', {});
    end

    % ticks and grids
    grid(ax, 'on')
    xtickangle(ax, 0);
    set(ax, 'TickDir', 'both');

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    % plot title
    if (~isempty(plotTitle))
        title(ax, plotTitle, 'Color', ff.titleColor, 'FontSize', ff.titleFontSize, 'FontWeight', ff.titleFontWeight);
    end

end
