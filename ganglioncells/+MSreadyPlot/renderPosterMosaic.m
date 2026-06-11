function renderPosterMosaic(hFig, theMosaicOBJ, plotTitle, ff, varargin)

    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('plotRFoutlines', true, @islogical);
    p.addParameter('identifyPooledCones', true, @islogical);
    p.addParameter('identifyInputCones', true, @islogical);
    p.addParameter('domainVisualizationLimits', [], @isnumeric);
    p.addParameter('domainVisualizationTicks', [],  @(x)(isempty(x)||isstruct(x)));
    p.addParameter('backgroundColor', [], @(x)((ischar(x)&&(strcmp(x,'none')))||isempty(x)||((isvector(x))&&(numel(x) == 3))));
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    identifyPooledCones = p.Results.identifyPooledCones;
    identifyInputCones = p.Results.identifyInputCones;
    plotRFoutlines = p.Results.plotRFoutlines;
    backgroundColor = p.Results.backgroundColor;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;

    
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    ax = theAxes{1,1};

    if (isempty(domainVisualizationLimits))
        xMin =  theMosaicOBJ.eccentricityDegs(1) - 0.5*theMosaicOBJ.sizeDegs(1);
        xMax =  theMosaicOBJ.eccentricityDegs(1) + 0.5*theMosaicOBJ.sizeDegs(1);
        yMin =  theMosaicOBJ.eccentricityDegs(2) - 0.25;
        yMax =  theMosaicOBJ.eccentricityDegs(2) + 0.25;
        domainVisualizationLimits = [xMin xMax yMin yMax];
    end

    if isempty(domainVisualizationTicks)
        xMin =  theMosaicOBJ.eccentricityDegs(1) - 0.5*theMosaicOBJ.sizeDegs(1);
        xMax =  theMosaicOBJ.eccentricityDegs(1) + 0.5*theMosaicOBJ.sizeDegs(1);
        domainVisualizationTicks = struct(...
            'x', [xMin theMosaicOBJ.eccentricityDegs(1) xMax], ...
            'y', [yMin theMosaicOBJ.eccentricityDegs(2) yMax]);
    end

    theMosaicOBJ.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'identifiedConeApertureThetaSamples', 32, ...
        'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
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
    grid(ax, 'off');

    box(ax, 'on');

    xtickangle(ax, 0);
    set(ax, 'TickDir', 'out', 'TickLength', ff.tickLength);

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    % plot title
    if (~isempty(plotTitle))
        title(ax, plotTitle, 'Color', ff.titleColor, 'FontSize', ff.titleFontSize, 'FontWeight', ff.titleFontWeight);
    end

end
