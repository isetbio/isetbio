function sceneOrOpticalImage(theSceneOrOpticalImage, ...
    thePresentationDisplay, ax, domainVisualizationLimits, ...
    domainVisualizationTicks, varargin)

    p = inputParser;
    p.addParameter('fontSize', 16, @isscalar);
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('plotTitle', '', @ischar);

    p.parse(varargin{:});
    fontSize = p.Results.fontSize;
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    plotTitle = p.Results.plotTitle;

    % Compute the RGB settings for the display
    displayLinearRGBToXYZ = displayGet(thePresentationDisplay, 'rgb2xyz');
    displayXYZToLinearRGB = inv(displayLinearRGBToXYZ);
    
    % Extract the XYZ image representation
    if (strcmp(theSceneOrOpticalImage.type, 'scene'))
        % its a scene
        xyzImage = sceneGet(theSceneOrOpticalImage, 'xyz');
        pixelSizeDegs = sceneGet(theSceneOrOpticalImage,'w angular resolution');
    else
        % its an optical image.
        % Increase the photons for better visibility
        %retinalIlluminanceBoostFactor = 1;
        %theSceneOrOpticalImage.data.photons = theSceneOrOpticalImage.data.photons * retinalIlluminanceBoostFactor;
        xyzImage = oiGet(theSceneOrOpticalImage, 'xyz');
        pixelSizeDegs = oiGet(theSceneOrOpticalImage,'w angular resolution');
    end

    displaySettingsImage = xyz2srgb(xyzImage);

    xPixels = size(xyzImage,2);
    yPixels = size(xyzImage,1);
    x = 1:xPixels;
    y = 1:yPixels;
    x = x-mean(x);
    y = y-mean(y);
    x = x*pixelSizeDegs;
    y = y*pixelSizeDegs;

   
    image(ax,x,y,displaySettingsImage);
    axis(ax, 'image');
    set(ax, 'CLim', [0 1]);

    if (~isempty(domainVisualizationLimits))
        set(ax, 'XLim',domainVisualizationLimits(1:2), 'YLim', domainVisualizationLimits(3:4));
    end

    if (~isempty(domainVisualizationTicks))
        set(ax, 'XTick', domainVisualizationTicks.x, ...
                'XTickLabel', sprintf('%.2f\n', domainVisualizationTicks.x), ...
                'YTick', domainVisualizationTicks.y, ...
                'YTickLabel', sprintf('%.2f\n', domainVisualizationTicks.y))
    end

    set(ax, 'FontSize', fontSize);
    if (~noXLabel)   
        xlabel(ax, 'space, x (degs)');
    end
    if (~noYLabel)   
        ylabel(ax, 'space, y (degs)');
    end
    if (~isempty(plotTitle))
        title(ax, plotTitle);
    end
end
