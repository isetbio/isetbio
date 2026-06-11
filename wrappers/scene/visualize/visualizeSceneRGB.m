function ax = visualizeSceneRGB(spatialSupport, spatialSupportUnits, ...
    RGBsettings, meanLuminance, meanChromaticity, sceneName, varargin)
p = inputParser;
p.addParameter('axesHandle', []);
p.addParameter('noTitle', false, @islogical);
p.addParameter('noYLabel', false, @islogical);
p.addParameter('noYTicks', false, @islogical);
p.addParameter('avoidAutomaticRGBscaling', false, @islogical);
% Parse input
p.parse(varargin{:});

% Method to visualize the scene as RGB
spatialSupportX = squeeze(spatialSupport(1,:,1));
spatialSupportY = squeeze(spatialSupport(:,1,2));
[~,midRow] = min(abs(spatialSupportY)); 


if (isempty(p.Results.axesHandle))
    figure(); clf;
    ax = subplot(1,1,1);
else
    ax = p.Results.axesHandle;
end

if (p.Results.avoidAutomaticRGBscaling)
    % Add a white and a black pixel to aid in rendering 
    %RGBsettings = 0.5*(RGBsettings / max(RGBsettings(:)));
    RGBsettings(1,1,1:3) = 0;
    RGBsettings(end,end,1:3) = 1.0;
end

image(ax,spatialSupportX, spatialSupportY, lrgb2srgb(RGBsettings));

axis(ax, 'image');
xtickformat('%0.2f'); ytickformat('%0.2f');
set(ax, 'XTick', max(spatialSupportX)*[-1 -0.5 0 0.5 1]);
if (p.Results.noYTicks)
    set(ax, 'YTick', []);
else
    set(ax, 'YTick', max(spatialSupportY)*[-1 -0.5 0 0.5 1]);
end

xlabel(ax,sprintf('space (%s)', spatialSupportUnits));
if (~p.Results.noYLabel)
    ylabel(ax,sprintf('space (%s)', spatialSupportUnits));
end

% Label plot
if (~p.Results.noTitle)
    if (isempty(meanLuminance))
        title(ax,sprintf('%s\nmean chroma: (%2.3f,%2.3f)', ...
            sceneName, meanChromaticity(1), meanChromaticity(2)));
    else
        title(ax,sprintf('%s\nmean luminance: %2.1f cd/m2;\nmean chroma: (%2.3f,%2.3f)', ...
            sceneName, meanLuminance, meanChromaticity(1), meanChromaticity(2)));
    end
end

set(ax, 'FontSize', 16);
drawnow;
end