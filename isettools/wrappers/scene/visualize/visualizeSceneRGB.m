function ax = visualizeSceneRGB(spatialSupport, spatialSupportUnits, ...
    RGBsettings, meanLuminance, meanChromaticity, sceneName, varargin)
p = inputParser;
p.addParameter('axesHandle', []);
p.addParameter('noTitle', false, @islogical);
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

image(ax,spatialSupportX, spatialSupportY, RGBsettings);
axis(ax, 'image');
xtickformat('%0.2f'); ytickformat('%0.2f');
set(ax, 'XTick', max(spatialSupportX)*[-1 -0.5 0 0.5 1]);
set(ax, 'YTick', max(spatialSupportY)*[-1 -0.5 0 0.5 1]);
xlabel(ax,sprintf('space (%s)', spatialSupportUnits));
ylabel(ax,sprintf('space (%s)', spatialSupportUnits));

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