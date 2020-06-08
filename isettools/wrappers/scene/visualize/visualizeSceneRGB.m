function visualizeSceneRGB(spatialSupport, spatialSupportUnits, ...
    RGBsettings, meanLuminance, meanChromaticity, sceneName)
% Method to visualize the scene as RGB
spatialSupportX = squeeze(spatialSupport(1,:,1));
spatialSupportY = squeeze(spatialSupport(:,1,2));
[~,midRow] = min(abs(spatialSupportY)); 
figure(); clf;
image(spatialSupportX, spatialSupportY, RGBsettings);
hold on;
plot([0 0], [spatialSupportY(1) spatialSupportY(end)], 'r-');
plot([spatialSupportX(1) spatialSupportX(end)], [0 0],'r-');
axis 'xy'; axis 'image';
xtickformat('%0.2f'); ytickformat('%0.2f');
set(gca, 'XTick', max(spatialSupportX)*[-1 -0.5 0 0.5 1]);
set(gca, 'YTick', max(spatialSupportY)*[-1 -0.5 0 0.5 1]);
xlabel(sprintf('\\it space (%s)', spatialSupportUnits));
ylabel(sprintf('\\it space (%s)', spatialSupportUnits));

% Label plot
if (isempty(meanLuminance))
    title(sprintf('%s\nmean chroma: (%2.3f,%2.3f)', ...
        sceneName, meanChromaticity(1), meanChromaticity(2)));
else
    title(sprintf('%s\nmean luminance: %2.1f cd/m2;\nmean chroma: (%2.3f,%2.3f)', ...
        sceneName, meanLuminance, meanChromaticity(1), meanChromaticity(2)));
end
set(gca, 'FontSize', 16);
drawnow;
end