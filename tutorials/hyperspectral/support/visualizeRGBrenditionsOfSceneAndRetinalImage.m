function visualizeRGBrenditionsOfSceneAndRetinalImage(theScene, theOI, visualizedXRangeDegs, visualizedYRangeDegs, figNo)
% Visualize RGB renditions of the input scene and of its retinal image
%
% 7/25/18  npc  Wrote it
%

    hFig = figure(figNo);
    set(hFig, 'Position', [10 10 700 900]);
    
    % Plot the scene at the top
    sceneRGB = sceneGet(theScene, 'rgb');
    spatialSizeDegs = [sceneGet(theScene, 'wAngular') sceneGet(theScene, 'hAngular') ]; 
    [rows, cols, ~] = size(sceneRGB);
    xSupportDegs = supportInDegs(cols, spatialSizeDegs(1));
    ySupportDegs = supportInDegs(rows, spatialSizeDegs(2)); 
    
    subplot('Position', [0.01 0.56, 0.98, 0.42]);
    image(xSupportDegs, ySupportDegs, sceneRGB)
    axis 'image'
    set(gca, 'XLim', visualizedXRangeDegs, 'YLim', visualizedYRangeDegs);
    set(gca, 'FontSize', 16, 'YTick', []);
    title('scene');
    
    % Plot the retinal image at the bottom
    retinalImageRGB = oiGet(theOI, 'rgb');
    spatialSizeDegs = [oiGet(theOI, 'wAngular') oiGet(theOI, 'hAngular') ]; 
    [rows, cols, ~] = size(retinalImageRGB);
    xSupportDegs = supportInDegs(cols, spatialSizeDegs(1));
    ySupportDegs = supportInDegs(rows, spatialSizeDegs(2)); 
    
    subplot('Position', [0.01 0.06, 0.98, 0.42]);
    image(xSupportDegs, ySupportDegs, retinalImageRGB)
    axis 'image'
    set(gca, 'XLim', visualizedXRangeDegs, 'YLim', visualizedYRangeDegs);
    set(gca, 'FontSize', 16, 'YTick', []);
    xlabel('space (degs)', 'FontWeight', 'bold');
    title('retinal image');
end
