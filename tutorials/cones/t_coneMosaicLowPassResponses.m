% Demonstrates how to low-pass a mosaic's response
%
% 
% NPC, ISETBIO Team, 2016  
%

function t_coneMosaicLowPassResponses()

    % Set up scene
    scene = sceneCreate; 
    scene = sceneSet(scene, 'h fov', 1.0);
    oi = oiCreate('human'); 
    oi = oiCompute(oi, scene);
   
    % Set up mosaic
    cMosaicOBJ = coneMosaic();
    cMosaicOBJ.setSizeToFOV([sceneGet(scene, 'h fov'), sceneGet(scene, 'v fov')]);
    cMosaicOBJ.noiseFlag = 'none';
    [isomerizationsMap, currentsMap] = cMosaicOBJ.compute(oi,'currentFlag', true);
    
    % Low-pass the isomerizations map
    lowPassSpaceConstantInMicrons = [5 5 10];
    [isomerizationsMapLowPassed, Lmap, Mmap, Smap] = cMosaicOBJ.lowPassMosaicResponse(isomerizationsMap, lowPassSpaceConstantInMicrons);
    
    % Visualize results
    climRange = [min([min(isomerizationsMap(:)) min(isomerizationsMapLowPassed(:))]) max([max(isomerizationsMap(:)) max(isomerizationsMapLowPassed(:))])];
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 500 1024 500]);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 3, ...
           'heightMargin',   0.08, ...
           'widthMargin',    0.05, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.01, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.03);
       
    subplot('Position', subplotPosVectors(1,1).v)
    imshow(sceneGet(scene, 'RGB'))
    axis 'image'
    title('input scene');
    
    subplot('Position', subplotPosVectors(1,2).v)
    imagesc(isomerizationsMap);
    set(gca, 'CLim', climRange);
    axis 'image'
    title(sprintf('original response (max = %2.1f)', max(isomerizationsMap(:))));
    
    subplot('Position', subplotPosVectors(1,3).v)
    imagesc(isomerizationsMapLowPassed);
    set(gca, 'CLim', climRange);
    axis 'image'
    title(sprintf('lowpassed response (max = %2.1f)', max(isomerizationsMapLowPassed(:))));
    
    subplot('Position', subplotPosVectors(2,1).v)
    imagesc(Lmap);
    axis 'image'
    title(sprintf('L demosaiced response (lowpassed, tau = %2.1f um)', lowPassSpaceConstantInMicrons(1)));
    
    subplot('Position', subplotPosVectors(2,2).v)
    imagesc(Mmap);
    axis 'image'
    title(sprintf('M demosaiced response (lowpassed, tau = %2.1f um)', lowPassSpaceConstantInMicrons(2)));
   
    subplot('Position', subplotPosVectors(2,3).v)
    imagesc(Smap);
    axis 'image'
    title(sprintf('S demosaiced response (lowpassed, tau = %2.1f um)', lowPassSpaceConstantInMicrons(3)));
    
    colormap(gray(1024));
end