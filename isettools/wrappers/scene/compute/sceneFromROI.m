% Generate scene by cropping an ROI from a source scene at a desired eccentricity
%
% Description:
%    Generate a new scene by cropping an ROI from a source scene at a desired 
%    eccentricity. The new scene is zero centered.
% 
% See Also:
%   wrappers/coneMosaic/coneMosaicHexRegForDesiredEcc
%   ISETBIO LiveScripts repository: ls_computingWithEccentricityVaryingHexMosaics.mlx

% History:
%    7/9/19  NPC  ISETBIO Team, Copyright 2019

function roiScene = sceneFromROI(scene, roiRect)

    % retrieve the spatial support of the scene(in millimeters)
    spatialSupportMilliMeters = sceneGet(scene, 'spatial support', 'mm');

    viewingDistance = sceneGet(scene, 'distance');
    spatialSupportDegs = 2 * atand(spatialSupportMilliMeters/1e3/2/viewingDistance);
    spatialSupportDegsX = squeeze(spatialSupportDegs(1,:,1));
    spatialSupportDegsY = squeeze(spatialSupportDegs(:,1,2));
    
    xLimits = roiRect.xo + roiRect.width/2*[-1 1];
    yLimits = roiRect.yo + roiRect.height/2*[-1 1];
    cols = find((spatialSupportDegsX >= xLimits(1)) & (spatialSupportDegsX <= xLimits(2)));
    rows = find((spatialSupportDegsY >= yLimits(1)) & (spatialSupportDegsY <= yLimits(2)));
    
    scenePhotonRate = sceneGet(scene, 'photons');
    roiScenePhotonRate = scenePhotonRate(rows,cols,:);
    
    
    % Generate new scene
    roiScene = scene;
    
    % zero photons everywhere
    scenePhotonRate = scenePhotonRate*0;
    
    % blit roiScenePhotonRate to cental part
    y1 = round(size(scenePhotonRate,1)/2 - numel(rows)/2);
    x1 = round(size(scenePhotonRate,2)/2 - numel(cols)/2);
    
    yy = y1 + (1:numel(rows));
    xx = x1 + (1:numel(cols));
    scenePhotonRate(yy,xx,:) = roiScenePhotonRate;
    
    % Replace photons
    roiScene = sceneSet(roiScene, 'photons', scenePhotonRate);
end
