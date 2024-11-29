function t_cMosaicSuperimposeArbitraryRGBimageOnTopOfConeMosaic
%% Illustrate how to superimpose an arbitrary RGB image on a @cMosaic

%   10/09/24  NPC      Wrote it.

    % Stimulus scene
    imageFOVdegs = 1.0;
    pixelsPerCheck = 256;
    numberOfChecks = 4;

    theStimulusScene = sceneCreate('checkerboard', pixelsPerCheck, numberOfChecks);
    theStimulusScene = sceneSet(theStimulusScene, 'fov', imageFOVdegs);

    % Generate a @cMosaic
    c = cMosaic();

    % Generate default OI
    theOI = oiCreate;

    % Compute the retinal image of the stimulus
    theStimulusRetinalImage = oiCompute(theOI,theStimulusScene);

    % Visualize the retinal image on top of the mosaic
    c.visualize(...
        'withSuperimposedOpticalImage', theStimulusRetinalImage, ...
        'withSuperimposedOpticalImageAlpha', 0.3);

    % Alter the RGB content of the retinal image
    theOpticalImageRGB = oiGet(theStimulusRetinalImage, 'rgbimage');

    % Make it gray scale
    theGrayScaleOpticalImage = sum(theOpticalImageRGB,3);
    theGrayScaleRGBOpticalImage = [];
    theGrayScaleRGBOpticalImage(:,:,1) = theGrayScaleOpticalImage;
    theGrayScaleRGBOpticalImage(:,:,2) = theGrayScaleOpticalImage;
    theGrayScaleRGBOpticalImage(:,:,3) = theGrayScaleOpticalImage;

    % Use the appropriate spatial support
    theSpatialSupportMeters = oiGet(theStimulusRetinalImage, 'spatial support');
    spatialSupportMicrons = squeeze(theSpatialSupportMeters(1,1:end,1)) * 1e6;

    % Alternatively change the spatial support
    spatialSupportMicrons = linspace(-200,200,size(theGrayScaleRGBOpticalImage,1));

    % Visualize the altered retinal image on top of the mosaic
    c.visualize(...
        'withSuperimposedRGBopticalImage', theGrayScaleRGBOpticalImage, ...
        'withSuperimposedRGBopticalImageSpatialSupportMicrons', spatialSupportMicrons, ...
        'withSuperimposedRGBopticalImageAlpha', 0.5);

end