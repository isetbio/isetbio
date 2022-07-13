% Demo usage of computing and using custom cone mosaics
%
% Description:
%    Shows how to generate a completely custom @cMosaic object in which the
%    position, type, and aperture size of each and every cone is specified
%    by the user. This functionality enables the modeling of real mosaics,
%    for example mosaic patches characterized via AO experiments.
%

% History:
%    07/13/22  NPC  ISETBIO Team, Copyright 2022 Wrote it.

function t_cMosaicCustomConeData()

    % All specifications are in units of retinal microns
    thePositionUnits = 'microns';

    % Specify a custom retinal magnification, here 250 microns/deg
    theCustomRetinalMagnification = 250;

    % Specify the (x,y) positions of all cones (in microns)
    % Here we are laying cones along a spiral path with cone aperture increasing
    % with distance from the center.
    angles = 0:15:3000;
    minConeDiameter = 2;
    for iCone = 1:numel(angles)
        iAngle = angles(iCone);
        radius = minConeDiameter + iCone*0.7;
        thePositions(iCone,:) = radius * [ cosd(iAngle) sind(iAngle)];
        theConeApertureDiameters(iCone) = minConeDiameter + (iCone*minConeDiameter*0.01)^1.2;
        switch (mod(iCone,11))
            case {1,2,4,5,7,8,9}
                theTypes(iCone) = cMosaic.LCONE_ID;  % an L-cone
            case {0,3,6}
                theTypes(iCone) = cMosaic.MCONE_ID;  % an M-cone
            case 10
                theTypes(iCone) = cMosaic.SCONE_ID;  % an S-cone
        end
    end
    

    % Generate struct with custom cone data
    theCustomConeDataStruct = struct(...
        'positionUnits', thePositionUnits, ...
        'positions', thePositions, ...
        'types', theTypes,...
        'lightGatheringApertureDiameters', theConeApertureDiameters ...                      
        );

     % Generate a @cMosaic from the completely custom cone data and retinal
     % magnification factor
     theConeMosaic = cMosaic(...
         'coneData', theCustomConeDataStruct, ...
         'micronsPerDegree', theCustomRetinalMagnification);

     % Generate a test stimulus
     theOI = generateTestStimulusAndItsOpticalImage;

     % Compute the cone excitation to this optical image
     theActivation = theConeMosaic.compute(theOI);

     % Visualize the retinal image of the test stimulus, the generated
     % @cMosaic, and its activation by the test stimulus

     hFig = figure(1); clf;
     set(hFig, 'Position', [10 10 1500 500], 'Color', [1 1 1]);
     ax = subplot(1,3,1);
     visualizeOpticalImage(theOI, ...
         'axesHandle', ax, ...
         'crossHairsAtOrigin', true, ...
         'displayRadianceMaps', false);
     set(ax, 'XLim', [-0.5 0.5], 'XTick', [-0.5 0 0.5], ...
             'YLim', [-0.5 0.5], 'YTick', [-0.5 0 0.5]);
     title(ax, 'optical image');

     ax = subplot(1,3,2);
     % Visualize the cone mosaic
     theConeMosaic.visualize(...
         'figureHandle', hFig, ...
         'axesHandle', ax, ...
         'crossHairsOnMosaicCenter', true, ...
         'domainVisualizationLimits', [-0.5 0.5 -0.5 0.5], ...
         'domainVisualizationTicks', struct('x', -0.5:0.5:0.5, 'y', -0.5:0.5:0.5), ...
         'plotTitle', 'custom @cMosaic');

     ax = subplot(1,3,3);
     % Visualize the cone mosaic activation
     theConeMosaic.visualize(...
         'figureHandle', hFig, ...
         'axesHandle', ax, ...
         'activation', theActivation, ...
         'crossHairsOnMosaicCenter', true, ...
         'domainVisualizationLimits', [-0.5 0.5 -0.5 0.5], ...
         'domainVisualizationTicks', struct('x', -0.5:0.5:0.5, 'y', -0.5:0.5:0.5), ...
         'horizontalActivationColorBarInside', true, ...
         'plotTitle', 'mosaic activation');


end

function theOI = generateTestStimulusAndItsOpticalImage
    %% Generate a grating scene (0.5 deg FOV, 20 c/deg)
    pixelsNum = 512;
    fovDegs = 1;
    stimFreqCyclesPerDeg = 1.25;
    
    parms.freq = stimFreqCyclesPerDeg*fovDegs;
    parms.contrast = 1;
    parms.ph = -pi/2;
    parms.ang = 0;
    parms.row = pixelsNum;
    parms.col = pixelsNum;
    parms.GaborFlag = 0;
    scene = sceneCreate('harmonic', parms);
    scene = sceneSet(scene, 'fov', fovDegs);
    
    %% Compute the optical image
    theOI = oiCreate;
    theOI = oiCompute(scene, theOI);

end
