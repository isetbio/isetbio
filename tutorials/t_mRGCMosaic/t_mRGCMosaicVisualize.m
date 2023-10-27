%% Introduction to the midget RGC mosaic (mRGCMosaic) object.
%
% Description:
%    Demonstrates
%        - how to visualize different aspects
%


% History:
%    10/25/23  NPC  Wrote it.

function t_mRGCMosaicVisualize
    %% Close all figures
    close all;

    %% Display available mRGCMosaics
    rgcMosaicType = 'ONcenterMidgetRGC';
    mRGCMosaic.availableComputeReadyMosaics(rgcMosaicType);

    %% Specify the desired eccentricity of the precomputed mRGC mosaic
    % Choose the x-eccentricity from one of the available mosaics,
    % displayed above
    % (e.g., -16.0 to load the mosaic 'mRGCMosaicEcDegs(-10.0_0.0)_SizeDegs(6.0_3.0)...'
    horizontalEccDegs = 7.0;

    %% Load the precomputed mRGCMosaic
    theMRGCMosaic = MosaicPoolingOptimizer.loadPreComputedMRGCMosaic(horizontalEccDegs);

    %% Visualize the mRGCMosaic
    % Lets visualize it together with its input cone mosaic.
    % In this visualization, the gray contours identify the mRGC RF centers, 
    % whereas cones connected to mRGCs are displayed by colored disks.
    theMRGCMosaic.visualize(...
        'identifyInputCones', ~true, ...
        'identifyPooledCones', ~true, ...
        'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
        'plotTitle', 'full mosaic');

    %% Crop the mRGCMosaic
    sizeDegs = [0.75 0.25];
    eccDegs = [6 0.5];
    theMRGCMosaic.cropToSizeAtEccentricity(sizeDegs, eccDegs);

    theMRGCMosaic.visualize(...
        'identifyInputCones', true, ...
        'identifyPooledCones', true, ...
        'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
        'centerSubregionContourSamples', 24, ...
        'plotTitle', 'cropped mosaic');

    %% Visualize the retinal cone pooling for a couple neurons
    % Lets visualize the retinal cone pooling of a couple of neurons.
    % For the first neuron visualization, find a unit with 2 L-cone
    % inputs in its RF center, and located near eccentricityDegs+[-0.5 0]; degs. 
    % 

    targetRGCposition = eccDegs + [-0.1 0.1];
    targetCenterConesNum = 2;
    targetCenterConeMajorityType = cMosaic.MCONE_ID;
    theRGCindex(1) = theMRGCMosaic.visualizeRetinalConePoolingRFmapNearPosition(...
       targetRGCposition, targetCenterConesNum, ...
       targetCenterConeMajorityType, ...
       'tickSeparationArcMin', 3);

    targetCenterConesNum = 3;
    targetCenterConeMajorityType = cMosaic.LCONE_ID;
    theRGCindex(2) = theMRGCMosaic.visualizeRetinalConePoolingRFmapNearPosition(...
       targetRGCposition, targetCenterConesNum, ...
       targetCenterConeMajorityType, ...
       'tickSeparationArcMin', 3);

    theMRGCMosaic.visualize(...
        'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
        'identifyInputCones', true, ...
        'identifyPooledCones', true, ...
        'labelRGCsWithIndices', theRGCindex, ...
        'labeledRGCsColor', [0 0 1], ...
        'labeledRGCsLineWidth', 4.0, ...
        'centerSubregionContourSamples', 24, ...
        'backgroundColor', [1 1 1]);

end

