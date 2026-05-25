% Demo using cMosaic.oiEnsembleGenerate() to compute under optics at
% several positions within a single cone mosaic
%
% Description:
% Demonstrate usage of cMosaic.oiEnsembleGenerate(): how to compute cone
% mosaic activation for a large mosaic, with optics sampled at a grid
% of positions within the mosaic and the mosaic's activation computed by merging
% weighted sums of activations to multiple retinal images of a single scene
% via the sampled optics
%
% See Also:
%   t_cMosaicOffAxisDistortion
%   t_cMosaicEccVaryingOptics

% History:
%    05/18/26  NPC  ISETBIO Team, Copyright 2026 Wrote it.

%% Initialize
ieInit;
clear;
close all;


% Load an fMRI stimulus image
theImageIndex = 1;
[spatialSupportXDegs, spatialSupportYDegs, theSelectedImage] = loadSampleImages(theImageIndex);


% Generate a cone mosaic that extends over the entire stimulus
% area (and beyond)
mosaicEccDegs = [5.5 0];
theConeMosaic = cMosaic(...
        'sizeDegs', [4 18], ...
        'eccentricityDegs',  mosaicEccDegs);



% Since we are dealing with a large stimulus, we will need optics 
% sampled at a number of positions
% Specify an oiSamplingGridStruct for the optics 
% This will generate a regular hexagonal grid with specified 
% spacing, height, width at the specified eccentricity 

oiSamplingGridData = struct(...
    'eccentricityDegs', mosaicEccDegs, ...
    'widthDegs', 3, ...
    'heightDegs', 16, ...
    'spacingDegs', 2, ...                % sampling grid spacing
    'weightingType', 'raised cosine' ... % choose between {'Gaussian', 'raised cosine'};
    );

% Or if we just wanted to use the optics at one position only (center
% of cone mosaic)
% oiSamplingGridData = theConeMosaic.eccentricityDegs;

% Generate an ISETBio scene for the selected gray-scale image
[theScene, theBackgroundScene] = ISETBioSceneFromGrayScaleImageOnDisplay(theConeMosaic.wave, ...
        spatialSupportXDegs, spatialSupportYDegs, theSelectedImage);


% Pick a Polans subject ranking (1..10)
PolansSubjectRankOrder = 4;

% Optics subject
rankedSujectIDs = PolansOptics.constants.subjectRanking;
testSubjectID = rankedSujectIDs(PolansSubjectRankOrder);


% Whether to subtract the central refraction if needed
subtractCentralRefractionAsNeeded = true;


if (subtractCentralRefractionAsNeeded)
    % correct for central refraction
    subtractCentralRefraction = ...
        PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);
else
    % Or, alternatively do not correct for central refraction
    subtractCentralRefraction = false;
end

% Generate the optics at the desired sampling grid
[oiEnsemble, psfEnsemble, ~, oiSamplingGridDegs, theMergingWeights] = theConeMosaic.oiEnsembleGenerate(...
    oiSamplingGridData, ...
    'zernikeDataBase', 'Polans2015', ...
    'subjectID', testSubjectID, ...
    'subtractCentralRefraction', subtractCentralRefraction, ...
    'zeroCenterPSF', false, ...
    'withZeroedPistonAndTiltZernikeCoefficients', false, ...
    'wavefrontSpatialSamples', 601, ...
    'pupilDiameterMM', 3.0, ...
    'refractiveErrorDiopters', 0.0, ...
    'visualizedSamplingGrid', ~true);

% Visualize the cone mosaic, marking with white 'x' the positions where the
% optics are sampled
hFig = figure(10); clf;
set(hFig, 'Position', [10 10 700 1300], 'Color', [1 1 1]);
ax = subplot('Position', [0.05 0.05 0.94 0.94]);
theConeMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'plotTitle', ' ');
hold(ax, 'on');
plot(ax, oiSamplingGridDegs(:,1), oiSamplingGridDegs(:,2), 'wx', 'MarkerSize', 12, 'LineWidth', 1.5);


% Visualize the PSFs of the sampled optics
visualizeThePSFs(theConeMosaic, psfEnsemble, oiSamplingGridDegs);

if (isstruct(oiSamplingGridData))
    % Multiple position optics

    % Visualize the merging weights with which weighted averages of cone
    % mosaic activations from the same scene but different optics are computed
    for oiPosIndex = 1:size(theMergingWeights,1)
        visualizeTheMergingWeights(theConeMosaic, theMergingWeights(oiPosIndex,:), oiSamplingGridDegs);
    end

    % Compute cone mosaic activations to the retinal images of the scene computed for
    % each OI in the oiEnsemble, and merge the comuted activations using theMergingWeights
    multiOImergedConeMosaicActivation = theConeMosaic.computeForOIensemble(...
        oiEnsemble, theMergingWeights, theScene);

    % Compute cone mosaic activations to the retinal images of the background scene computed for
    % each OI in the oiEnsemble, and merge the comuted activations using theMergingWeights
    multiOImergedConeMosaicBackgroundActivation = theConeMosaic.computeForOIensemble(...
        oiEnsemble, theMergingWeights, theBackgroundScene);

    % Compute cone modulations
    coneMosaicModulations = (multiOImergedConeMosaicActivation - multiOImergedConeMosaicBackgroundActivation)./multiOImergedConeMosaicBackgroundActivation;
else
    % Single position optics
    theOI = oiEnsemble{1};

    % Compute the retinal image of the scene under this OI
    theRetinalImage = oiCompute(theOI, theScene, 'pad value','mean');

    % Compute the cone mosaic activation for the current OI
    theConeMosaicActivation = theConeMosaic.compute(...
        theRetinalImage, ...
        'opticalImagePositionDegs', [0 0]);

    % Compute the retinal image of the background scene under this OI
    theRetinalImage = oiCompute(theOI, theBackgroundScene, 'pad value','mean');

    % Compute the cone mosaic activation for the current OI
    theConeMosaicBackgroundActivation = theConeMosaic.compute(...
        theRetinalImage, ...
        'opticalImagePositionDegs', [0 0]);

    % Compute cone modulations
    coneMosaicModulations = (theConeMosaicActivation - theConeMosaicBackgroundActivation)./theConeMosaicBackgroundActivation;
end

% Visualized the loaded weightedConeMosaicModulations
visualizeConeMosaicActivation(theConeMosaic, coneMosaicModulations);


%
% Support routines for visualization and data loading
%

function hFig = visualizeTheMergingWeights(theConeMosaic, theMergingWeights, oiSamplingGridDegs)
    hFig = figure(10); clf;
    set(hFig, 'Position', [10 10 700 1300], 'Color', [1 1 1]);
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);
    theConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'activation', reshape(theMergingWeights, [1 1 numel(theMergingWeights)]), ...
        'activationRange', [0 1], ...
        'plotTitle', ' ');
    hold(ax, 'on');
    plot(ax, oiSamplingGridDegs(:,1), oiSamplingGridDegs(:,2), 'wx', 'MarkerSize', 12, 'LineWidth', 1.5);
end

function  hFig = visualizeThePSFs(theConeMosaic, psfEnsemble, oiSamplingGridDegs)

    opticsSamplingPositionsNum = size(oiSamplingGridDegs, 1);
    hFig = figure(3002); clf;
    if (theConeMosaic.sizeDegs(1)>theConeMosaic.sizeDegs(2))
        % wide mosaic
        figureWidthPixels = 1200;
        aspectRatio = theConeMosaic.sizeDegs(1)/theConeMosaic.sizeDegs(2);
        figureHeightPixels = round(figureWidthPixels/aspectRatio);
    else
        figureHeightPixels = 1200;
        aspectRatio = theConeMosaic.sizeDegs(2)/theConeMosaic.sizeDegs(1);
        figureWidthPixels = round(figureHeightPixels/aspectRatio);
    end

    set(hFig, 'Position', [10 10 figureWidthPixels figureHeightPixels], 'Color', [1 1 1]);

    for oiPos = 1:opticsSamplingPositionsNum
        width = 0.2;
        axPosition(1) = 0.5*(1+(oiSamplingGridDegs(oiPos,1) - theConeMosaic.eccentricityDegs(1))/(0.52*theConeMosaic.sizeDegs(1)))-0.5*width;
        axPosition(2) = 0.5*(1+(oiSamplingGridDegs(oiPos,2) - theConeMosaic.eccentricityDegs(2))/(0.52*theConeMosaic.sizeDegs(2))) - 0.5*width;
        axPosition(3:4) = width;

        ax = axes('Position', axPosition);
        set(ax, 'Color', [0 0 0]);
        thePSF = psfEnsemble{oiPos};
        [~, wIdx] = min(abs(thePSF.supportWavelength-550));
        wavePSF = squeeze(thePSF.data(:,:,wIdx));
        zLevels = 0.1:0.1:0.9;
        % half a degree
        xyRangeArcMin = 10*[-1 1];
        PolansOptics.renderPSF(ax, ...
            thePSF.supportX, thePSF.supportY, wavePSF/max(wavePSF(:)), ...
            xyRangeArcMin, zLevels,  gray(1024), [0 0 0], ...
            'plotTitle',  sprintf('%2.1f,%2.1f', oiSamplingGridDegs(oiPos,1), oiSamplingGridDegs(oiPos,2)));
        box(ax, 'on')
        set(ax, 'XColor', 0.2*[1 1 1], 'YColor', 0.2*[1 1 1]);
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        xlabel(ax, '');
        ylabel(ax, '');
        colormap(ax, gray(1024));
        drawnow;
    end
    pause(0.5)
end

function visualizeConeMosaicActivation(theConeMosaic, theConeMosaicActivation)

    hFig = figure(3000); clf;
    set(hFig, 'Position', [10 10 430 1300]);
    set(hFig, 'Color', [1 1 1]);
    ax = subplot('Position', [0.06 0.07 0.93 0.93]);
    theConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'visualizedConeAperture', 'lightCollectingArea5sigma', ...
        'visualizedConeApertureThetaSamples', 12, ...
        'fontSize', 24, ...
        'activation', reshape(theConeMosaicActivation, [1 1 numel(theConeMosaicActivation)]), ...
        'plotTitle', ' ');

end


function [theScene, theBackgroundScene] = ISETBioSceneFromGrayScaleImageOnDisplay(...
    wavelengthSupport, spatialSupportXDegs, spatialSupportYDegs, theSampleImage)

    % Generate presentation display
    stimulusResolutionDegs = 1/100;
    viewingDistanceMeters = 1.0;
    thePresentationDisplay = visualStimulusGenerator.presentationDisplay(...
            wavelengthSupport, ...
            stimulusResolutionDegs, ...
            viewingDistanceMeters, ...
            'displayType', 'CRT-Sony-HorwitzLab', ...
            'bitDepth', 20, ...
            'meanLuminanceCdPerM2', 100, ...
            'luminanceHeadroom', 0.5);

    grayScaleImage = double(theSampleImage);
    grayScaleImage = grayScaleImage / max(grayScaleImage(:));
    RGBimage = repmat(grayScaleImage, [1 1 3]);

    % Generate a gamma corrected RGB image (RGBsettings) that we can pop in the
    % isetbio scene straightforward
    RGBsettings = (ieLUTLinear(RGBimage, displayGet( thePresentationDisplay, 'inverse gamma'))) / displayGet(thePresentationDisplay, 'nLevels');

    % Generate scene corresponding to the test stimulus on the presentation display
    format = 'rgb';
    meanLuminance = []; % EMPTY, so that mean luminance is determined from the rgb settings values we pass
    theScene = sceneFromFile(flipud(RGBsettings), format, meanLuminance, thePresentationDisplay);

    % Set the desired FOV 
    theScene = sceneSet(theScene, 'h fov', max(spatialSupportXDegs)-min(spatialSupportXDegs));
    
    % Background scene (only used to compute cone modulations from cone excitations)
    RGBimage = RGBimage*0+0.5;
    RGBsettings = (ieLUTLinear(RGBimage, displayGet(thePresentationDisplay, 'inverse gamma'))) / displayGet(thePresentationDisplay, 'nLevels');
    theBackgroundScene = sceneFromFile(flipud(RGBsettings), format, meanLuminance, thePresentationDisplay);
    theBackgroundScene = sceneSet(theBackgroundScene, 'h fov', max(spatialSupportXDegs)-min(spatialSupportXDegs));
end


function [spatialSupportXDegs, spatialSupportYDegs, theSelectedImage, imageHeightDegs] = loadSampleImages(theImageIndex)

    load('fMRIsampleImages', 'sample_images');
    
    imageHeightDegs = 18;
    rowsNum = size(sample_images,1);
    colsNum = size(sample_images,2);

    
    pixelSizeDegs = imageHeightDegs/rowsNum;
    spatialSupportYpixels = 1:rowsNum;
    spatialSupportXpixels = 1:colsNum;
    
    spatialSupportXpixels = spatialSupportXpixels-mean(spatialSupportYpixels);
    spatialSupportYpixels = spatialSupportYpixels-mean(spatialSupportYpixels);
    spatialSupportXDegs = spatialSupportXpixels * pixelSizeDegs;
    spatialSupportYDegs = spatialSupportYpixels * pixelSizeDegs;

    theSelectedImage = sample_images(:,:, theImageIndex);
end

