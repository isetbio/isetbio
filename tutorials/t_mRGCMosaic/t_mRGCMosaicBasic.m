%% Introduction to the midget RGC mosaic (mRGCMosaic) object.
%
% Description:
%    Demonstrates
%        - creation of a midget RGC mosaic,
%        - how to compute with it, and
%        - how to visualize different aspects of the mRGCMosaic
%        - how to visualize its response
%


% History:
%    01/27/23  NPC  ISETBIO Team, Copyright 2023 Wrote it.

function t_mRGCMosaicBasic

%% Specify params of the source computeReadyMosaic 

mosaicParams = struct(...
    'eccDegs', [0 0], ...
    'sizeDegs',  [2 2], ...
    'rgcType', 'ONcenterMidgetRGC');

opticsParams = struct(...
            'ZernikeDataBase', 'Polans2015', ...
            'examinedSubjectRankOrder', 6, ...
            'pupilDiameterMM', 3.0, ...
            'analyzedEye', 'right eye', ...
            'refractiveErrorDiopters', 0.0, ...
            'positionDegs', [] ...
        );

retinalRFmodelParams = struct(...
    'conePoolingModel', 'arbitraryCenterConeWeights_doubleExpH1cellIndex4SurroundWeights' ...
    );

% Load the source compute-ready mosaic
theComputeReadyRGCMosaic = loadSourceComputeReadyRGCMosaic(mosaicParams, opticsParams, retinalRFmodelParams);

% Generate a smaller compute-ready mosaic by cropping the source
% compute-ready mosaic
eccentricityDegs = [0.65 0.25];
sizeDegs = [0.4 0.4];

% Instantiate a compute-ready mRGCMosaic from the sourceMidgetRGCMosaic
% Here we are using part of the sourceMidgetRGCMosaic, centered at (x,y) = (1,0.5), 
% with width = 0.4 degs and height = 0.2 degs

domainVisualizationLimits = [...
    mosaicParams.eccDegs(1)-0.5*mosaicParams.sizeDegs(1) ...
    mosaicParams.eccDegs(1)+0.5*mosaicParams.sizeDegs(1) ...
    mosaicParams.eccDegs(2)-0.5*mosaicParams.sizeDegs(2) ...
    mosaicParams.eccDegs(2)+0.5*mosaicParams.sizeDegs(2) ...
    ];

domainVisualizationTicks = struct(...
    'x', linspace(mosaicParams.eccDegs(1)-0.5*mosaicParams.sizeDegs(1), mosaicParams.eccDegs(1)+0.5*mosaicParams.sizeDegs(1), 3), ...
    'y', linspace(mosaicParams.eccDegs(2)-0.5*mosaicParams.sizeDegs(2), mosaicParams.eccDegs(2)+0.5*mosaicParams.sizeDegs(2), 3));

hFig = figure(1); clf;
ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);

theComputeReadyRGCMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax1, ...
    'domainVisualizationLimits', domainVisualizationLimits, ...
    'domainVisualizationTicks', domainVisualizationTicks, ...
    'identifyInputCones', true, ...
    'identifyPooledCones', true, ...
    'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
    'plotTitle', 'original');

theComputeReadyRGCMosaic.cropToSizeAtEccentricity(...
        sizeDegs, eccentricityDegs, ...
        'name', 'my small off-center mRGC mosaic', ...
        'beVerbose', true, ...
        'visualizeSpatialRelationshipToSourceMosaic', true);


theComputeReadyRGCMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax2, ...
    'domainVisualizationLimits', domainVisualizationLimits, ...
    'domainVisualizationTicks', domainVisualizationTicks, ...
    'identifyInputCones', true, ...
    'identifyPooledCones', true, ...
    'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
    'plotTitle', sprintf('cropped at %2.2f,%2.2f with size %2.2f,%2.2f', eccentricityDegs(1), eccentricityDegs(2), sizeDegs(1), sizeDegs(2)));

pause;


% Retrive the native optics for this mRGCmosaic
[theOI, thePSF, theZcoeffs] = theMRGCMosaic.nativeOptics;


spatialFrequenciesTested = [4]; %[0.25 0.5 1 2 4 6 8 12 16 20 24 32 48 64];
orientationTested = 0;
[theConeMosaicResponses, coneMosaicNullResponse] = computeConeMosaicResponses(...
    spatialFrequenciesTested, orientationTested, ...
    theMRGCMosaic.inputConeMosaic, theOI);


stimFramesNum = size(theConeMosaicResponses,2);
nConesNum = size(theConeMosaicResponses,3);

frameDurationSeconds = 1/30;
theConeMosaicResponseTemporalSupportSeconds = (0:(stimFramesNum-1))*frameDurationSeconds;

iTrial = 1; nTrialsNum = 1;
for iSF = 1:numel(spatialFrequenciesTested)
    theSpatioTemporalConeMosaicResponse = squeeze(theConeMosaicResponses(iSF,:,:));
    
    theSpatioTemporalConeMosaicResponse = reshape(theSpatioTemporalConeMosaicResponse, [nTrialsNum stimFramesNum, nConesNum]);
 
    [theMRGCMosaicResponses(iTrial, iSF,:,:), theMRGCresponseTemporalSupportSeconds] = theMRGCMosaic.compute( ...
             theSpatioTemporalConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds);
end

activationRange = [min(theMRGCMosaicResponses(:)) max(theMRGCMosaicResponses(:))];

hFig = figure(1);
clf;
ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);

oneCenterConeRGCindices = find(theMRGCMosaic.centerSubregionConesNums == 1);
twoCenterConeRGCindices = find(theMRGCMosaic.centerSubregionConesNums == 2);
visualizedRGCindices = [oneCenterConeRGCindices(1:2) twoCenterConeRGCindices(1:2)];


% Visualize the centers of all RGCs, identifying 2 1-cone center RGCs and 2, 2-cone center RGCs
theMRGCMosaic.visualize(...,
    'figureHandle', hFig, ...
    'axesHandle', ax1, ...
    'identifyInputCones', true, ...
    'labelRGCsWithIndices', visualizedRGCindices, ...
    'backgroundColor', [0.7 0.7 0.7]);

% Visualize the mRGCmosaic activation, frame-by-frame
for iTrial = 1:nTrialsNum
    for iSF = 1:numel(spatialFrequenciesTested)
    for iTimeBin = 1:numel(theMRGCresponseTemporalSupportSeconds)
        theMRGCMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax2, ...
            'identifyInputCones', true, ...
            'activation', squeeze(theMRGCMosaicResponses(iTrial, iSF, iTimeBin,:)), ...
            'activationRange', activationRange, ...
            'plotTitle', sprintf('%2.2f c/deg, %d msec', spatialFrequenciesTested(iSF), theMRGCresponseTemporalSupportSeconds(iTimeBin)));
        pause
    end
    end
end






theMRGCMosaic.multifocalRTVFgrids
theMRGCMosaic.multifocalRTVFopticsParams

% Visualize the RFs of the identified  RGCs
theMRGCMosaic.visualizeRFs(visualizedRGCindices);

end


% ======== HELPER FUNCTIONS ========

function theSourceRGCMosaic = loadSourceComputeReadyRGCMosaic(mosaicParams, opticsParams, retinalRFmodelParams)

    % Load the compute-ready RGC mosaic
    [computeReadyMosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);
    
    sourceRGCMosaicFileName = fullfile(resourcesDirectory, computeReadyMosaicFileName);

    fprintf('Will try to load %s ... \n', sourceRGCMosaicFileName)
    fprintf('from %s ... \n', resourcesDirectory);

    % Check that the mosaic directory exists
    assert(isfolder(resourcesDirectory), sprintf('Mosaic directory (''%s'') not found.', resourcesDirectory));

    % Check that the mosaic file exists
    assert(isfile(sourceRGCMosaicFileName), sprintf('Mosaic file (''%s'') not found.', sourceRGCMosaicFileName));

    % Mosaic file found, so load the data
    load(sourceRGCMosaicFileName, 'theComputeReadyMRGCmosaic');
    theSourceRGCMosaic = theComputeReadyMRGCmosaic;
    clear 'theComputeReadyMRGCmosaic';
    fprintf('Loaded source mosaic.\n');
end

function [theConeMosaicResponses, coneMosaicNullResponse] = computeConeMosaicResponses(...
    spatialFrequenciesTested, orientationTested, ...
    theInputConeMosaic, theOI)

    fprintf('Computing cone mosaic responses\n')
    [stimParams, thePresentationDisplay] = setupVisualSTFmappingExperiment(...
        spatialFrequenciesTested, orientationTested, ...
        theInputConeMosaic.sizeDegs, ...
        theInputConeMosaic.wave);

    % Allocate memory
    conesNum = size(theInputConeMosaic.coneRFpositionsDegs,1);

    % Single precision responses
    theConeMosaicResponses = zeros(...
        numel(stimParams.spatialFrequenciesTested), ...
        numel(stimParams.spatialPhasesDegs), ...
        conesNum, ...
        'single');

    % Empty the noiseFreeAbsorptionsCountNull
    noiseFreeAbsorptionsCountNull = [];

    for iFreq = 1:numel(stimParams.spatialFrequenciesTested)
        theCurrentStimParams = stimParams;
        theCurrentStimParams.spatialFrequencyCPD = stimParams.spatialFrequenciesTested(iFreq);
                
        % Generate spatial modulation patterns for each stimulus frame
        theDriftingGratingSpatialModulationPatterns = rfMappingStimulusGenerator.driftingGratingFrames(theCurrentStimParams);
    
        % Generate scenes for the different spatial phases
        [theDriftingGratingFrameScenes, theNullStimulusScene] = ...
           rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
                        thePresentationDisplay, theCurrentStimParams, theDriftingGratingSpatialModulationPatterns, ...
                        'validateScenes', false);
    
         % Compute the cone mosaic responses to the null scene OI
         if (isempty(noiseFreeAbsorptionsCountNull))
            % Compute the optical image of the null scene
            theOI = oiCompute(theNullStimulusScene, theOI);

            % Compute the cone mosaic null responses
            coneMosaicNullResponse = theInputConeMosaic.compute(theOI, 'nTrials', 1);
         end

         % Allocate memory for all frames
         theFrameResponses = zeros(numel(theCurrentStimParams.spatialPhasesDegs), conesNum, 'single');

         parfor iFrame = 1:numel(theCurrentStimParams.spatialPhasesDegs)
            % Get scene corresponding to this stimulus frame
            theFrameScene = theDriftingGratingFrameScenes{iFrame};

            % Compute the optical image of the frame scene
            theCurrentOI = oiCompute(theFrameScene, theOI);

            % Compute the cone mosaic responses
            [noiseFreeAbsorptionsCount, noisyAbsorptionsCountInstances] = ...
                        theInputConeMosaic.compute(theCurrentOI, 'nTrials', 1);

            theFrameResponses(iFrame,:) = single(noiseFreeAbsorptionsCount(1,1,:));
         end

         % Save theFrameResponses
         theConeMosaicResponses(iFreq,:,:) = theFrameResponses;

    end % iFreq
    fprintf('Done computing cone mosaic responses\n')
end


function [stimParams, theDisplay] = setupVisualSTFmappingExperiment(...
        spatialFrequenciesTested, orientationTested, ...
        coneMosaicSizeDegs, wavelengthSupport)

    viewingDistanceMeters = 4;
    stimulusPixelsNum = 512*2;
    coneContrasts = [1 1 0];

    % Generate a presentation display with a desired resolution
    sceneFOVdegs = coneMosaicSizeDegs;
    retinalImageResolutionDegs = max(sceneFOVdegs)/stimulusPixelsNum;

    % At least 6 samples / period
    maxSF = 1/(2*3*retinalImageResolutionDegs);
    if (max(spatialFrequenciesTested) > maxSF)
        fprintf('Max SF examined (%2.2f c/deg) is too high for this FOV (%2.2f degs) and pixels num (%d). (SFmax: %2.2f c/deg)\n', ...
            max(spatialFrequenciesTested), max(sceneFOVdegs), stimulusPixelsNum, maxSF);
        idx = find(spatialFrequenciesTested <= maxSF);
        spatialFrequenciesTested = spatialFrequenciesTested(idx);
        if (maxSF > max(spatialFrequenciesTested))
            spatialFrequenciesTested(numel(spatialFrequenciesTested)+1) = maxSF;
        end

        fprintf('Will only measure the STF up to %2.2f c/deg.\n', max(spatialFrequenciesTested));
    end

    stimSizeDegs = max(sceneFOVdegs);
    pixelSizeDegs = retinalImageResolutionDegs;

    theDisplay = rfMappingStimulusGenerator.presentationDisplay(...
            wavelengthSupport, pixelSizeDegs, ...
            viewingDistanceMeters);

    spatialSupportDegs = rfMappingStimulusGenerator.spatialSupport(...
        stimSizeDegs, pixelSizeDegs);
        
    
    % Stim params for the STF mapping
    stimParams = struct(...
            'backgroundLuminanceCdM2', 50.0, ...
            'backgroundChromaticity', [0.301 0.301], ...
            'coneContrasts', coneContrasts, ...
            'contrast', 0.75, ...
            'orientationsTested', orientationTested, ...
            'spatialFrequenciesTested', spatialFrequenciesTested, ...
            'orientationDegs', 0, ...
            'spatialFrequencyCPD', spatialFrequenciesTested(1), ...
            'spatialPhaseIncrementDegs', 30, ...
            'pixelSizeDegs', pixelSizeDegs, ...
            'stimSizeDegs', stimSizeDegs, ...
            'spatialMask', ones(numel(spatialSupportDegs), numel(spatialSupportDegs)), ...
            'wavelengthSupport', displayGet(theDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(theDisplay, 'viewing distance') ...
            );


    [~, stimParams.spatialPhasesDegs] = ...
        rfMappingStimulusGenerator.driftingGratingFrames(stimParams);
   
end
