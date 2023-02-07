%% Introduction to the midget RGC mosaic (mRGCMosaic) object.
%
% Description:
%    Demonstrates
%        - creation o a midget RGC mosaic,
%        - how to compute with it, and
%        - how to visualize different aspects of the mRGCMosaic
%        - how to visualize its response
%


% History:
%    01/27/23  NPC  ISETBIO Team, Copyright 2023 Wrote it.

function t_mRGCMosaicBasic

%% Load the source midgetRGCMosaic from the database of prebaked mRGCMosaics
% Choose the prebaked mRGCMosaic that was generated at an eccentricity of (0,0) 
% extending over a 3x3 deg, with cone pooling weights tuned for the Polans subject
% with rank 6 and a pupil size of 3.0 mm

mosaicCenterParams = struct(...
    'positionDegs',[0 0], ...
    'sizeDegs',  [3 3], ...        
    'whichEye', 'right eye');

opticsParams = struct(...
            'ZernikeDataBase', 'Polans2015', ...
            'subjectRankOrder', 6, ...
            'pupilDiameterMM', 3.0 ...
        );

rfModelParams = struct(...
    'H1cellIndex', 1 ...
    );

% Load it
theSourceMidgetRGCMosaic = loadSourceMidgetRGCMosaic(mosaicCenterParams, rfModelParams, opticsParams);

% The SourceMidgetRGCMosaic contains all the information that was used to
% derive the weights. It is not compute-ready.

% Instantiate a compute-ready mRGCMosaic from the sourceMidgetRGCMosaic
% Here we are using part of the sourceMidgetRGCMosaic, centered at (x,y) = (1,0.5), 
% with width = 0.4 degs and height = 0.2 degs

theMRGCMosaic = mRGCMosaic(theSourceMidgetRGCMosaic, ...
        'eccentricityDegs', [1 0.5], ...
        'sizeDegs', [0.5 0.5], ...
        'name', 'my small off-center mRGC mosaic', ...
        'beVerbose', true, ...
        'visualizeSpatialRelationshipToSourceMosaic', true);

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

function theSourceMidgetRGCMosaic = loadSourceMidgetRGCMosaic(mosaicCenterParams, rfModelParams, opticsParams)

    dropboxDir = midgetRGCMosaicInspector.localDropboxPath;
    frozenMidgetRGCMosaicsDir = 'productionMidgetRGCMosaics/frozenMosaics';
    directoryPath = fullfile(dropboxDir,frozenMidgetRGCMosaicsDir);

    
    sourceMidgetRGCMosaicFileName = ...
        sprintf('MRGCmosaic_Ecc_%2.1f_%2.1f_sizeDegs_%2.1f_%2.1f_%s_Rank_%d_Pupil_%2.1f_H1cellIndex_%d_Frozen.mat', ...
        mosaicCenterParams.positionDegs(1), mosaicCenterParams.positionDegs(2), ...
        mosaicCenterParams.sizeDegs(1), mosaicCenterParams.sizeDegs(2), ...
        opticsParams.ZernikeDataBase, ...
        opticsParams.subjectRankOrder, ...
        opticsParams.pupilDiameterMM, ...
        rfModelParams.H1cellIndex);

    theFileName = fullfile(directoryPath, sourceMidgetRGCMosaicFileName);
    fprintf('Will try to load %s ... \n', sourceMidgetRGCMosaicFileName)
    fprintf('from %s ... \n', directoryPath);

    % Check that the mosaic directory exists
    assert(isfolder(directoryPath), sprintf('Mosaic directory (''%s'') not found.', directoryPath));

    % Check that the mosaic file exists
    assert(isfile(theFileName), sprintf('Mosaic file (''%s'') not found.', theFileName));

    % Mosaic file found, so load the data
    load(theFileName, 'theMidgetRGCmosaic');
    theSourceMidgetRGCMosaic = theMidgetRGCmosaic;
    clear 'theMidgetRGCmosaic';
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
            theOI  = oiCompute(theNullStimulusScene, theOI);

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
