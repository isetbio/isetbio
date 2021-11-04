function t_cMosaicSpatioTemporalStimulation
% Demo computation with an off-axis @cMosaic and optics
%
% Description:
%    Shows how to generate and use the new cone mosaic class, @cMosaic.
%    to compute responses to a spatiotemporal stimulus.
%
% See Also:
%   t_cMosaicBasic
%   t_cMosaicOffAxis
%   t_cMosaicBenchMark

% *************************************************************************
% NOTE: This tutorial requires that ISETBioCSFGenerator 
% (available at https://github.com/isetbio/ISETBioCSFGenerator.git)
% is on the MATLAB path
% *************************************************************************

% History:
%    08/08/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

    % Mosaic size and eccentricity
    mosaicSizeDegs = [1 1]*0.5;
    mosaicEcc = [0 0];
    
    
    % Stimulus params
    meanLuminancdCdPerM2 = 100;
    % This mean chromaticity maximizes the achievable acrhomatic contrast
    meanChromaticityXY = [0.31 0.34];
    chromaDir = [1 1 1];
    stimContrast = 1.0/4;

    % The stimulus size in degrees
    stimFOVdegs = max(mosaicSizeDegs);

    % The stimulus presentation mode            
    presentationMode = 'drifted';

    % Frame duration for each spatial phase (20 msec)
    % The cone mosaic integration time has to be smaller or equal to this
    displayFrameDurationSeconds = 25/1000;

    % Spatial frequency in c/deg
    spatialFrequency = 10;
    
    % Temporal frequency in Hz
    temporalFrequencyHz = 4;

    % Stimulus duration in seconds
    stimulusDurationSeconds = 0.5;

    % How motion is sampled, 45 degs = 8 spatial phases/period
    spatialPhaseAdvanceDegs = 360 * displayFrameDurationSeconds * temporalFrequencyHz;

    % Instantiate a scene engine for drifting gratings
    driftingGratingSceneEngine = createGratingScene(chromaDir(:), spatialFrequency, ...
                    'meanLuminanceCdPerM2', meanLuminancdCdPerM2, ...
                    'meanChromaticityXY', meanChromaticityXY, ...
                    'duration', stimulusDurationSeconds, ...
                    'temporalFrequencyHz', temporalFrequencyHz, ...
                    'spatialPhaseAdvanceDegs', spatialPhaseAdvanceDegs, ...
                    'fovDegs', stimFOVdegs, ...
                    'spatialEnvelopeRadiusDegs', stimFOVdegs, ...
                    'minPixelsNumPerCycle', 8, ...
                    'pixelsNum', 256, ...
                    'spatialEnvelope', 'rect', ...
                    'spectralSupport', 400:10:700, ...
                    'presentationMode', presentationMode ...
                    );

    % Generating drifting grating sequence
    [theDriftingGratingSequence, theStimulusTemporalSupportSeconds] = ...
        driftingGratingSceneEngine.compute(stimContrast);


    % Visualize the stimulus
    driftingGratingSceneEngine.visualizeSceneSequence(...
            theDriftingGratingSequence, theStimulusTemporalSupportSeconds);
  

    % Generate mosaic centered at target eccentricity
    fprintf('\t Computing mosaic\n');
    theConeMosaic = cMosaic(...
                'sizeDegs', mosaicSizeDegs, ...     % SIZE in degs
                'eccentricityDegs', mosaicEcc, ...  % ECC in degs
                'opticalImagePositionDegs', 'mosaic-centered', ...
                'integrationTime', displayFrameDurationSeconds ...
                );

    
    fprintf('\t Computing optics\n');
    % Select ranking of displayed subject
    rankedSujectIDs = PolansOptics.constants.subjectRanking;
    subjectRankOrder = 4;
    testSubjectID = rankedSujectIDs(subjectRankOrder);
    subtractCentralRefraction = ...
        PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);

    % Generate optics appropriate for the mosaic's eccentricity
    oiEnsemble = theConeMosaic.oiEnsembleGenerate(mosaicEcc, ...
            'zernikeDataBase', 'Polans2015', ...
            'subjectID', testSubjectID, ...
            'pupilDiameterMM', 3.0, ...
            'subtractCentralRefraction', subtractCentralRefraction);

    % Extract the optics
    theOptics = oiEnsemble{1};
   

    % Compute the sequence of optical images corresponding to the drifting grating
    fprintf('Computing the optical image sequence');
    
    framesNum = numel(theDriftingGratingSequence);
    theListOfOpticalImages = cell(1, framesNum);
    for frame = 1:framesNum
        theListOfOpticalImages{frame} = ...
            oiCompute(theDriftingGratingSequence{frame}, theOptics);
    end
    % Generate an @oiSequence object from the list of computed optical images
    theOIsequence = oiArbitrarySequence(theListOfOpticalImages, theStimulusTemporalSupportSeconds);


    % Compute the spatiotemporal cone-mosaic activation (mean response + 4
    % noisy response instances)
    fprintf('Computing cone mosaic response to %2.1f c/deg\n', spatialFrequency);
    [coneMosaicSpatiotemporalActivation, noisyAbsorptionInstances, ~, ~, responseTimeAxis] = ...
        theConeMosaic.compute(theOIsequence, 'nTrials', 4);

    % Visualize the mean spatiotemporal response of the cone mosaic and a single
    % noisy spatiotemporal instance
    hFig = figure(1);
    set(hFig, 'Position', [10 10 600 1200]);
    ax1 = subplot(2,1,1);
    ax2 = subplot(2,1,2);
    for k = 1:size(noisyAbsorptionInstances,1)
        for t = 1:size(noisyAbsorptionInstances,2)
            theConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax1, ...
                'activation', coneMosaicSpatiotemporalActivation(1,t,:), ...
                'plotTitle', sprintf('noise-free response\n(t = %2.0fmsec)', 1000*responseTimeAxis(t)));
            theConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax2, ...
                'activation', squeeze(noisyAbsorptionInstances(k,t,:)), ...
                'plotTitle', sprintf('noisy response instance\n(t = %2.0fmsec)', 1000*responseTimeAxis(t)));
            
        end
    end


    