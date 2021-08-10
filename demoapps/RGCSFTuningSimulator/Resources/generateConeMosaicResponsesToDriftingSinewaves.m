function generateConeMosaicResponsesToDriftingSinewaves
% Generate the data used by the RGCSFTuningSimulator
%

% History:
%    08/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

    % Eccentricities examined
    mosaicEccentricities = [...
        0 0; ...
        0.2 0; ...
        1.0 0];
    
    % Corresponding mosaic sizes
    mosaicSizes = [ ...
        0.4; ...
        0.4; ...
        0.5];
    
    % SFs examined
    examinedSFs = [0.5 1 2 4 8 12 16 20 25 35 60 90];
    
    % Visualization params
    visualizeStimulusSequence = ~true;
    dataSourceDir = strrep(isetRootPath, 'isettools','demoapps/RGCSFTuningSimulator/Resources');
    
    % Action !
    for mosaicEccIndex = 1:1 % size(mosaicEccentricities,1)
        
        mosaicEcc = mosaicEccentricities(mosaicEccIndex,:);
        mosaicSizeDegs = mosaicSizes(mosaicEccIndex) * [1 1];
        
        % Generate mosaic
        theConeMosaic = [];
        coneMosaicSpatiotemporalActivation = single([]);
        
        for subjectRankOrder = -1:10
            % Generate optics
            theOptics = [];
            
            for iSF = 1:numel(examinedSFs)

                spatialFrequency = examinedSFs(iSF);
                if (iSF == 1)
                    [coneMosaicSpatiotemporalActivation(iSF,:,:), baserateActivation, ...
                        temporalSupportSeconds, theConeMosaic, theOptics] = ...
                        doIt(theConeMosaic, theOptics, subjectRankOrder, spatialFrequency, ...
                        mosaicSizeDegs, mosaicEcc, true, visualizeStimulusSequence);
                else
                    coneMosaicSpatiotemporalActivation(iSF,:,:) = ...
                        doIt(theConeMosaic, theOptics, subjectRankOrder, spatialFrequency, ...
                        mosaicSizeDegs, mosaicEcc, false, visualizeStimulusSequence);
                end

            end

            % Extract the PSF data at 550
            optics = oiGet(theOptics, 'optics');
            waves = opticsGet(optics, 'wave');
            [~,idx] = min(abs(550-waves));
            targetWavelength = waves(idx);
            psfData = struct('psf', [], 'psfSupportDegsX', [], 'psfSupportDegsY', []);
            
            psfData.psf = opticsGet(optics,'psf data',targetWavelength);
            psfData.psf  = psfData.psf  / max(psfData.psf (:));
            
            psfSupportMicrons = opticsGet(optics,'psf support','um');
            psfSupportDegsX = psfSupportMicrons{1}/theConeMosaic.micronsPerDegree;
            psfSupportDegsY = psfSupportMicrons{2}/theConeMosaic.micronsPerDegree;
            psfData.psfSupportDegsX = psfSupportDegsX(1,:) + theConeMosaic.eccentricityDegs(1);
            psfData.psfSupportDegsY = psfSupportDegsY(:,1) + theConeMosaic.eccentricityDegs(2);
            
            % EXPORT
            dataSourceFile = sprintf('PolansSubject%d_ecc_%2.1f_%2.1f_ConeMosaicDriftingSinewaveResponses.mat', subjectRankOrder, mosaicEcc(1), mosaicEcc(2));
            save(fullfile(dataSourceDir, dataSourceFile), ...
                'examinedSFs', ...
                'coneMosaicSpatiotemporalActivation', 'baserateActivation', ...
                'temporalSupportSeconds', 'theConeMosaic', 'psfData', '-v7.3');
        end
    end

end


function [coneMosaicSpatiotemporalActivation, baserateActivation, ...
        temporalSupportSeconds, theConeMosaic, theOptics] = ...
        doIt(theConeMosaic, theOptics, subjectRankOrder, spatialFrequency, mosaicSizeDegs, mosaicEcc, ...
        computeBaseRateActivation, visualizeStimulusSequence)

    % Stimulus params
    meanLuminancdCdPerM2 = 100;
    % This mean chromaticity maximizes the achievable acrhomatic contrast
    meanChromaticityXY = [0.31 0.34];
    chromaDir = [1 1 1];
    stimContrast = 1.0;

    % The stimulus size in degrees (large enough to cover the RGC RF surrounds
    % at the mosaic's max eccentricity)
    maxMosaicEcc = sqrt(sum(mosaicEcc.^2,2)) + 0.5*max(mosaicSizeDegs);
    extraDegsForRGCSurround = 2.0 * ...
          RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(maxMosaicEcc);
    stimFOVdegs = max(mosaicSizeDegs) + 2.0*extraDegsForRGCSurround;

    % The stimulus presentation mode            
    presentationMode = 'drifted';

    % Frame duration for each spatial phase (20 msec)
    % The cone mosaic integration time has to be smaller or equal to this
    displayFrameDurationSeconds = 25/1000;

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

    if (computeBaseRateActivation)
    	% Generate the null stimulus (zero contrast)
        theNullStimulusSequence  = driftingGratingSceneEngine.compute(0);
        theNullStimulus = theNullStimulusSequence{1};
        clear 'theNullStimulusSequence';
    end
    

    if (visualizeStimulusSequence)
        % Visualize the drifting grating sequence
        driftingGratingSceneEngine.visualizeSceneSequence(...
            theDriftingGratingSequence, theStimulusTemporalSupportSeconds);
    end    

    if (isempty(theConeMosaic))
        % Generate mosaic centered at target eccentricity
        fprintf('\t Computing mosaic\n');
        theConeMosaic = cMosaic(...
                'sizeDegs', mosaicSizeDegs, ...     % SIZE in degs
                'eccentricityDegs', mosaicEcc, ...  % ECC in degs
                'opticalImagePositionDegs', 'mosaic-centered', ...
                'integrationTime', displayFrameDurationSeconds...
                );
    end
    
    if (isempty(theOptics))
        fprintf('\t Computing optics\n');
        pupilDiameterMM = 3.0;
        if (subjectRankOrder == -1)
            theOI = oiCreate('wvf human', pupilDiameterMM,[],[], theConeMosaic.micronsPerDegree);
            theOptics = ptb.oiSetPtbOptics(theOI,'opticsModel', 'DeltaFunction');
        else
            %% Select ranking of displayed subject
            if (subjectRankOrder > 0)
                rankedSujectIDs = PolansOptics.constants.subjectRanking;
                testSubjectID = rankedSujectIDs(subjectRankOrder);
                subtractCentralRefraction = ...
                    PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);
            else
                testSubjectID = 0;
                subtractCentralRefraction = false;
            end

            % Generate optics appropriate for the mosaic's eccentricity
            oiEnsemble = theConeMosaic.oiEnsembleGenerate(mosaicEcc, ...
                    'zernikeDataBase', 'Polans2015', ...
                    'subjectID', testSubjectID, ...
                    'pupilDiameterMM', pupilDiameterMM, ...
                    'subtractCentralRefraction', subtractCentralRefraction);

            % Extract the optics
            theOptics = oiEnsemble{1};
        end
    end
    
    
    if (computeBaseRateActivation)
        fprintf('Computing cone mosaic response to the NULL stimulus\n');
        % Compute the null optical image
        theNullOpticalImage = oiCompute(theNullStimulus, theOptics);
        % Compute the null stimulus cone mosaic activation (baserateActivation)
        baserateActivation = single(theConeMosaic.compute(theNullOpticalImage));
    else
        baserateActivation = [];
    end
    
    fprintf('Computing cone mosaic response to %2.1f c/deg\n', spatialFrequency);

    % Compute the sequence of optical images corresponding to the drifting grating
    framesNum = numel(theDriftingGratingSequence);
    theListOfOpticalImages = cell(1, framesNum);
    for frame = 1:framesNum
        theListOfOpticalImages{frame} = ...
            oiCompute(theDriftingGratingSequence{frame}, theOptics);
    end

    
    % Generate an @oiSequence object from the list of computed optical images
    theOIsequence = oiArbitrarySequence(theListOfOpticalImages, theStimulusTemporalSupportSeconds);

    % Compute the spatiotemporal cone-mosaic activation
    [coneMosaicSpatiotemporalActivation, ~, ~, ~, temporalSupportSeconds] = ...
        theConeMosaic.compute(theOIsequence);
    
    % Single precision to save space
    coneMosaicSpatiotemporalActivation = single(coneMosaicSpatiotemporalActivation);
end


    