function testMRGCmosaic
    
    % RGC mosaic ecc and size
    mosaicEccDegs = [5 0]; mosaicSizeDegs = 0.2*[1 1]; whichEye = 'right';
        
    % Chromatic direction examined
    chromaDir = [1.0, 1.0, 0.0]';
    
    % Simulate no Poisson noise in the cone mosaic
    coneNoise = 'none';
    
    % Simulate Gaussian noise in mRGCs with sigma = 0.1 x max(response)
    mRGCnoise = 'random';
    mRGCnoiseFactor = 0.1;
    
    % Save directory
    datasaveDir = '/Volumes/SSDdisk/SpatialTransferData';
    
    % Save filename
    datasaveFile = fullfile(datasaveDir, ...
            sprintf('spatialTransferData_L%2.1f_M_%2.1f_ecc_%2.1f.mat', ...
            chromaDir(1), chromaDir(2), mosaicEccDegs(1)));
        
    % Compute spatial transfer function
    recomputeReponses = ~true;
    
    
    if (recomputeReponses)    
        % Examined spatial frequencies (c/deg)
        examinedSFs = logspace(log10(0.1), log10(50), 16);
        
        % Frame duration for each spatial phase (20 msec)
        % The cone mosaic integration time has to be smaller or equal to this
        frameDurationSeconds = 25/1000;
        coneMosaicIntegrationTime = frameDurationSeconds;
        
        % Specify spatial and temporal parameters of the drifting grating sequence
        presentationMode = 'drifted';

        % Temporal frequency in Hz
        temporalFrequencyHz = 4;    

        % Stimulus duration in seconds
        stimulusDurationSeconds = 0.5;

        % How motion is sampled, 45 degs = 8 spatial phases/period
        spatialPhaseAdvanceDegs = 360 * frameDurationSeconds * temporalFrequencyHz;

        % Generate a drifting grating sequence with 70% contrast
        stimContrast = 0.1;
          
       
         % Stimulus size: cover all of the cone mosaic
        maxEccDegs = max(mosaicEccDegs) + max(0.5*mosaicEccDegs);
        extraDegsForRGCSurround = 2.0 * ...
            RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(maxEccDegs);
        stimFOVdegs = max(mosaicSizeDegs) + extraDegsForRGCSurround;
        
    
        % Compute responses
        instancesNum = 64;
        
        % Generate a midget RGC mosaic with an integration time equal to
        % the frame duration. Shorter integration times are possible, but
        % will take longer to compute. Longer integration times are not
        % appropriate.
        theMidgetRGCmosaic = mRGCmosaic(mosaicEccDegs, mosaicSizeDegs, whichEye, ...
            'coneSpecificityLevel', 100, ...
            'viewTesselationMaps', ~true, ...
            'coneMosaicIntegrationTime', coneMosaicIntegrationTime);

        % Retrieve the coneMosaic that provides input to the mRGC mosaic
        theConeMosaic = theMidgetRGCmosaic.inputConeMosaic;
        
        % Set noise flags for the cone mosaic and for the mRGC mosaic
        % No Poisson noise in the cone responses
        theConeMosaic.noiseFlag = coneNoise;
        
        % Gaussian noise in the mRGC responses
        theMidgetRGCmosaic.noiseFlag = mRGCnoise;
        
        % Sigma of Gaussian noise = 0.1 x max(response)
        theMidgetRGCmosaic.noiseFactor = mRGCnoiseFactor;
        
        % Visualize mRGC positions together with imported ecc-varying cone positions
        % and together with cone positions in the equivalent employed reg-hex mosaic
        %theMidgetRGCmosaic.visualizeInputPositions();

        % Visualize connections of cones to mRGC RF centers
        %theMidgetRGCmosaic.visualizeConeMosaicTesselation('degrees');

        %theMidgetRGCmosaic.visualizeSynthesizedParams();
        %
        %theMidgetRGCmosaic.visualizeConeWeights();

        optics = 'polans';
        if (strcmp(optics, 'default'))
            % Generate default optics
            theOptics = oiCreate();
        else
            % Generate Polans optics for subject 10 with pupil set to 3.0 mm
            PolansSubjectID = 10;
            pupilSizeMM = 4.0;

            theOptics = PolansOptics.oiForSubjectAtEccentricity(PolansSubjectID, mosaicEccDegs, ...
                    pupilSizeMM, theConeMosaic.wave, theConeMosaic.micronsPerDegree, ...
                    'wavefrontSpatialSamples', 501, ...
                    'noLCA', ~true, ...
                    'subtractCentralRefraction', true);
        end
        
        for iSF = 1:numel(examinedSFs)
        
            fprintf('Computing spatial frequency data %d of %d\n', iSF, numel(examinedSFs));
   
            % Instantiate a scene engine for drifting gratings
            driftingGratingSceneEngine = createGratingScene(chromaDir, examinedSFs(iSF), ...
                'duration', stimulusDurationSeconds, ...
                'temporalFrequencyHz', temporalFrequencyHz, ...
                'spatialPhaseAdvanceDegs', spatialPhaseAdvanceDegs, ...
                'fovDegs', stimFOVdegs, ...
                'spatialEnvelopeRadiusDegs', stimFOVdegs, ...
                'minPixelsNumPerCycle', 8, ...
                'spatialEnvelope', 'square', ...
                'presentationMode', presentationMode ...
                );
            
            % Generating drifting grating sequence
            [theDriftingGratingSequence, theStimulusTemporalSupportSeconds] = driftingGratingSceneEngine.compute(stimContrast);

            % Visualize the drifting grating sequence
            driftingGratingSceneEngine.visualizeSceneSequence(theDriftingGratingSequence, theStimulusTemporalSupportSeconds);
            
            % Compute the sequence of optical images corresponding to the drifting grating
            framesNum = numel(theDriftingGratingSequence);
            theListOfOpticalImages = cell(1, framesNum);
            for frame = 1:framesNum
                theListOfOpticalImages{frame} = oiCompute(theDriftingGratingSequence{frame}, theOptics);
            end

            % Generate an @oiSequence object from the list of computed optical images
            theOIsequence = oiArbitrarySequence(theListOfOpticalImages, theStimulusTemporalSupportSeconds);
        
            % Zero eye movements
            eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
            emPaths = zeros(instancesNum, eyeMovementsNum, 2);
        
            % Compute cone mosaic responses
            theConeMosaicResponses = theConeMosaic.computeForOISequence(theOIsequence, ...
                'emPaths', emPaths, ...   % the emPaths
                'currentFlag', false ...  % no photocurrent
            );
    
            % Compute the mRGC mosaic responses
            [theMRGCMosaicResponses, temporalSupportSeconds] = theMidgetRGCmosaic.compute( ...
                    theConeMosaicResponses, theConeMosaic.timeAxis);
            
            % Save responses
            neuralPipelineResponses{iSF} = struct(...
                'coneMosaicTimeAxis', single(theConeMosaic.timeAxis), ...
                'coneMosaicResponses', single(theConeMosaicResponses), ...
                'mRGCMosaicTimeAxis', single(temporalSupportSeconds), ...
                'mRGCMosaicResponses', single(theMRGCMosaicResponses) ...
                );
            
            % Visualize responses
%             theMidgetRGCmosaic.visualizeResponses(temporalSupportSeconds, ...
%                 theMRGCMosaicResponses, ...
%                 'stimulusTemporalSupportSeconds', theStimulusTemporalSupportSeconds,...
%                 'stimulusSceneSequence', theDriftingGratingSequence, ...
%                 'opticalSequence', theOIsequence, ...
%                 'coneMosaicResponse', theConeMosaicResponses);
    
        end % iSF
        
        % Save results
        save(datasaveFile, ...
            'theMidgetRGCmosaic', ...
            'examinedSFs', 'temporalFrequencyHz', ...
            'chromaDir', 'stimContrast', ...
            'neuralPipelineResponses', '-v7.3');
    else
        % Load results
        load(datasaveFile, ...
            'theMidgetRGCmosaic', ...
            'examinedSFs', 'temporalFrequencyHz', ...
            'chromaDir', 'stimContrast', ...
            'neuralPipelineResponses');
        % Plot results
        plotSpatialTransferFunction(theMidgetRGCmosaic, neuralPipelineResponses, ...
                examinedSFs, temporalFrequencyHz, ...
                chromaDir, stimContrast);
    end
end

function plotSpatialTransferFunction(theMidgetRGCmosaic, neuralPipelineResponses, ...
                examinedSFs, temporalFrequencyHz, chromaDir, stimContrast)
            
    responseTimeAxis = neuralPipelineResponses{1}.mRGCMosaicTimeAxis;
    rgcsNum = size(neuralPipelineResponses{1}.mRGCMosaicResponses,2);
    
    % Compute response range and mean responses
    minResponse = Inf;
    maxResponse =-Inf;
    meanResponseMatrix = zeros(numel(examinedSFs), rgcsNum, numel(responseTimeAxis));
    
    for iSF = 1:numel(examinedSFs)
        % Mean over iterations
        meanResponses = squeeze(mean(neuralPipelineResponses{iSF}.mRGCMosaicResponses,1));
        minResponse = min([minResponse min(meanResponses(:))]);
        maxResponse = max([maxResponse max(meanResponses(:))]);
        meanResponseMatrix(iSF,:,:) = meanResponses;
    end % iSF
    
    % Rearrange dimensions [cell, sf, time]
    meanResponseMatrix = permute(meanResponseMatrix, [2 1 3]);
    
    % Initial spatial transfer function matrix
    spatialTransferFunctionMatrix = zeros(rgcsNum, numel(examinedSFs), 2);
    
    visualizeSinusoidalFits = ~true;
    if (visualizeSinusoidalFits)
        figure(1); clf;
    end
    parfor iSF = 1:numel(examinedSFs)
        for iRGC = 1:rgcsNum
            cellResponseTimeSeries = squeeze(meanResponseMatrix(iRGC,iSF,:));
            [responseTimeAxisHR, ...
             fittedSinusoid, ...
             spatialTransferFunctionMatrix(iRGC,iSF,:)] = fit.sinusoidToSingleCellResponse(...
                    responseTimeAxis, cellResponseTimeSeries', temporalFrequencyHz);
            if (visualizeSinusoidalFits)
                plot(responseTimeAxis, cellResponseTimeSeries, 'ks-'); hold on;
                plot(responseTimeAxisHR, fittedSinusoid, 'r-');
                set(gca, 'YLim', [minResponse maxResponse]);
                hold off;
                drawnow;
            end
            
        end
    end % iSF
    
    spatialTransferFunctionsGain = squeeze(spatialTransferFunctionMatrix(:,:,1));
    spatialTransferFunctionsPhase = squeeze(spatialTransferFunctionMatrix(:,:,2));
    
    fittedSpatialTransferFunctionData.x = logspace(log10(examinedSFs(1)), log10(examinedSFs(end)), 60);
    [mappedDoGmodelParams, fittedSpatialTransferFunctionData.y] = fit.spatialTransferFunction(...
        examinedSFs, spatialTransferFunctionsGain, fittedSpatialTransferFunctionData.x);
    
    % Visualize the spatial tranfer function gains
    theMidgetRGCmosaic.visualizeResponseMatrix(examinedSFs, spatialTransferFunctionsGain, ...
        'fittedResponses', fittedSpatialTransferFunctionData);
    
    theMidgetRGCmosaic.visualizeCorrespondenceBetweenMappedAndSynthesizedModelParams(mappedDoGmodelParams);
    
end



