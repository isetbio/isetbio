function computeVisualSTFsOfComputeReadyMidgetRGCMosaic(...
    computeReadyMosaicFilename, coneMosaicSTFresponsesFilename, ...
    mRGCMosaicSTFresponsesFilename)

    % Load the compute-ready MRGC mosaic
    load(computeReadyMosaicFilename, 'theComputeReadyMRGCmosaic');

    % Load the cone mosaic STF responses 
    load(coneMosaicSTFresponsesFilename, 'theConeMosaicSTFresponses', 'theConeMosaicNullResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts');

    fprintf('MRGC mosaic STF responses will be saved to %s \n', mRGCMosaicSTFresponsesFilename);

    % Transform cone excitation responses to cone modulation responses
    coneIndicesWithZeroNullResponse = find(theConeMosaicNullResponses== 0);
    normalizingResponses = 1./theConeMosaicNullResponses;
    normalizingResponses(coneIndicesWithZeroNullResponse) = 0;
    normalizingResponses = reshape(normalizingResponses, [1 1 numel(normalizingResponses)]);

    for iOri = 1:numel(orientationsTested)
        theConeMosaicSTFresponses(iOri,:,:,:) = ...
            bsxfun(@times, bsxfun(@minus, squeeze(theConeMosaicSTFresponses(iOri,:,:,:)), theConeMosaicNullResponses), ...
                 normalizingResponses);
    end

    nTimePoints = size(theConeMosaicSTFresponses,3);
    nTrials = 1;
    nCones = size(theConeMosaicSTFresponses,4);

    frameDurationSeconds = 20/1000;
    theConeMosaicResponseTemporalSupportSeconds = (0:1:(nTimePoints-1))*frameDurationSeconds;
    
    for iOri = 1:numel(orientationsTested)
        for iSF = 1:numel(spatialFrequenciesTested)
            theConeMosaicResponse = squeeze(theConeMosaicSTFresponses(iOri, iSF,:,:));
            theConeMosaicResponse = reshape(theConeMosaicResponse, [nTrials nTimePoints nCones]);

            fprintf('Computing mRGCMosaic response to %2.0f degs %2.2f c/degs\n', ...
                orientationsTested(iOri), spatialFrequenciesTested(iSF));

            % Compute !
            [theMRGCMosaicResponse, theMRGCresponseTemporalSupportSeconds] = ...
                theComputeReadyMRGCmosaic.compute(theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds);

            % Save responses across all Orientations x Spatial frequencies
            if (iOri == 1) && (iSF == 1)
                theMRGCMosaicSTFresponses = zeros(...
                    numel(orientationsTested), ...
                    numel(spatialFrequenciesTested), ...
                    numel(theMRGCresponseTemporalSupportSeconds), ...
                    theComputeReadyMRGCmosaic.rgcsNum, 'single');
            end
            theMRGCMosaicSTFresponses(iOri, iSF, :,:) = single(squeeze(theMRGCMosaicResponse(1, :,:)));
        end
    end

    % Save the computed mRGC  STF responses 
    save(mRGCMosaicSTFresponsesFilename, 'theMRGCMosaicSTFresponses', ...
        'theMRGCresponseTemporalSupportSeconds', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts', '-v7.3');


    fprintf('MRGC mosaic STF responses were saved to %s. \n', mRGCMosaicSTFresponsesFilename);

end
