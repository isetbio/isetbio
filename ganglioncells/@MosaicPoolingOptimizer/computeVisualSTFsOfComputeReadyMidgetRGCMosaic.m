function computeVisualSTFsOfComputeReadyMidgetRGCMosaic(...
    computeReadyMosaicFilename, coneMosaicSTFresponsesFilename, ...
    mRGCMosaicSTFresponsesFilename)

    % Load the compute-ready MRGC mosaic
    load(computeReadyMosaicFilename, 'theComputeReadyMRGCmosaic');

    % Load the cone mosaic STF responses 
    load(coneMosaicSTFresponsesFilename, 'theConeMosaicSTFresponses', 'theConeMosaicNullResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts');

    fprintf('Will generate mRGCMosaic STF responses based on coneMosaic STF responses in %s\n', coneMosaicSTFresponsesFilename);
    fprintf('The computed mRGCMosaic STF responses will be saved to %s \n', mRGCMosaicSTFresponsesFilename);

    % Transform cone excitation responses to cone modulation responses
    coneIndicesWithZeroNullResponse = find(theConeMosaicNullResponses == 0);
    normalizingResponses = 1./theConeMosaicNullResponses;
    normalizingResponses(coneIndicesWithZeroNullResponse) = 0;
    normalizingResponses = reshape(normalizingResponses, [1 1 numel(normalizingResponses)]);

    % Compute cone mosaic **modulation** responses from the cone mosaic
    % **excitation** responses
    theConeMosaicModulationSTFresponses = 0 * theConeMosaicSTFresponses;
    for iOri = 1:numel(orientationsTested)
        theConeMosaicModulationSTFresponses(iOri,:,:,:) = ...
            bsxfun(@times, bsxfun(@minus, squeeze(theConeMosaicSTFresponses(iOri,:,:,:)), theConeMosaicNullResponses), ...
                 normalizingResponses);
    end
    clear 'theConeMosaicSTFresponses'


    % Compute the visual Rc from summing cone responses pooled by the RF center
    visualRcDegsEstimates = zeros(1, theComputeReadyMRGCmosaic.rgcsNum);
    parfor iRGC = 1:theComputeReadyMRGCmosaic.rgcsNum

        inputConeIndices = find(squeeze(theComputeReadyMRGCmosaic.rgcRFcenterConePoolingMatrix(:,iRGC ))>0.001);
        inputConeWeights = full(theComputeReadyMRGCmosaic.rgcRFcenterConePoolingMatrix(inputConeIndices, iRGC));
        theRFCenterResponsesAcrossAllOrientationsAndSpatialFrequencies = sum(bsxfun(@times, ...
            theConeMosaicModulationSTFresponses(:,:,:,inputConeIndices), ...
            reshape(inputConeWeights, [1 1 1 numel(inputConeIndices)])),4);

        theRFcenterVisualSTF = MosaicPoolingOptimizer.optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies( ...
               orientationsTested, spatialFrequenciesTested, ...
               theRFCenterResponsesAcrossAllOrientationsAndSpatialFrequencies);
        
        fprintf('Computing estimate of the visual RF center characteristic radius for RGC %d of %d \n', ...
            iRGC, theComputeReadyMRGCmosaic.rgcsNum);

        visualRcDegsEstimates(iRGC) = computeVisualRcDegsEstimate(...
            theComputeReadyMRGCmosaic, inputConeIndices, ...
            theRFcenterVisualSTF, spatialFrequenciesTested);
    end


    fprintf('Computing visual STFs for all RGCs in the mosaic ... \n')
    nTimePoints = size(theConeMosaicModulationSTFresponses,3);
    nTrials = 1;
    nCones = size(theConeMosaicModulationSTFresponses,4);

    frameDurationSeconds = 20/1000;
    theConeMosaicResponseTemporalSupportSeconds = (0:1:(nTimePoints-1))*frameDurationSeconds;
    
    % Set noise flag to none. We dont need noisy responses now
    lastNoiseFlag = theComputeReadyMRGCmosaic.noiseFlag;
    theComputeReadyMRGCmosaic.noiseFlag = 'none';

    for iOri = 1:numel(orientationsTested)
        for iSF = 1:numel(spatialFrequenciesTested)
            theConeMosaicResponse = squeeze(theConeMosaicModulationSTFresponses(iOri, iSF,:,:));
            theConeMosaicResponse = reshape(theConeMosaicResponse, [nTrials nTimePoints nCones]);

            fprintf('Computing mRGCMosaic response to %2.0f degs %2.2f c/degs\n', ...
                orientationsTested(iOri), spatialFrequenciesTested(iSF));

            % Compute !
            [theMRGCMosaicResponse, ~, theMRGCresponseTemporalSupportSeconds] = ...
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
        end % iSF
    end % iOri

    % Set noise flag to what was before.
    theComputeReadyMRGCmosaic.noiseFlag =  lastNoiseFlag;


    % Save the computed mRGC  STF responses 
    save(mRGCMosaicSTFresponsesFilename, ...
        'visualRcDegsEstimates', ...
        'theMRGCMosaicSTFresponses', ...
        'theMRGCresponseTemporalSupportSeconds', ...
        'orientationsTested', 'spatialFrequenciesTested', ...
        'spatialPhasesDegs', 'coneContrasts', '-v7.3');


    fprintf('MRGC mosaic STF responses were saved to %s. \n', mRGCMosaicSTFresponsesFilename);

end


function visualRcDegs = computeVisualRcDegsEstimate(theRGCMosaic, inputConeIndices, ...
    theRFcenterVisualSTF, spatialFrequenciesTested)

    % An estimate of the anatomical RcDegs
    anatomicalRcDegs = sqrt(numel(inputConeIndices)) * ...
                       theRGCMosaic.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
                       mean(theRGCMosaic.inputConeMosaic.coneApertureDiametersDegs(inputConeIndices));

    % Some initial estimate of the visual Rc (this is clearly off as it
    % does not account the effect of optics), but it is better than
    % nothing. The Gaussian STF model is only a 2-parameter model (gain,
    % Rc), so it should be able to find the Rc without getting stuct to a
    % local minimum
    initialRcDegs = anatomicalRcDegs;
    
    multiStartsNumDoGFit = 64;
    % Fit the visual STF with a DoG model
    fittedParamsStruct = MosaicPoolingOptimizer.fitGaussianToSubregionSTF(...
                      spatialFrequenciesTested, ...
                      theRFcenterVisualSTF, ...
                      initialRcDegs, ...
                      [], ...
                      multiStartsNumDoGFit);

    visualRcDegs = fittedParamsStruct.finalValues(2);
end


