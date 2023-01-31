function computeMosaicLMnonOpponentSTFsFromConeMosaicLMnonOpponentSTFs(mosaicCenterParams, rfModelParams, opticsParams, useParfor)

    % Generate the frozen mosaic filename
    frozenMosaicFileName = midgetRGCMosaicInspector.frozenMosaicFileName(...
        mosaicCenterParams, rfModelParams.H1cellIndex, opticsParams);

    % Load the frozen midget RGC mosaic
    load(frozenMosaicFileName, 'theMidgetRGCmosaic');

    % Ask the user which optics position to use for the computation
    opticsPositionDegs = midgetRGCMosaicInspector.selectOpticsPosition(theMidgetRGCmosaic);

    % Generate the responses filename
    responsesFileName = midgetRGCMosaicInspector.responsesFileName(...
        frozenMosaicFileName, opticsPositionDegs);
           

    % Generate the input cone mosaics responses filename
    coneMosaicResponsesFileName = midgetRGCMosaicInspector.coneMosaicResponsesFileName(...
        frozenMosaicFileName, opticsPositionDegs);

    % Load cone mosaic responses from disk
    load(coneMosaicResponsesFileName, 'theConeMosaicResponses', 'coneMosaicNullResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts', 'opticsPositionDegs');

    % Allocate memory
    rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsMicrons,1);

    % Single precision responses
    theMidgetRGCMosaicResponses = zeros(...
        numel(orientationsTested), ...
        numel(spatialFrequenciesTested), ...
        numel(spatialPhasesDegs), ...
        rgcsNum, ...
        'single');
   

    disp('Allocated memory');

    % Initialize theInputDataStruct for computation with input cone mosaic responses
    theInputDataStruct = midgetRGCMosaic.inputDataStruct(midgetRGCMosaic.CONE_MOSAIC_RESPONSE_COMPUTE_INPUT_DATA_TYPE);

    % Pre-compute the normalizing cone mosaic responses 
    % Find any cones with zero null response
    coneMosaicNullResponses = squeeze(coneMosaicNullResponses);
    coneIndicesWithZeroNullResponse = find(coneMosaicNullResponses == 0);

    normalizingResponses = 1./coneMosaicNullResponses;
    normalizingResponses(coneIndicesWithZeroNullResponse) = 0;
    normalizingResponses = reshape(normalizingResponses, [1 1 1 numel(normalizingResponses)]);
    coneMosaicNullResponses = reshape(coneMosaicNullResponses, [1 1 1 numel(coneMosaicNullResponses)]);

    % Modulations from excitations
    theConeMosaicResponses = ...
                bsxfun(@times, bsxfun(@minus, theConeMosaicResponses, coneMosaicNullResponses), normalizingResponses);


    for iOri = 1:numel(orientationsTested)
        orientationDegs = orientationsTested(iOri);

        fprintf('Computing midget RGC mosaic STFs by pooling pre-computed cone mosaic responses to the %d degs orientation patterns.\n', ...
                orientationDegs);

        for iFreq = 1:numel(spatialFrequenciesTested)
            % Retrieve cone mosaic responses for all frames of this stimulus
            theConeMosaicSpatioTemporalResponsesToThisStimulusFrames = squeeze(theConeMosaicResponses(iOri, iFreq,:,:));
            
            % Only 1 trial so reshape to reflect that
            theConeMosaicSpatioTemporalResponsesToThisStimulusFrames = reshape(...
                theConeMosaicSpatioTemporalResponsesToThisStimulusFrames, ...
                [1 size(theConeMosaicSpatioTemporalResponsesToThisStimulusFrames,1) size(theConeMosaicSpatioTemporalResponsesToThisStimulusFrames,2)]);

            % Update theCurrentInputDataStruct
            theCurrentInputDataStruct = theInputDataStruct;
            theCurrentInputDataStruct.inputConeMosaicActivationTemporalSupport = 1:numel(spatialPhasesDegs);
            theCurrentInputDataStruct.inputConeMosaicSpatioTemporalActivation = theConeMosaicSpatioTemporalResponsesToThisStimulusFrames;

            % Compute
            r = theMidgetRGCmosaic.compute(theCurrentInputDataStruct, ...
                    'nTrials', 1);

            % Save responses
            trialNo = 1;
            theMidgetRGCMosaicResponses(iOri, iFreq, :,:) = r(trialNo,:,:);

        end % iFreq
    end % iOri


    % Save response data to disk
    save(responsesFileName, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts', 'opticsPositionDegs', '-v7.3');

    fprintf('Saved computed mRGC responses to %s\n', responsesFileName);
end
