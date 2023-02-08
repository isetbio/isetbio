function [midgetRGCresponses, midgetRGCresponseTemporalSupport] = computeResponsesByPoolingConeResponses(...
    obj, coneMosaicResponses, coneMosaicResponseTemporalSupport)
    
    assert(ndims(coneMosaicResponses) == 3, ...
        'midgetRGCMosaic.computeResponsesByPoolingConeResponses: coneMosaicResponses must have 3 dimensions [nTrials x nTimePoints x nCones]');
    
    assert(size(coneMosaicResponses,2) == numel(coneMosaicResponseTemporalSupport), ...
        'midgetRGCMosaic.computeResponsesByPoolingConeResponses: size(coneMosaicResponses,2) (%d) does not equal the length of temporal support (%d)', ...
        size(coneMosaicResponses,2), numel(coneMosaicResponseTemporalSupport));

    nTrials = size(coneMosaicResponses,1);
    nTimePoints = size(coneMosaicResponses,2);
    

    if (isempty(obj.rgcRFcenterConePoolingMatrix))
        mRGCsNum = size(obj.rgcRFcenterConeConnectivityMatrix,2);
    else
        mRGCsNum = size(obj.rgcRFcenterConePoolingMatrix,2);
    end

    % No temporal filtering for now, so same output temporal support same
    % as the input
    midgetRGCresponseTemporalSupport = coneMosaicResponseTemporalSupport;
    nTimePointsForMidgetRGCResponses = numel(midgetRGCresponseTemporalSupport);

    % Allocate memory for the responses
    midgetRGCresponses = zeros(nTrials, nTimePointsForMidgetRGCResponses, mRGCsNum);
    

    if (isempty(obj.rgcRFcenterConePoolingMatrix))
        fprintf(2,'The center and surround cone pooling matrices have not yet been set (no center/surround RF and no overlap) !!\n');
        fprintf(2,'Using the rgcRFcenterConeConnectivityMatrix instead to compute RF center responses only !!\n');

        % Compute the response of each mRGC
        parfor iRGC = 1:mRGCsNum
            % Retrieve the center cone indices & weights
            centerConnectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
            centerConeIndices = find(centerConnectivityVector > 0.0001);
            centerConeWeights = reshape(centerConnectivityVector(centerConeIndices), [1 1 numel(centerConeIndices)]);
    
            % Sum the weighted cone responses pooled by the RF center only
            centerActivation = sum(bsxfun(@times, coneMosaicResponses(1:nTrials,1:nTimePoints,centerConeIndices), centerConeWeights),3);
            responses(:,:,iRGC) = centerActivation;
        end

    else

        % Compute the response of each mRGC
        parfor iRGC = 1:mRGCsNum
            % Retrieve the center cone indices & weights
            centerConnectivityVector = full(squeeze(obj.rgcRFcenterConePoolingMatrix(:, iRGC)));
            centerConeIndices = find(centerConnectivityVector > 0.0001);
            centerConeWeights = reshape(centerConnectivityVector(centerConeIndices), [1 1 numel(centerConeIndices)]);
    
            % Retrieve the surround cone indices & weights
            surroundConnectivityVector = full(squeeze(obj.rgcRFsurroundConePoolingMatrix(:, iRGC)));
            surroundConeIndices = find(surroundConnectivityVector > 0.0001);
            surroundConeWeights = reshape(surroundConnectivityVector(surroundConeIndices), [1 1 numel(surroundConeIndices)]);
    
            % Sum the weighted cone responses pooled by the RF center and RF surround
            centerActivation = sum(bsxfun(@times, coneMosaicResponses(1:nTrials,1:nTimePoints,centerConeIndices), centerConeWeights),3);
            surroundActivation = sum(bsxfun(@times, coneMosaicResponses(1:nTrials,1:nTimePoints,surroundConeIndices), surroundConeWeights),3);
    
            % Temporal filtering of centerActivation and surroundActivation
            % would go here
            if (nTimePoints > 1)
                %for iTrial = 1:nTrials
                %centerActivation(iTrial,:) = temporalPooling(squeeze(centerActivation(iTrial,:)), coneMosaicResponseTemporalSupport, centerImpulseResponse);
                %surroundActivation(iTrial,:)  = temporalPooling(squeeze(surroundActivation(iTrial,:)), coneMosaicResponseTemporalSupport, surroundImpulseResponse);
                %end
            end

            % Composite response: centerActivation - surroundActivation
            midgetRGCresponses(:,:,iRGC) = centerActivation - surroundActivation;
        end
    end

end

function response = temporalPooling(inputResponse, inputResponseTemporalSupport, impulseResponse)
    %inpulseResponse.temporalSupport
    %impulseResponse.value
end
