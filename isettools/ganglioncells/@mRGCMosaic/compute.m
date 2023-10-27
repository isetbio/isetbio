% Method to compute the spatiotemporal response of the mRGCMosaic given the response of its input cone
% mosaic
function [noiseFreeMRGCresponses, noisyMRGCresponseInstances, responseTemporalSupportSeconds] = compute(obj, ...
            theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds, varargin)

    p = inputParser;
    p.addParameter('nTrials', [], @isscalar);
    p.addParameter('timeResolutionSeconds', [], @(x)(isempty(x))||(isscalar(x)));
    p.addParameter('seed', [], @isnumeric);

    % Parse input
    p.parse(varargin{:});
    mRGCMosaicNoisyResponseInstancesNum = p.Results.nTrials;
    timeResolutionSeconds = p.Results.timeResolutionSeconds;
    noiseSeed = p.Results.seed;

    % Validate input: ensure theConeMosaicResponse is a 3D matrix
    assert(ndims(theConeMosaicResponse) == 3, ...
        'The theConeMosaicResponse must have 3 dimensions [nTrials x nTimePoints x nCones]');
    
    % Validate input: ensure that the temporal dimensions of the cone
    % mosaic response and its temporal support match
    assert(size(theConeMosaicResponse,2) == numel(theConeMosaicResponseTemporalSupportSeconds), ...
        'The size(theConeMosaicResponsth,2) (%d) does not equal the length of temporal support (%d)', ...
        size(theConeMosaicResponse,2), numel(theConeMosaicResponseTemporalSupportSeconds));


    if (numel(theConeMosaicResponseTemporalSupportSeconds)>1)
        inputTimeResolutionSeconds = theConeMosaicResponseTemporalSupportSeconds(2)-theConeMosaicResponseTemporalSupportSeconds(1);
    else
        inputTimeResolutionSeconds = 1.0;
    end

    % Find out the # of trials to compute
    coneMosaicResponseInstancesNum = size(theConeMosaicResponse,1);
    if (isempty(mRGCMosaicNoisyResponseInstancesNum))
        % No 'nTrials' passed, so we will create as many instances as there
        % are cone mosaic noisy intances
        mRGCMosaicNoisyResponseInstancesNum = coneMosaicResponseInstancesNum;
    else
        % 'nTrials' for mRGCMosaic responses is passed
        if (coneMosaicResponseInstancesNum > 1)
            % We also have more than 1 input cone mosaic response.
            % Must be noisy cone mosaic responses instances. 
            % Ensure that mRGCMosaicNoisyResponseInstancesNum == coneMosaicResponseInstancesNum
            assert(mRGCMosaicNoisyResponseInstancesNum == coneMosaicResponseInstancesNum, ...
               '''nTrials'' (%d) does not match the trials of the input cone mosaic response (%d)', ...
               mRGCMosaicNoisyResponseInstancesNum,coneMosaicResponseInstancesNum);
        else
            % A single trial of input cone mosaic response, must be the noise-free cone mosaic response. Replicate it 
            theConeMosaicResponse = repmat(theConeMosaicResponse, [mRGCMosaicNoisyResponseInstancesNum 1 1]);
        end
    end
    nTrials = mRGCMosaicNoisyResponseInstancesNum;


    inputTimePoints = numel(theConeMosaicResponseTemporalSupportSeconds);
    if (isempty(timeResolutionSeconds))
        timeResolutionSeconds = inputTimeResolutionSeconds;
        responseTemporalSupportSeconds = theConeMosaicResponseTemporalSupportSeconds;
    else
        responseTemporalSupportSeconds = - inputTimeResolutionSeconds/2 + ...
                                                (responseTemporalSupportSeconds(1)) : ...
                                                timeResolutionSeconds : ...
                                                (responseTemporalSupportSeconds(end)+inputTimeResolutionSeconds/2);
    end

   
    % Allocate memory for the computed responses
    noiseFreeMRGCresponses = zeros(nTrials, numel(responseTemporalSupportSeconds), obj.rgcsNum);
    
    % Delta function center impulse response with a length of 200 mseconds
    theImpulseResponseTemporalSupport = 0:timeResolutionSeconds:0.2;
    theRFcenterImpulseResponse = generateCenterTemporalImpulseResponse(theImpulseResponseTemporalSupport);
    theRFsurroundImpulseResponse = generateSurroundTemporalImpulseResponse(theImpulseResponseTemporalSupport);

    if (isempty(obj.rgcRFgains))
        fprintf(2,'the mRGCMosaic.rgcRFgains property has not been set: will employ the ''1/integrated center cone weights'' method\n');
    end

    % Compute the response of each mRGC
    parfor iRGC = 1:obj.rgcsNum
        % Retrieve the center cone indices & weights
        centerConnectivityVector = full(squeeze(obj.rgcRFcenterConePoolingMatrix(:, iRGC)));
        centerConeIndices = find(centerConnectivityVector > 0.0001);
        centerConeWeights = reshape(centerConnectivityVector(centerConeIndices), [1 1 numel(centerConeIndices)]);
        
        % Retrieve the surround cone indices & weights
        surroundConnectivityVector = full(squeeze(obj.rgcRFsurroundConePoolingMatrix (:, iRGC)));
        surroundConeIndices = find(surroundConnectivityVector > 0.0001);
        surroundConeWeights = reshape(surroundConnectivityVector(surroundConeIndices), [1 1 numel(surroundConeIndices)]);

        % Spatially pool the weighted cone responses to the RF center
        centerSpatiallyIntegratedActivations = sum(bsxfun(@times, theConeMosaicResponse(1:nTrials,1:inputTimePoints,centerConeIndices), centerConeWeights),3);

        % Spatially pool the weighted cone responses to the RF surround
        surroundSpatiallyIntegratedActivations = sum(bsxfun(@times, theConeMosaicResponse(1:nTrials,1:inputTimePoints, surroundConeIndices), surroundConeWeights),3);

        if (numel(theConeMosaicResponseTemporalSupportSeconds)>1)
            % Temporally filter the center responses
            centerSpatiallyIntegratedActivations = temporalFilter(centerSpatiallyIntegratedActivations, ...
                theConeMosaicResponseTemporalSupportSeconds, ...
                responseTemporalSupportSeconds, ...
                theRFcenterImpulseResponse);
    
            % Temporally filter the surround responses
            surroundSpatiallyIntegratedActivations = temporalFilter(surroundSpatiallyIntegratedActivations, ...
                theConeMosaicResponseTemporalSupportSeconds, ...
                responseTemporalSupportSeconds, ...
                theRFsurroundImpulseResponse);
        end

        % Response gain
        if (isempty(obj.rgcRFgains))
            responseGain = 1.0 / sum(centerConeWeights);
        else
            responseGain = obj.rgcRFgains(iRGC);
        end

        % Composite respose
        noiseFreeMRGCresponses(:,:,iRGC) = responseGain * (centerSpatiallyIntegratedActivations - surroundSpatiallyIntegratedActivations);
    end % parfor

    % Check noiseFlag. If empty, set it to 'random'
    if (isempty(obj.noiseFlag))
        fprintf('Warning: The mRGCMosaic.noiseFlag not set before calling the compute() method. Setting it to ''random''.');
        obj.noiseFlag = 'random';
    end

    if (strcmp(obj.noiseFlag, 'none'))
        noisyMRGCresponseInstances = [];
        return;
    end

    % Generate noisy instances
    noisyMRGCresponseInstances = obj.noisyInstances(noiseFreeMRGCresponses, ...
        'seed', noiseSeed);
end


function theResponses = temporalFilter(theInputResponses, ...
                theInputResponsesTemporalSupportSeconds, ...
                theOutputResponsesTemporalSupportSeconds, ...
                theImpulseResponse)

    [nTrials, nInputTimePoints] = size(theInputResponses);
    nOutputPoints = numel(theOutputResponsesTemporalSupportSeconds);

    % Allocate memory
    theResponses = zeros(nTrials, nOutputPoints);
    
    % Scaling factor to account for difference in binwidth
    theOutputResponseResolution = theOutputResponsesTemporalSupportSeconds(2)-theOutputResponsesTemporalSupportSeconds(1);
    theInputResponseResolution = theInputResponsesTemporalSupportSeconds(2)-theInputResponsesTemporalSupportSeconds(1);
    theScalingFactor = theOutputResponseResolution/theInputResponseResolution;


    for iTrial = 1:nTrials
        theInputResponse = theInputResponses(iTrial,:);

        % Interpolation of input response to same timescale as the impulse
        % response
        theInterpolatedInputResponse = theScalingFactor * interp1(...
            theInputResponsesTemporalSupportSeconds, ...
            theInputResponse, ...
            theOutputResponsesTemporalSupportSeconds, 'nearest', 'extrap');

        % Temporal filter
        theResponse = conv(theInterpolatedInputResponse, theImpulseResponse);

        % Decimate to the # of output time points
        theResponses(iTrial,:) = theResponse(1:nOutputPoints);

        debugFilter = ~true;
        if (debugFilter)
            figure(1); clf;
            
            theImpulseResponseTemporalSupportSeconds = theOutputResponsesTemporalSupportSeconds(1) + (0:(numel(theImpulseResponse)-1))*theOutputResponseResolution;
            plot(theInputResponsesTemporalSupportSeconds, theInputResponse, 'gs', 'MarkerSize', 20);
            hold on;
            plot(theOutputResponsesTemporalSupportSeconds, theInterpolatedInputResponse, 'k.-', 'MarkerSize', 12, 'MarkerFaceColor', [0 0 0]);
            plot(theOutputResponsesTemporalSupportSeconds, theResponses(iTrial,:), 'b-');
            plot(theImpulseResponseTemporalSupportSeconds, theImpulseResponse, 'r.', 'LineWidth', 1.5);
            drawnow
        end
    end
end


function theImpulseResponse = generateCenterTemporalImpulseResponse(temporalSupport)
    theImpulseResponse = zeros(1,numel(temporalSupport));
    theImpulseResponse(1:1) = 1;
    %theImpulseResponse(4:8) = -0.2;
end

function theImpulseResponse = generateSurroundTemporalImpulseResponse(temporalSupport)
    % Just a delay for now
    theImpulseResponse = zeros(1,numel(temporalSupport));
    theImpulseResponse(1:1) = 1;

    %theImpulseResponse(3+(1:2)) = 1;
    %theImpulseResponse(3+(4:8)) = -0.2;
end