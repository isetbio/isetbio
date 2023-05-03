% Method to compute the spatiotemporal response of the mRGCMosaic given the response of its input cone
% mosaic
function [noiseFreeMRGCresponses, noisyMRGCresponseInstances, responseTemporalSupportSeconds] = compute(obj, ...
            theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds, varargin)

    p = inputParser;
    p.addParameter('timeResolutionSeconds', [], @(x)(isempty(x))||(isscalar(x)));

    % Parse input
    p.parse(varargin{:});
    timeResolutionSeconds = p.Results.timeResolutionSeconds;
    

    % Parse input
    assert(ndims(theConeMosaicResponse) == 3, ...
        'The theConeMosaicResponse must have 3 dimensions [nTrials x nTimePoints x nCones]');
    
    assert(size(theConeMosaicResponse,2) == numel(theConeMosaicResponseTemporalSupportSeconds), ...
        'The size(theConeMosaicResponsth,2) (%d) does not equal the length of temporal support (%d)', ...
        size(theConeMosaicResponse,2), numel(theConeMosaicResponseTemporalSupportSeconds));

    if (numel(theConeMosaicResponseTemporalSupportSeconds)>1)
        inputTimeResolutionSeconds = theConeMosaicResponseTemporalSupportSeconds(2)-theConeMosaicResponseTemporalSupportSeconds(1);
    else
        inputTimeResolutionSeconds = 1.0;
    end

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
    nTrials = size(theConeMosaicResponse,1);
    noiseFreeMRGCresponses = zeros(nTrials, numel(responseTemporalSupportSeconds), obj.rgcsNum);
    
    % Delta function center impulse response with a length of 200 mseconds
    theImpulseResponseTemporalSupport = 0:timeResolutionSeconds:0.2;
    theRFcenterImpulseResponse = generateCenterTemporalImpulseResponse(theImpulseResponseTemporalSupport);
    theRFsurroundImpulseResponse = generateSurroundTemporalImpulseResponse(theImpulseResponseTemporalSupport);


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

        % Composite response gain: inversely proportional to sum(centerWeights)
        responseGain = 1.0 / sum(centerConeWeights);

        % Composite respose
        noiseFreeMRGCresponses(:,:,iRGC) = responseGain * (centerSpatiallyIntegratedActivations - surroundSpatiallyIntegratedActivations);
    end % parfor

    % vMembrane additive noise
    fprintf('Computing noisy mRGC response instances with vMembrane noise std = %2.3f\n', obj.vMembraneGaussianNoiseSigma);
    noisyMRGCresponseInstances = noiseFreeMRGCresponses + ...
            obj.vMembraneGaussianNoiseSigma * randn(size(noiseFreeMRGCresponses));
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