% Method to compute the spatiotemporal response of the mRGCMosaic given the response of its input cone mosaic
function [noiseFreeMRGCresponses, noisyMRGCresponseInstances, responseTemporalSupportSeconds] = compute(obj, ...
            theInputConeMosaicResponse, theInputConeMosaicResponseTemporalSupportSeconds, varargin)

    p = inputParser;
    p.addParameter('nTrials', [], @isscalar);
    p.addParameter('timeResolutionSeconds', [], @(x)(isempty(x))||(isscalar(x)));
    p.addParameter('nonLinearitiesList', [], @(x)(isempty(x))||(iscell(x))||(isstruct(x)));
    p.addParameter('seed', [], @isnumeric);


    % Parse input
    p.parse(varargin{:});
    mRGCMosaicNoisyResponseInstancesNum = p.Results.nTrials;
    timeResolutionSeconds = p.Results.timeResolutionSeconds;
    nonLinearitiesList = p.Results.nonLinearitiesList;
    noiseSeed = p.Results.seed;

    % Validate input: ensure theConeMosaicResponse is a 3D matrix
    assert(ndims(theInputConeMosaicResponse) == 3, ...
        'theInputConeMosaicResponse must have 3 dimensions [nTrials x nTimePoints x nCones]');
    
    % Validate input: ensure that the temporal dimensions of the cone
    % mosaic response and its temporal support match
    assert(size(theInputConeMosaicResponse,2) == numel(theInputConeMosaicResponseTemporalSupportSeconds), ...
        'The size(theInputConeMosaicResponse,2) (%d) does not equal the length of temporal support (%d)', ...
        size(theInputConeMosaicResponse,2), numel(theInputConeMosaicResponseTemporalSupportSeconds));


    if (isstruct(nonLinearitiesList))
        theSingleNonLinearity = nonLinearitiesList;
        nonLinearitiesList = cell(1,1)
        nonLinearitiesList{1} = theSingleNonLinearity;
    end

    % Nonlinearities
    if (numel(nonLinearitiesList)>0)
        fprintf('Will apply %d RGC nonlinearities with params:', numel(nonLinearitiesList))
        for iNonLinearityIdx = 1:numel(nonLinearitiesList)
            nlParamsStruct = nonLinearitiesList{iNonLinearityIdx};
            fprintf('#%d:\n', iNonLinearityIdx);
            fprintf('\ttype: %s\n', nlParamsStruct.type);
            fprintf('\tapplied to: %s\n', nlParamsStruct.sourceSignal);
            fprintf('\twithParams: \n');

            paramNames = fieldnames(nlParamsStruct.params);
            for iParam = 1:numel(paramNames)
                theParamName = paramNames{iParam};
                fprintf('\t\t %s: %f\n', theParamName, nlParamsStruct.params.(theParamName));
            end % iParam
        end
    end % Output static nonlinearity

    if (numel(theInputConeMosaicResponseTemporalSupportSeconds)>1)
        inputTimeResolutionSeconds = theInputConeMosaicResponseTemporalSupportSeconds(2)-theInputConeMosaicResponseTemporalSupportSeconds(1);
    else
        inputTimeResolutionSeconds = 1.0;
    end

    % Find out the # of trials to compute
    inputConeMosaicResponseInstancesNum = size(theInputConeMosaicResponse,1);
    if (isempty(mRGCMosaicNoisyResponseInstancesNum))
        % No 'nTrials' passed, so we will create as many instances as there
        % are cone mosaic noisy intances
        mRGCMosaicNoisyResponseInstancesNum = inputConeMosaicResponseInstancesNum;
    else
        % 'nTrials' for mRGCMosaic responses is passed
        if (coneMosaicResponseInstancesNum > 1)
            % We also have more than 1 input cone mosaic response.
            % Must be noisy cone mosaic responses instances. 
            % Ensure that mRGCMosaicNoisyResponseInstancesNum == inputConeMosaicResponseInstancesNum
            assert(mRGCMosaicNoisyResponseInstancesNum == inputConeMosaicResponseInstancesNum, ...
               '''nTrials'' (%d) does not match the trials of the input cone mosaic response (%d)', ...
               mRGCMosaicNoisyResponseInstancesNum, inputConeMosaicResponseInstancesNum);
        else
            % A single trial of input cone mosaic response, must be the noise-free cone mosaic response. Replicate it 
            theConeMosaicResponse = repmat(theInputConeMosaicResponse, [mRGCMosaicNoisyResponseInstancesNum 1 1]);
        end
    end
    nTrials = mRGCMosaicNoisyResponseInstancesNum;


    inputTimePoints = numel(theInputConeMosaicResponseTemporalSupportSeconds);
    if (isempty(timeResolutionSeconds))
        timeResolutionSeconds = inputTimeResolutionSeconds;
        responseTemporalSupportSeconds = theInputConeMosaicResponseTemporalSupportSeconds;
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

    % Compute the response of each mRGC
    parfor iRGC = 1:obj.rgcsNum
        % Retrieve the center cone indices & weights
        centerConnectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
        centerConeIndices = find(centerConnectivityVector > mRGCMosaic.minCenterWeightForInclusionInComputing);
        centerConeWeights = reshape(centerConnectivityVector(centerConeIndices), [1 1 numel(centerConeIndices)]);
        
        % Spatially pool the weighted cone responses to the RF center
        centerSpatiallyIntegratedActivations = sum(bsxfun(@times, theInputConeMosaicResponse(1:nTrials,1:inputTimePoints,centerConeIndices), centerConeWeights),3);

        if (~isempty(obj.rgcRFsurroundConeConnectivityMatrix))
            % Retrieve the surround cone indices & weights
            surroundConnectivityVector = full(squeeze(obj.rgcRFsurroundConeConnectivityMatrix(:, iRGC)));
            surroundConeIndices = find(surroundConnectivityVector > mRGCMosaic.minSurroundWeightForInclusionInComputing);
            surroundConeWeights = reshape(surroundConnectivityVector(surroundConeIndices), [1 1 numel(surroundConeIndices)]);

            % Spatially pool the weighted cone responses to the RF surround
            surroundSpatiallyIntegratedActivations = sum(bsxfun(@times, theInputConeMosaicResponse(1:nTrials,1:inputTimePoints, surroundConeIndices), surroundConeWeights),3);
        else
            surroundSpatiallyIntegratedActivations = centerSpatiallyIntegratedActivations * 0;
        end

        if (numel(theInputConeMosaicResponseTemporalSupportSeconds)>1)
            % Temporally filter the center responses
            centerSpatiallyIntegratedActivations = temporalFilter(centerSpatiallyIntegratedActivations, ...
                theInputConeMosaicResponseTemporalSupportSeconds, ...
                responseTemporalSupportSeconds, ...
                theRFcenterImpulseResponse);
    
            if (~isempty(obj.rgcRFsurroundConeConnectivityMatrix))
                % Temporally filter the surround responses
                surroundSpatiallyIntegratedActivations = temporalFilter(surroundSpatiallyIntegratedActivations, ...
                    theInputConeMosaicResponseTemporalSupportSeconds, ...
                    responseTemporalSupportSeconds, ...
                    theRFsurroundImpulseResponse);
            end
        end

        % Apply any (center/surround) component response nonlinearities
        if (numel(nonLinearitiesList)>0)
            for iNonLinearityIdx = 1:numel(nonLinearitiesList)
                nlParamsStruct = nonLinearitiesList{iNonLinearityIdx};
                if (strcmp(nlParamsStruct.sourceSignal, 'centerComponentResponse'))
                    centerSpatiallyIntegratedActivations = applyNonLinearity(centerSpatiallyIntegratedActivations, nlParamsStruct);
                end
            end % for iNonLinearityIdx

            for iNonLinearityIdx = 1:numel(nonLinearitiesList)
                nlParamsStruct = nonLinearitiesList{iNonLinearityIdx};
                if (strcmp(nlParamsStruct.sourceSignal, 'surroundComponentResponse'))
                    surroundSpatiallyIntegratedActivations = applyNonLinearity(surroundSpatiallyIntegratedActivations, nlParamsStruct);
                end
            end % for iNonLinearityIdx
        end % Input nonlinearities


        % Composite respose
        theNoiseFreeCompositeResponse = obj.responseGains(iRGC) * (centerSpatiallyIntegratedActivations - surroundSpatiallyIntegratedActivations);


        % Apply any composite response nonlinearities
        if (numel(nonLinearitiesList)>0)
            for iNonLinearityIdx = 1:numel(nonLinearitiesList)
                nlParamsStruct = nonLinearitiesList{iNonLinearityIdx};
                if (strcmp(nlParamsStruct.sourceSignal, 'compositeResponse'))
                    theNoiseFreeCompositeResponse = applyNonLinearity(theNoiseFreeCompositeResponse, nlParamsStruct);
                end
            end % for iNonLinearityIdx
        end % Output static nonlinearity

        % Store the response
        noiseFreeMRGCresponses(:,:,iRGC) = theNoiseFreeCompositeResponse;
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
    noisyMRGCresponseInstances = obj.noisyResponseInstances(noiseFreeMRGCresponses, ...
        'seed', noiseSeed);
end


function theResponse = applyNonLinearity(theResponse, nlParamsStruct)
    switch (nlParamsStruct.type)

        case 'Naka Rushton'
            fprintf('Will apply %s nonLinearity to %s with bias: %f, c50: %f, exponent: %f, super-saturation exponent: %f, and gain: %f\n', ...
                nlParamsStruct.type, ...
                nlParamsStruct.sourceSignal, ...
                nlParamsStruct.params.bias, ...
                nlParamsStruct.params.c50, ...
                nlParamsStruct.params.n, ...
                nlParamsStruct.params.s, ...
                nlParamsStruct.params.k);

                % Apply bias
                theResponse = theResponse + nlParamsStruct.params.bias;

                % Half-wave rectify
                theResponse(theResponse<0) = 0;

                % Pass through Naka Rushton activation function
                Rn = theResponse.^(nlParamsStruct.params.n);
                Rsn = theResponse.^(nlParamsStruct.params.s * nlParamsStruct.params.n);
                c50sn = nlParamsStruct.params.c50 ^ (nlParamsStruct.params.s * nlParamsStruct.params.n);
                theResponse = nlParamsStruct.params.k * Rn ./ (Rsn + c50sn);

        otherwise
            fprintf('Unknown %s non-linearity: ''%s''.', nlParamsStruct.source, nlParamsStruct.type)
    end % switch
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