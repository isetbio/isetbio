% Method to compute the spatiotemporal response of the mRGCMosaic given the response of its input cone mosaic
function [noiseFreeMRGCresponses, noisyMRGCresponseInstances, responseTemporalSupportSeconds, noiseFreeLinearMRGCresponses] = compute(obj, ...
            theInputConeMosaicResponse, theInputConeMosaicResponseTemporalSupportSeconds, varargin)

    p = inputParser;
    p.addParameter('withCenterSurroundTemporalFilters', [], @(x)(isempty(x)||isstruct(x)));
    p.addParameter('nTrials', [], @isscalar);
    p.addParameter('nonLinearitiesList', [], @(x)(isempty(x))||(iscell(x))||(isstruct(x)));
    p.addParameter('computeResponsesOnlyForCellsWithIndices', [], @isnumeric);
    p.addParameter('deactivatedCenter', false, @islogical);
    p.addParameter('deactivatedSurround', false, @islogical);
    p.addParameter('flipLinearResponsePolarityForCellsWithIndices', [], @isnumeric);
    p.addParameter('seed', [], @isnumeric);
    p.addParameter('beVerbose', false, @islogical);

    % Parse input
    p.parse(varargin{:});
    centerSurroundTemporalFilters = p.Results.withCenterSurroundTemporalFilters;
    mRGCMosaicNoisyResponseInstancesNum = p.Results.nTrials;
    nonLinearitiesList = p.Results.nonLinearitiesList;
    deactivatedCenter = p.Results.deactivatedCenter;
    deactivatedSurround = p.Results.deactivatedSurround;
    flipLinearResponsePolarityForCellsWithIndices = p.Results.flipLinearResponsePolarityForCellsWithIndices;
    computeResponsesOnlyForCellsWithIndices = p.Results.computeResponsesOnlyForCellsWithIndices;
    noiseSeed = p.Results.seed;
    beVerbose = p.Results.beVerbose;

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
        nonLinearitiesList = cell(1,1);
        nonLinearitiesList{1} = theSingleNonLinearity;
    end

    % Nonlinearities
    if (numel(nonLinearitiesList)>0)
        if (beVerbose)
            fprintf('Will apply %d RGC nonlinearities with params:', numel(nonLinearitiesList));
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
        end % beVerbose
    end % if (numel(nonLinearitiesList)>0)

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



    inputTimePoints = numel(theInputConeMosaicResponseTemporalSupportSeconds);

    % If we have centerSurroundTemporalFilters, examine their time base and
    % interpolate theInputConeMosaicResponse to match that of the filters

    theInterpolatedInputConeMosaicTemporalSupportSeconds = [];
    if (~isempty(centerSurroundTemporalFilters))
        
        % Validate inner retina temporal filter data
        assert(isstruct(centerSurroundTemporalFilters), ...
            '''centerSurroundTemporalFilters'' must be a struct');

        assert(isfield(centerSurroundTemporalFilters, 'temporalSupportSeconds'), ...
            '''centerSurroundTemporalFilters'' must have a ''temporalSupportSeconds'' field');

        assert(isfield(centerSurroundTemporalFilters, 'centerImpulseResponseFunction'), ...
            '''centerSurroundTemporalFilters'' must have a ''centerImpulseResponseFunction'' field');

        assert(isfield(centerSurroundTemporalFilters, 'surroundImpulseResponseFunction'), ...
            '''centerSurroundTemporalFilters'' must have a ''surroundImpulseResponseFunction'' field');

        assert(numel(centerSurroundTemporalFilters.centerImpulseResponseFunction) == numel(centerSurroundTemporalFilters.temporalSupportSeconds), ...
            '''centerSurroundTemporalFilters'' centerImpulseResponseFunction must have the same length as its temporal support');

        assert(numel(centerSurroundTemporalFilters.surroundImpulseResponseFunction) == numel(centerSurroundTemporalFilters.temporalSupportSeconds), ...
            '''centerSurroundTemporalFilters'' surroundImpulseResponseFunction must have the same length as its temporal support');
        
        % Ensure center and surround temporal filters have unit volume
        theVolume = sum(abs(centerSurroundTemporalFilters.centerImpulseResponseFunction(:)));
        centerSurroundTemporalFilters.centerImpulseResponseFunction = centerSurroundTemporalFilters.centerImpulseResponseFunction / theVolume;
        theVolume = sum(abs(centerSurroundTemporalFilters.surroundImpulseResponseFunction(:)));
        centerSurroundTemporalFilters.surroundImpulseResponseFunction = centerSurroundTemporalFilters.surroundImpulseResponseFunction / theVolume;

        % New timebase for interpolated input responses
        if (numel(theInputConeMosaicResponseTemporalSupportSeconds)>1)
            dTinnerRetina = centerSurroundTemporalFilters.temporalSupportSeconds(2)-centerSurroundTemporalFilters.temporalSupportSeconds(1);
            dtInput = theInputConeMosaicResponseTemporalSupportSeconds(2)-theInputConeMosaicResponseTemporalSupportSeconds(1);
            if (abs(dTinnerRetina-dtInput)*1e3 > 10*eps)
                fprintf('Will interpolate inputs to new time base with dT = %2.1f msec (from dT = %2.1f msec)\n', dTinnerRetina*1e3, dtInput*1e3);
                theInterpolatedInputConeMosaicTemporalSupportSeconds = ...
                    theInputConeMosaicResponseTemporalSupportSeconds(1):dTinnerRetina:theInputConeMosaicResponseTemporalSupportSeconds(end);
            end
        end

    end


    % Compute the output response temporal support
    responseTemporalSupportSeconds = theInputConeMosaicResponseTemporalSupportSeconds;

    dtRaw = [];
    dtInterpolated = [];

    if (~isempty(centerSurroundTemporalFilters))
        if (~isempty(theInterpolatedInputConeMosaicTemporalSupportSeconds))
           temporallyFilteredActivationsLength = numel(theInterpolatedInputConeMosaicTemporalSupportSeconds) + numel(centerSurroundTemporalFilters.temporalSupportSeconds)-1;
           dtRaw = theInputConeMosaicResponseTemporalSupportSeconds(2)-theInputConeMosaicResponseTemporalSupportSeconds(1);
           dtInterpolated = theInterpolatedInputConeMosaicTemporalSupportSeconds(2)-theInterpolatedInputConeMosaicTemporalSupportSeconds(1);
           responseTemporalSupportSeconds = (0:1:(temporallyFilteredActivationsLength-1)) * dtInterpolated;
        end
    end
    

    % Allocate memory for the computed responses
    noiseFreeMRGCresponses = zeros(mRGCMosaicNoisyResponseInstancesNum, numel(responseTemporalSupportSeconds), obj.rgcsNum);
    noiseFreeLinearMRGCresponses = zeros(mRGCMosaicNoisyResponseInstancesNum, numel(responseTemporalSupportSeconds), obj.rgcsNum);

    % Compute the response of each mRGC
    parfor iRGC = 1:obj.rgcsNum

        if (~isempty(computeResponsesOnlyForCellsWithIndices))
            if (~ismember(iRGC, computeResponsesOnlyForCellsWithIndices))
                fprintf('Skipping computation of mRGC #%d response.\n', iRGC);
                continue;
            end
        end

        % Retrieve the center cone indices & weights
        centerConnectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
        centerConeIndices = find(centerConnectivityVector > mRGCMosaic.minCenterWeightForInclusionInComputing);
        centerConeWeights = reshape(centerConnectivityVector(centerConeIndices), [1 1 numel(centerConeIndices)]);
        
        if (deactivatedCenter)
            centerConeWeights = 0*centerConeWeights;
        end

        % Spatially pool the weighted cone responses to the RF center
        centerSpatiallyIntegratedActivations = sum(bsxfun(@times, theInputConeMosaicResponse(:,1:inputTimePoints,centerConeIndices), centerConeWeights),3);

        if (~isempty(obj.rgcRFsurroundConeConnectivityMatrix))
            % Retrieve the surround cone indices & weights
            surroundConnectivityVector = full(squeeze(obj.rgcRFsurroundConeConnectivityMatrix(:, iRGC)));
            surroundConeIndices = find(surroundConnectivityVector > mRGCMosaic.minSurroundWeightForInclusionInComputing);
            surroundConeWeights = reshape(surroundConnectivityVector(surroundConeIndices), [1 1 numel(surroundConeIndices)]);

            if (deactivatedSurround)
                surroundConeWeights = 0*surroundConeWeights;
            end

            % Spatially pool the weighted cone responses to the RF surround
            surroundSpatiallyIntegratedActivations = sum(bsxfun(@times, theInputConeMosaicResponse(:,1:inputTimePoints, surroundConeIndices), surroundConeWeights),3);
        else
            surroundSpatiallyIntegratedActivations = centerSpatiallyIntegratedActivations * 0;
        end


        % Apply inner retina tempora filters to the spatially integrated
        % center/surround responses
        if (~isempty(centerSurroundTemporalFilters))

            % Interpolate responses
            if (~isempty(theInterpolatedInputConeMosaicTemporalSupportSeconds))
                
                nTrials = size(centerSpatiallyIntegratedActivations,1);
                
                temporallyFilteredActivationsLength = numel(theInterpolatedInputConeMosaicTemporalSupportSeconds) + numel(centerSurroundTemporalFilters.temporalSupportSeconds)-1;
                centerSpatiallyIntegratedTemporallyFilteredActivations = zeros(nTrials, temporallyFilteredActivationsLength);
                surroundSpatiallyIntegratedTemporallyFilteredActivations = zeros(nTrials, temporallyFilteredActivationsLength);

                for iTrial = 1:nTrials

                    centerSpatiallyIntegratedTemporallyFilteredActivations(iTrial,:) = applyTemporalFilter(...
                        theInputConeMosaicResponseTemporalSupportSeconds, ...
                        squeeze(centerSpatiallyIntegratedActivations(iTrial,:)), ...
                        theInterpolatedInputConeMosaicTemporalSupportSeconds, ...
                        centerSurroundTemporalFilters.centerImpulseResponseFunction, ...
                        dtInterpolated/dtRaw);

                    surroundSpatiallyIntegratedTemporallyFilteredActivations(iTrial,:) = applyTemporalFilter(...
                        theInputConeMosaicResponseTemporalSupportSeconds, ...
                        squeeze(surroundSpatiallyIntegratedActivations(iTrial,:)), ...
                        theInterpolatedInputConeMosaicTemporalSupportSeconds, ...
                        centerSurroundTemporalFilters.surroundImpulseResponseFunction, ...
                        dtInterpolated/dtRaw);
                end

                % Replace unfiltered with filtered responses 
                centerSpatiallyIntegratedActivations = centerSpatiallyIntegratedTemporallyFilteredActivations;
                surroundSpatiallyIntegratedActivations =  surroundSpatiallyIntegratedTemporallyFilteredActivations;
                
            end
        end




        % Composite response before any non-linearities are applied
        theNoiseFreeCompositeLinearResponse = obj.responseGains(iRGC) * (centerSpatiallyIntegratedActivations - surroundSpatiallyIntegratedActivations);

        % Flip composite response sign if so specified (simulating OFF mosaic)
        if (~isempty(find(flipLinearResponsePolarityForCellsWithIndices == iRGC)))
            theNoiseFreeCompositeLinearResponse = -theNoiseFreeCompositeLinearResponse;
        end

        noiseFreeLinearMRGCresponses(:,:,iRGC) = theNoiseFreeCompositeLinearResponse;


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


        % Flip composite response sign if so specified (simulating OFF mosaic)
        if (~isempty(find(flipLinearResponsePolarityForCellsWithIndices == iRGC)))
            theNoiseFreeCompositeResponse = -theNoiseFreeCompositeResponse;
        end


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


function theNonLinearResponse = applyNonLinearity(theLinearResponse, nlParamsStruct)
    switch (nlParamsStruct.type)

        case 'Naka Rushton'
            fprintf('Will apply %s nonLinearity to %s with RinBaseline: %f, c50: %f, exponent: %f, and super-saturation exponent: %f\n', ...
                    nlParamsStruct.type, ...
                    nlParamsStruct.sourceSignal, ...
                    nlParamsStruct.params.RinHalfWaveRectifierThreshold, ...
                    nlParamsStruct.params.Rin50, ...
                    nlParamsStruct.params.n, ...
                    nlParamsStruct.params.s);
            
            % Half-wave rectification
            theRectifiedResponse = theLinearResponse - nlParamsStruct.params.RinHalfWaveRectifierThreshold;
            theRectifiedResponse(theRectifiedResponse < 0) = 0;

            % Exponent cannot be less than 0
            nlParamsStruct.params.n  = max([0 nlParamsStruct.params.n]);

            % Super-saturation cannot be less than 1
            nlParamsStruct.params.s = max([1 nlParamsStruct.params.s]);

            % Rin50 cannot be less than 0
            nlParamsStruct.params.Rin50 = max([0 nlParamsStruct.params.Rin50]);

            % Pass through Naka Rushton activation function
            Rn = theRectifiedResponse .^ (nlParamsStruct.params.n);

            sn = nlParamsStruct.params.s * nlParamsStruct.params.n;
            Rsn = theRectifiedResponse.^ sn;
            R50sn = nlParamsStruct.params.Rin50 ^ sn;

            theNonLinearResponse = Rn ./ (Rsn + R50sn);

            % Determine gain to bring the nlParamsStruct.params.RoutGainAdjustmentPercentile of the
            % theRectifiedResponse in same range as theNonLinearResponse
            p = prctile(theRectifiedResponse, nlParamsStruct.params.RoutGainAdjustmentPercentile);
            idx = find(theRectifiedResponse<=p);
            gain = theRectifiedResponse(idx) / theNonLinearResponse(idx);

            theNonLinearResponse = theNonLinearResponse * gain;
        otherwise
            fprintf('Unknown %s non-linearity: ''%s''.', nlParamsStruct.source, nlParamsStruct.type)
    end % switch
end


function theFilteredResponse = applyTemporalFilter(theResponseTemporalSupportSeconds, theResponse, ...
              theInterpolatedResponseTemporalSupportSeconds, theFilterImpulseResponseFunction, ...
              theScalingFactor)

    % Interpolate to desired timebase, but keep same volume
    theResponse = interp1(theResponseTemporalSupportSeconds, theResponse, theInterpolatedResponseTemporalSupportSeconds);
                    
    % Filter center integrated response
    theFilteredResponse = theScalingFactor * conv(theResponse, theFilterImpulseResponseFunction);
end
