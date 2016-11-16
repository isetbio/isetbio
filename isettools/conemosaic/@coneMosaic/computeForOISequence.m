function [absorptions, absorptionsTimeAxis, varargout] = computeForOISequence(obj, oiSequence, varargin)
% Compute cone absorptions and (optionally) cone photocurrent for a
% sequence of optical images (@oiSequence). It is also possible to run this
% for a multiple eye movement paths,
%
%    [absorptions, absorpionsTimeAxis, [current, currentTimeAxis]] = ...
%         cMosaic.compute(oiSequence, varargin);
%
% This method is typically called by coneMosaic.compute(oiSequence)
%
% Inputs:
%   oiSequence  - an @oiSequence object
%
% Optional inputs:
%   emPaths      - N eye movement paths, each with M eye movements
%                  organized in an [N x M x 2] matrix
%   currentFlag  - logical, whether to compute photocurrent
%   newNoise     - logical, whether to use new random seed
%
% Outputs:
%   absorptions          - cone photon absorptions (photon counts in integrationTime)
%   absorptionsTimeAxis  - time axis for the absorptions signal
%   photocurrent         - cone photocurrent
%   photocurrentTimeAxis - time axis for photocurrent signal
%
% Examples:
%
% See also: coneMosaic.compute
%
% NPC ISETBIO Team 2016

%% Parse inputs
p = inputParser;
p.addRequired('oiSequence', @(x)isa(x, 'oiSequence'));
p.addParameter('emPaths', [], @isnumeric);
p.addParameter('currentFlag', false, @islogical);
p.addParameter('newNoise', true, @islogical);
p.parse(oiSequence, varargin{:});

oiSequence = p.Results.oiSequence;
emPaths = p.Results.emPaths;
currentFlag = p.Results.currentFlag;
newNoise = p.Results.newNoise;
oiTimeAxis = oiSequence.oiTimeAxis;

if (oiSequence.length ~= numel(oiTimeAxis))
    error('oiTimeAxis and oiSequence must have equal length\n');
end

if (isempty(emPaths))
    emPaths = reshape(obj.emPositions, [1 size(obj.emPositions)]);
end

if (isempty(emPaths))
    error('Either supply an ''emPaths'' key-value pair, or preload coneMosaic.emPositions');
end

%% Compute eye movement time axis
instancesNum = size(emPaths,1);
eyeMovementsNum = size(emPaths,2);
eyeMovementTimeAxis = oiTimeAxis(1) + (0:1:(eyeMovementsNum-1)) * obj.integrationTime;

%% Compute OIrefresh
oiRefreshInterval = oiTimeAxis(2)-oiTimeAxis(1);

% Save default integration time
defaultIntegrationTime = obj.integrationTime;

% Only allocate memory for the non-null cones in a 3D matrix [instances x
% time x numel(nonNullConesIndices)]
nonNullConesIndices = find(obj.pattern>1);
absorptions = zeros(instancesNum, numel(eyeMovementTimeAxis), numel(nonNullConesIndices), 'single');

if (oiRefreshInterval >= defaultIntegrationTime)
    % There are two main time sampling scenarios.  This one is when the oi
    % update rate is SLOWER than the cone integration time which is also
    % equal to the eye movement update rate
    
    % SCENARIO
    %              |----- oiRefreshInterval ----|----- oiRefreshInterval ----|----- oiRefreshInterval ----|
    %      oi      |                            |                            |                            |
    %   sequence   |                            |<-----tFrameStart           |                            |
    % _____________.____________________________|           tFrameEnd------->|____________________________|
    %     eye      |          |          |          |          |          |          |          |          |
    %   movement   |          |          |-intTime->|-intTime->|-intTime->|          |          |          |
    %   sequence   |          |          |          |          |          |          |          |          |
    % ------------------------------------******xxxx|**********|**********|--------------------------------|
    %                                       p1   p2     full        full
    %                   partial absorption_/      \_partial absorption  \_ full absorption
    %
    
    % Loop over the optical images
    for oiIndex = 1:oiSequence.length
        
        % Current oi time limits
        tFrameStart = oiTimeAxis(oiIndex);
        tFrameEnd   = tFrameStart + oiRefreshInterval;
        
        % Find eye movement indices withing the oi limits
        indices = find( (eyeMovementTimeAxis >  tFrameStart-defaultIntegrationTime) & ...
            (eyeMovementTimeAxis <= tFrameEnd-defaultIntegrationTime+eps) );
        
        if (isempty(indices))
            %  time samples in
            % the mosaic than we have oi samples.  That should be OK.
            disp('Fewer Eye movement time samples than oi samples')
            break;
            %error('empty indices');
        end
        % the first eye movement requires special treatment as it may have started before the current frame,
        % so we need to compute partial absorptions over the previous frame and over the current frame
        idx = indices(1);
        integrationTimeForFirstPartialAbsorption = tFrameStart-eyeMovementTimeAxis(idx);
        integrationTimeForSecondPartialAbsorption = eyeMovementTimeAxis(idx)+defaultIntegrationTime-tFrameStart;
        
        % Partial absorptions (p1 in graph above) with previous oi
        % (across all instances)
        if (oiIndex > 1)
            % Update the @coneMosaic with the partial integration time
            obj.integrationTime = integrationTimeForFirstPartialAbsorption;
            % Compute partial absorptions
            emSubPath = reshape(squeeze(emPaths(1:instancesNum,idx,:)), [instancesNum 2]);
            absorptionsDuringPreviousFrame = obj.compute(...
                oiSequence.frameAtIndex(oiIndex-1), ...
                'emPath', emSubPath, ...
                'newNoise', newNoise, ...
                'currentFlag', false ...
                );
        else
            absorptionsDuringPreviousFrame = zeros(size(obj.pattern,1), size(obj.pattern,2), instancesNum);
        end
        
        % Partial absorptions (p2 in graph above) with current oi
        % across all instances
        % Update the @coneMosaic with the partial integration time
        obj.integrationTime = integrationTimeForSecondPartialAbsorption;
        % Compute partial absorptions
        emSubPath = reshape(squeeze(emPaths(1:instancesNum,idx,:)), [instancesNum 2]);
        absorptionsDuringCurrentFrame =  obj.compute(...
            oiSequence.frameAtIndex(oiIndex), ...
            'emPath', emSubPath, ...
            'newNoise', newNoise, ...
            'currentFlag', false ...
            );
        
        % summed absorptions
        totalAbsorptions = absorptionsDuringPreviousFrame+absorptionsDuringCurrentFrame;
        
        % Only get the absorptions for the non-null cones
        totalAbsorptions = reshape(permute(totalAbsorptions, [3 1 2]), [instancesNum size(obj.pattern,1) * size(obj.pattern,2)]);
        totalAbsorptions = reshape(totalAbsorptions(:, nonNullConesIndices), [instancesNum 1 numel(nonNullConesIndices)]);
        
        % insert the sum of the two partial absorptions in the time series
        insertionIndex = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
        absorptions(1:instancesNum, insertionIndex, :) = single(totalAbsorptions);
        
        
        % Full absorptions with current oi and default integration time)
        if (numel(indices)>1)
            % Update the @coneMosaic with the default integration time
            obj.integrationTime = defaultIntegrationTime;
            % Compute absorptions for all remaining the OIs
            idx = indices(2:end);
            emSubPath = reshape(emPaths(1:instancesNum, idx,:), [instancesNum*numel(idx) 2]);
            absorptionsForRemainingEyeMovements = obj.compute(...
                oiSequence.frameAtIndex(oiIndex), ...
                'emPath', emSubPath, ...
                'newNoise', newNoise, ...
                'currentFlag', false ...
                );
            
            % Only get the absorptions for the non-null cones
            absorptionsForRemainingEyeMovements = reshape(permute(absorptionsForRemainingEyeMovements, [3 1 2]), [instancesNum*numel(idx) size(obj.pattern,1)*size(obj.pattern,2)]);
            absorptionsForRemainingEyeMovements = reshape(absorptionsForRemainingEyeMovements(:, nonNullConesIndices), [instancesNum numel(idx) numel(nonNullConesIndices)]);
            
            % insert in time series
            insertionIndices = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
            absorptions(1:instancesNum, insertionIndices,:,:) = single(absorptionsForRemainingEyeMovements);
            
        end
    end  % oiIndex
    
    % oiRefreshInterval > defaultIntegrationTime
else
    % There are two main time sampling scenarios.  This one is when the oi
    % update rate is FASTER than the cone integration time which is also equal
    % to the eye movement update rate
    
    % SCENARIO
    %     eye      |                             |                             |                             |
    %   movement   |                             |                             |                             |
    %   sequence   |                             |<-----emStart                |                             |
    % _____________._____________________________|                emEnd------->|_____________________________|
    %              |            |            |            |            |            |          |          |
    %      oi      |            |            |-oiRefresh->|-oiRefresh->|-oiRefresh->|          |          |
    %   sequence   |            |            |    ********|************|*******     |          |          |
    % ---------------------------------------------partial-----full-----partial-------------------------------
    %                                            absorption  absorption absorption
    
    % Loop over the eye movements
    for emIndex = 1:eyeMovementsNum
        
        % Current eye movement time limits
        emStart = eyeMovementTimeAxis(emIndex);
        emEnd   = emStart + defaultIntegrationTime;
        
        % sum the partial integration times
        actualIntegrationTime = 0;
        
        % Find oi indices withing the eye movement frame time limits
        indices = find( (oiTimeAxis > emStart-oiRefreshInterval) & ...
            (oiTimeAxis <= emEnd + eps) );
        
        if (isempty(indices))
            disp('Fewer eye movement time samples than oi samples')
            break;
            % error('empty indices');
        end
        
        % Partial absorptions during the ovelap with the OI that started before the emStart
        idx = indices(1);
        integrationTimeForFirstPartialAbsorption = oiTimeAxis(idx)+oiRefreshInterval-emStart;
        % Update the @coneMosaic with the partial integration time
        obj.integrationTime = integrationTimeForFirstPartialAbsorption;
        % Update the sum of partial integration times
        actualIntegrationTime = actualIntegrationTime + obj.integrationTime;
        % Compute absorptions
        emSubPath = reshape(emPaths(1:instancesNum, emIndex,:), [instancesNum 2]);
        absorptionsAccum = obj.compute(...
            oiSequence.frameAtIndex(idx), ...
            'emPath', emSubPath, ...
            'newNoise', newNoise, ...
            'currentFlag', false ...
            );
        
        % Next, compute full absorptions for the remaining OIs
        if (numel(indices)>1)
            for k = 2:numel(indices)
                % Update the @coneMosaic with the partial integration time
                if (k < numel(indices))
                    % full absorption during the duration of the OI
                    obj.integrationTime = oiRefreshInterval;
                else
                    % partial absorption during the last OI
                    idx = indices(end);
                    integrationTimeForLastPartialAbsorption = emEnd - oiTimeAxis(idx);
                    obj.integrationTime = integrationTimeForLastPartialAbsorption;
                end
                % Update the sum of partial integration times
                actualIntegrationTime = actualIntegrationTime + obj.integrationTime;
                % Compute absorptions
                emSubPath = reshape(emPaths(1:instancesNum, emIndex,:), [instancesNum 2]);
                absorptionsAccum = absorptionsAccum + obj.compute(...
                    oiSequence.frameAtIndex(indices(k)), ...
                    'emPath', emSubPath, ...
                    'newNoise', newNoise, ...
                    'currentFlag', false ...
                    );
            end  % for k
        end
        
        % Check the sum of partial absorptions equal the default absorption time
        if (abs(actualIntegrationTime-defaultIntegrationTime) > defaultIntegrationTime*0.0001) && (emIndex < numel(eyeMovementTimeAxis))
            error('Actual integration time (%3.5f) not equal to desired value (%3.5f) [emIndex: %d / %d]\n', actualIntegrationTime, defaultIntegrationTime, emIndex, numel(eyeMovementTimeAxis));
        end
        
        % Only get the absorptions for the non-null cones
        absorptionsAccum = reshape(permute(absorptionsAccum, [3 1 2]), [instancesNum size(obj.pattern,1) * size(obj.pattern,2)]);
        absorptionsAccum = reshape(absorptionsAccum(:, nonNullConesIndices), [instancesNum 1 numel(nonNullConesIndices)]);
        
        % insert to time series
        insertionIndices = round((eyeMovementTimeAxis(emIndex)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
        absorptions(1:instancesNum, insertionIndices, :) = single(absorptionsAccum);
    end % emIndex
end % oiRefreshInterval > defaultIntegrationTime

% Restore default integrationTime
obj.integrationTime = defaultIntegrationTime;

% Reload the full eye movement sequence for the first instance only
if (ndims(emPaths) == 3)
    obj.emPositions = squeeze(emPaths(1,:,:));
else
    obj.emPositions = emPaths;
end

if (isa(obj, 'coneMosaicHex'))
    % Reshape absorptions to correct dimensions [instances, cone_indices, time]
    absorptions = permute(absorptions, [1 3 2]);
    
    if (currentFlag)
        % Add one more absorption at the end
        absorptions = cat(3, absorptions, squeeze(absorptions(:,:,end)));
    end
else
    % Reshape absorptions to correct dimensions [instances, cone_rows, cone_cols, time]
    absorptions = permute(reshape(absorptions, [instancesNum numel(eyeMovementTimeAxis) size(obj.pattern,1) size(obj.pattern,2)]),[1 3 4 2]);
    
    if (currentFlag)
        % Add one more absorption at the end
        absorptions = cat(4, absorptions, absorptions(:,:,:,end));
    end
end

% Reload the absorptions signal from the last instance
% @BW Replaced the original code: squeeze(absorptions(end,:,:,:));
% The squeeze created problems for 1D cone mosaic data.  This method keeps
% the absoprtions a 3D matrix (row,col,time).
obj.absorptions = reshape(absorptions,obj.rows,obj.cols,size(absorptions,4));

% align absorptions time axis with respect to optical image sequence time axis
absorptionsTimeAxis = oiTimeAxis(1) + obj.absorptionsTimeAxis;

if (currentFlag)
    % compute the photocurrent time axis It is possible to call
    % coneMosaic.computeCurrent at a later time, rather than set this flag.
    % @BW:  I think it might be better to replace this code with that call
    % so that we only have one computeCurrent method.
    dtOS = obj.os.timeStep;
    osTimeAxis = absorptionsTimeAxis(1): dtOS : absorptionsTimeAxis(end);
    
    if (isa(obj, 'coneMosaicHex'))
        photocurrents = zeros(instancesNum, numel(nonNullConesIndices), numel(osTimeAxis), 'single');
        for instanceIndex = 1:instancesNum
            tmp = squeeze(absorptions(instanceIndex,:,:));
            
            % Resample to osTimeAxis
            tmp = coneMosaic.tResample(tmp, obj.pattern(nonNullConesIndices), absorptionsTimeAxis, osTimeAxis);
            
            % osCompute expects a 3D pRate, so make it so
            tmp = reshape(tmp, [size(tmp,1) 1 size(tmp,2)]);
            % Compute photocurrent from photonRate (tmp/dtOS)
            tmp = single(obj.os.osCompute(tmp/dtOS, obj.pattern(nonNullConesIndices), 'append', false));
            % Put it back in correct shape
            photocurrents(instanceIndex,:,:) = permute(tmp, [2 1 3]);
        end
    else
        photocurrents = zeros(instancesNum, size(obj.pattern,1), size(obj.pattern,2), numel(osTimeAxis), 'single');
        
        for instanceIndex = 1:instancesNum
            tmp = squeeze(absorptions(instanceIndex,:,:,:));
            
            % Resample to osTimeAxis (reshape to 2D for faster processing)
            % Reshape needed for spatial case where we have singleton dimensions
            tmp = reshape(tmp, [size(obj.pattern,1) size(obj.pattern,2) numel(absorptionsTimeAxis)]);
            tmp = coneMosaic.tResample(tmp, obj.pattern, absorptionsTimeAxis, osTimeAxis);
            tmp = reshape(tmp, [size(obj.pattern,1) size(obj.pattern,2) numel(osTimeAxis)]);
            
            % Compute photocurrent from photonRate (tmp/dtOS)
            photocurrents(instanceIndex,:,:,:) = single(obj.os.osCompute(tmp/dtOS, obj.pattern, 'append', false));
        end
    end
    
    % Remove the last absorption we inserted at the end
    if (isa(obj, 'coneMosaicHex'))
        absorptions = absorptions(:,:,1:(end-1));
    else
        absorptions = absorptions(:,:,:,1:(end-1));
    end
    
    % Reload the absorptions signal from the last instance
    sz = size(absorptions);
    lastInstance = sz(1);
    obj.absorptions = reshape(squeeze(absorptions(lastInstance,:,:,:)), [size(obj.pattern,1) size(obj.pattern,2) sz(end)]);
    
    % Re-align absorptions time axis with respect to optical image sequence time axis
    absorptionsTimeAxis = oiTimeAxis(1) + obj.absorptionsTimeAxis;
    
    % Return photocurrents
    varargout{1} = photocurrents;
    
    % Return the photocurrent time axis
    varargout{2} = osTimeAxis;
else
    varargout{1} = [];
    varargout{2} = [];
end % currentFlag

end