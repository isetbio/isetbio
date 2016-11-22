function [absorptions, photocurrents] = computeForOISequence(obj, oiSequence, varargin)
% Compute cone absorptions and (optionally) cone photocurrents for a
% sequence of optical images (@oiSequence). 
%
%    [absorptions, photocurrents] = cMosaic.compute(oiSequence, varargin);
%
% TO CLARIFY:
%
% It is also possible to run this for a multiple eye movement paths.  I am
% not entirely clear about how we are managing the number of trials and
% sequences here.  I think the number of trials is determined by the number
% of eye movement sequences.  The oiSequence itself is always just one
% current time series. (BW).  
%
% Also, NC is adding features into this and trapping cases. We need to
% integrate.
%
% This method is typically called by coneMosaic.compute(oiSequence)
%
% The returned absorptions has an extra dimension (the first one) so that
% we can calculate for multiple eye movement paths.  The absorptions from a
% single eye movement case would be
%
%    absorptions(thisTrial,:,:,:)
%
% The coneMosaic object (obj) always has the absorptions from the last
% trial and dimensions (row,col,time).
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
%   photocurrent         - cone photocurrent
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
p.addParameter('workerID', [], @isnumeric);
p.addParameter('workDescription', '', @ischar);
p.parse(oiSequence, varargin{:});

oiSequence  = p.Results.oiSequence;
emPaths     = p.Results.emPaths;
currentFlag = p.Results.currentFlag;
newNoise    = p.Results.newNoise;
workerID    = p.Results.workerID;
workDescription = p.Results.workDescription;
oiTimeAxis  = oiSequence.oiTimeAxis;

nTimes = numel(oiTimeAxis);
if (oiSequence.length ~= nTimes)
    error('oiTimeAxis and oiSequence must have equal length\n');
end

if (isempty(emPaths))
    emPaths = reshape(obj.emPositions, [1 size(obj.emPositions)]);
end

if (isempty(emPaths))
    error('Either supply an ''emPaths'' key-value pair, or preload coneMosaic.emPositions');
end

varargout{1} = [];
varargout{2} = [];

%% Save default integration time
defaultIntegrationTime = obj.integrationTime;

%% Compute eye movement time axis
nTrials    = size(emPaths,1);
eyeMovementsNum = size(emPaths,2);
eyeMovementTimeAxis = oiTimeAxis(1) + (0:1:(eyeMovementsNum-1)) * obj.integrationTime;

%% Compute OIrefresh
if (numel(oiTimeAxis) == 1)
    oiRefreshInterval = defaultIntegrationTime;
else
    oiRefreshInterval = oiTimeAxis(2)-oiTimeAxis(1);
end

% Only allocate memory for the non-null cones in a 3D matrix [instances x numel(nonNullConesIndices) x time]
nonNullConesIndices = find(obj.pattern>1);
absorptions = zeros(nTrials, numel(nonNullConesIndices), numel(eyeMovementTimeAxis), 'single');

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
        
        if (~isempty(workerID))
            % Update progress in command window
            displayProgress(workerID, workDescription, 0.5*oiIndex/oiSequence.length);
        end
            
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
            emSubPath = reshape(squeeze(emPaths(1:nTrials,idx,:)), [nTrials 2]);
            obj.absorptions = [];
            absorptionsDuringPreviousFrame = obj.compute(...
                oiSequence.frameAtIndex(oiIndex-1), ...
                'emPath', emSubPath, ...
                'newNoise', newNoise, ...
                'currentFlag', false ...
                );
        else
            absorptionsDuringPreviousFrame = zeros(size(obj.pattern,1), size(obj.pattern,2), nTrials);
        end
        
        % Partial absorptions (p2 in graph above) with current oi
        % across all instances
        % Update the @coneMosaic with the partial integration time
        obj.integrationTime = integrationTimeForSecondPartialAbsorption;
        % Compute partial absorptions
        emSubPath = reshape(squeeze(emPaths(1:nTrials,idx,:)), [nTrials 2]);
        absorptionsDuringCurrentFrame =  obj.compute(...
            oiSequence.frameAtIndex(oiIndex), ...
            'emPath', emSubPath, ...
            'newNoise', newNoise, ...
            'currentFlag', false ...
            );
        
        % summed absorptions
        absorptionsAllTrials = absorptionsDuringPreviousFrame+absorptionsDuringCurrentFrame;
        
        % Reformat and insert to time series  
        insertionIndices = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1; 
        reformatAbsorptionsAllTrialsMatrix(nTrials, numel(insertionIndices), size(obj.pattern,1), size(obj.pattern,2));
        absorptions(1:nTrials, :, insertionIndices) = absorptionsAllTrials;
      
        % Full absorptions with current oi and default integration time)
        if (numel(indices)>1)
            % Update the @coneMosaic with the default integration time
            obj.integrationTime = defaultIntegrationTime;
            % Compute absorptions for all remaining the OIs
            idx = indices(2:end);
            emSubPath = reshape(emPaths(1:nTrials, idx,:), [nTrials*numel(idx) 2]);
            absorptionsAllTrials = obj.compute(...
                oiSequence.frameAtIndex(oiIndex), ...
                'emPath', emSubPath, ...
                'newNoise', newNoise, ...
                'currentFlag', false ...
                );
            % Reformat and insert to time series     
            insertionIndices = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
            reformatAbsorptionsAllTrialsMatrix(nTrials, numel(insertionIndices), size(obj.pattern,1), size(obj.pattern,2));
            absorptions(1:nTrials, :, insertionIndices) = absorptionsAllTrials;
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
        
        if (~isempty(workerID))
            displayProgress(workerID, workDescription, 0.5*emIndex/eyeMovementsNum);
        end
            
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
        emSubPath = reshape(emPaths(1:nTrials, emIndex,:), [nTrials 2]);
        obj.absorptions = [];
        absorptionsAllTrials = obj.compute(...
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
                emSubPath = reshape(emPaths(1:nTrials, emIndex,:), [nTrials 2]);
                obj.absorptions = [];
                absorptionsAllTrials = absorptionsAllTrials + obj.compute(...
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
        
        % Reformat and insert to time series  
        insertionIndices = round((eyeMovementTimeAxis(emIndex)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
        reformatAbsorptionsAllTrialsMatrix(nTrials, numel(insertionIndices), size(obj.pattern,1), size(obj.pattern,2));
        absorptions(1:nTrials, :, insertionIndices) = absorptionsAllTrials;
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


if (~currentFlag) || (numel(eyeMovementTimeAxis) == 1)
   % If we don't compute the current, then we just return the absorptions from the last trial.
   if (isa(obj, 'coneMosaicHex')) 
       % Reshape to full hex pattern
       obj.absorptions = oneDhexTo2Dhex(squeeze(absorptions(nTrials,:,:)), obj.pattern);
   else
        % Reshape absorptions to correct dimensions [instances, cone_rows, cone_cols, time]
        absorptions = reshape(absorptions, [nTrials size(obj.pattern,1) size(obj.pattern,2) numel(eyeMovementTimeAxis)]);
        
        % If we don't compute the current, then we just return the absorptions from the last trial.
        obj.absorptions = reshape(squeeze(absorptions(nTrials,:,:)), [size(obj.pattern,1) size(obj.pattern,2) numel(eyeMovementTimeAxis)]);
   end
   % Return to caller
   return;
end


% Align absorptions time axis with respect to optical image sequence time
% axis; Note sure what this is.  Ask NC.
% absorptionsTimeAxis = oiTimeAxis(1) +  obj.timeAxis;

    
photocurrents = [];
if (currentFlag)
    
    if (isa(obj, 'coneMosaicHex'))
        photocurrents = zeros(nTrials, numel(nonNullConesIndices), numel(eyeMovementTimeAxis), 'single');
        for ii=1:nTrials
            obj.absorptions = squeeze(absorptions(ii,:,:));
            obj.computeCurrent;
            photocurrents(ii,:,:) = single(reshape(obj.current, [1 numel(nonNullConesIndices) numel(eyeMovementTimeAxis)]));
        end
    else
        photocurrents  = zeros(nTrials, obj.rows, obj.cols, numel(eyeMovementTimeAxis), 'single');
        % compute the photocurrent for all of the trials
        % coneMosaic.computeCurrent at a later time, rather than set this flag.
        % @BW:  I think it might be better to replace this code with that call
        % so that we only have one computeCurrent method.
        for ii=1:nTrials
            obj.absorptions = squeeze(absorptions(ii,:,:,:));
            obj.computeCurrent;
            photocurrents(ii,:,:,:) = single(reshape(obj.current, [1 obj.rows obj.cols numel(eyeMovementTimeAxis)]));
        end
    end
    
    % varargout{1} = photocurrents;
    
    % Time axis - No longer needed.
    %     dtOS = obj.os.timeStep;
    %     osTimeAxis = absorptionsTimeAxis(1): dtOS : absorptionsTimeAxis(end);
    %     varargout{2} = osTimeAxis;
end     


% Reload the absorptions signal from the last instance (again, since we
% destroyed obj.absorptions in the current computation)
if (isa(obj, 'coneMosaicHex')) 
   % Reshape to full hex pattern
   obj.absorptions = oneDhexTo2Dhex(squeeze(absorptions(nTrials,:,:)), obj.pattern);
else
    % Reshape absorptions to correct dimensions [instances, cone_rows, cone_cols, time]
    absorptions = reshape(absorptions, [nTrials size(obj.pattern,1) size(obj.pattern,2) numel(eyeMovementTimeAxis)]);

    % If we don't compute the current, then we just return the absorptions from the last trial.
    obj.absorptions = reshape(squeeze(absorptions(nTrials,:,:)), [size(obj.pattern,1) size(obj.pattern,2) numel(eyeMovementTimeAxis)]);
end
        

% Function to reformat absorptions
function reformatAbsorptionsAllTrialsMatrix(nTrials, timePointsNum, coneRows, coneCols)
    % Save all absorptions as singles to save space
    absorptionsAllTrials = single(absorptionsAllTrials);
    % Reshape to cones x instances x timePoints. Note the 3rd dimension of
    % absorptionsAllTrials is traditionally time, but when we are doing multiple instances, it is instances * time
    absorptionsAllTrials = reshape(absorptionsAllTrials, [coneRows*coneCols nTrials timePointsNum]);
    % Only get the absorptions for the non-null cones - this has an effect only for coneMosaicHex mosaics
    absorptionsAllTrials = absorptionsAllTrials(nonNullConesIndices,:,:);
    % Reshape to [instances x cones x timePoints]
    absorptionsAllTrials = permute(absorptionsAllTrials, [2 1 3]);
end
 
end

function twoDhex = oneDhexTo2Dhex(oneDhex, pattern)
    
    nonNullCones = pattern(find(pattern>1));

    iLsource = find(nonNullCones==2);
    iMsource = find(nonNullCones==3);
    iSsource = find(nonNullCones==4);
        
    iLdest = find(pattern==2);
    iMdest = find(pattern==3);
    iSdest = find(pattern==4);
        
    twoDhex = zeros(size(obj.pattern,1)*size(obj.pattern,2), size(oneDhex,2));
    twoDhex(iLdest,:) = oneDhex(iLsource,:);
    twoDhex(iMdest,:) = oneDhex(iMsource,:);
    twoDhex(iSdest,:) = oneDhex(iSsource,:); 
    twoDhex = reshape(twoDhex, size(obj.pattern,1), size(obj.pattern,2), size(oneDhex,2));
end

function displayProgress(workerID, workDescription, progress)
    
    maxStarsNum = 60;
    if (isnan(progress))
        fprintf('worker-%02d: %s |', workerID, workDescription);
        for k = 1:60
            fprintf('*');
        end
        fprintf('|');
    else
        fprintf('worker-%02d: %s |', workerID, workDescription);
        if (progress>1)
            progress = 1;
        end
        starsNum = round(maxStarsNum*progress);
        for k = 1:starsNum
            fprintf('*');
        end
    end
    fprintf('\n');
end

