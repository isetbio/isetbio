% Compute the pattern of cone absorptions and typically the photocurrent
% for a sequence of optica images (@oiSequence), and optionally a number of
% eye movement paths.
%
%    [absorptions, absorpionsTimeAxis, [current, currentTimeAxis]] = ...
%         cMosaic.compute(oiSequence, varargin);
%
% Inputs:
%   oiSequence  - an @oiSequence object
%
% Optional inputs:
%   emPaths      - N eye movement paths, each with M eye movements, organized in an [N x M x 2] matrix
%
%   currentFlag  - logical, whether to compute photocurrent
%   newNoise     - logical, whether to use new random seed
%
% Outputs:
%   absorptions          - cone photon absorptions (photon counts in integrationTime)
%   absorptionsTimeAxis  - time axis for the absorptions signal
%   photocurrent         - cone photocurrent
%   photocurrentTimeAxis - time axis for photocurrent signal 
%
% NPC ISETBIO Team 2016
%

function [absorptions, absorptionsTimeAxis, varargout] = computeForOISequence(obj, oiSequence, varargin)

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
    
    % Compute eye movement time axis
    instancesNum = size(emPaths,1);
    eyeMovementsNum = size(emPaths,2);
    eyeMovementTimeAxis = oiTimeAxis(1) + (0:1:(eyeMovementsNum-1)) * obj.integrationTime;
    
    % Compute OIrefresh
    oiRefreshInterval = oiTimeAxis(2)-oiTimeAxis(1);
    
    % Save default integration time
    defaultIntegrationTime = obj.integrationTime;

    % Allocate memory for absorptions in [instancesNum x time x cone_rows x
    % cone_cols] format. This is done because Matlab drops the last
    % dimension if it is singleton, i.e., when there is 1 eye movement
    % At then end we reshape it to [instancesNum x cone_rows x cone_cols x time]
    absorptions = zeros(instancesNum, numel(eyeMovementTimeAxis), size(obj.pattern,1), size(obj.pattern,2), 'single');
    
    if (oiRefreshInterval >= defaultIntegrationTime)

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
                error('empty indices');
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
            
            totalAbsorptions = absorptionsDuringPreviousFrame+absorptionsDuringCurrentFrame;
            % insert the sum of the two partial absorptions in the time series  
            insertionIndex = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
            absorptions(1:instancesNum, insertionIndex, :, :) = single(permute(reshape(totalAbsorptions, [1 size(totalAbsorptions)]), [4 1 2 3]));
            
            
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
                absorptionsForRemainingEyeMovements = reshape(absorptionsForRemainingEyeMovements, [size(absorptionsForRemainingEyeMovements,1) size(absorptionsForRemainingEyeMovements,2) instancesNum numel(idx)]);
                absorptionsForRemainingEyeMovements = permute(absorptionsForRemainingEyeMovements, [3 4 1 2]);
                
                % insert in time series  
                insertionIndices = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
                absorptions(1:instancesNum, insertionIndices,:,:) = single(absorptionsForRemainingEyeMovements); 
            end
        end  % oiIndex

    % oiRefreshInterval > defaultIntegrationTime
    else
        
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
                error('empty indices');
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
           
           % insert to time series  
           insertionIndices = round((eyeMovementTimeAxis(emIndex)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
           absorptions(1:instancesNum, insertionIndices, :, :) = single(permute(reshape(absorptionsAccum, [1 size(absorptionsAccum)]), [4 1 2 3]));
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
        absorptions = reshape(absorptions, [size(absorptions,1) size(absorptions,2) size(absorptions,3)*size(absorptions,4)]);
        nonNullConeIndices = find(obj.pattern > 1);
        absorptions = absorptions(:,:,nonNullConeIndices);
        % Reshape absorptions to correct dimensions [instances, cone_indices, time]
        absorptions = permute(absorptions, [1 3 2]);
        if (currentFlag)
            % Add one more absorption at the end
            absorptions = cat(3, absorptions, squeeze(absorptions(:,:,end)));
        end
    else
        % Reshape absorptions to correct dimensions [instances, cone_rows, cone_cols, time]
        absorptions = permute(absorptions, [1 3 4 2]);
        if (currentFlag)
            % Add one more absorption at the end
            absorptions = cat(4, absorptions, absorptions(:,:,:,end));
        end
    end

    % Reload the mean(over all instances) absorptions signal
    obj.absorptions = mean(absorptions,1);
    
    % align absorptions time axis with respect to optical image sequence time axis
    absorptionsTimeAxis = oiTimeAxis(1) + obj.absorptionsTimeAxis; 

    if (currentFlag)
        % compute the photocurrent time axis
        dtOS = obj.os.timeStep;
        osTimeAxis = absorptionsTimeAxis(1): dtOS : absorptionsTimeAxis(end);
        
        if (isa(obj, 'coneMosaicHex'))
            photocurrents = zeros(instancesNum, numel(nonNullConeIndices), numel(osTimeAxis), 'single');
            for instanceIndex = 1:instancesNum
                fprintf('Computing photocurrents for instance %d/%d\n', instanceIndex,instancesNum);
                tmp = squeeze(absorptions(instanceIndex,:,:));
                
                % Resample to osTimeAxis
                tmp = coneMosaic.tResample(tmp, obj.pattern(nonNullConeIndices), absorptionsTimeAxis, osTimeAxis);

                % osCompute expects a 3D pRate, so make it so
                tmp = reshape(tmp, [size(tmp,1) 1 size(tmp,2)]);
                % Compute photocurrent from photonRate (tmp/dtOS)
                tmp = single(obj.os.osCompute(tmp/dtOS, obj.pattern(nonNullConeIndices), 'append', false));
                % Put it back in correct shape
                photocurrents(instanceIndex,:,:) = permute(tmp, [2 1 3]);
            end
        else
            photocurrents = zeros(instancesNum, size(obj.pattern,1), size(obj.pattern,2), numel(osTimeAxis), 'single');
            for instanceIndex = 1:instancesNum
                fprintf('Computing photocurrents for instance %d/%d\n', instanceIndex,instancesNum);
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
        
        % Reload the mean(over all instances) absorptions signal
        obj.absorptions = mean(absorptions,1);
    
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