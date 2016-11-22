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
%   workerID     - if passed and is non-empty, display the progress for the current worker in the parpool pool
%
% Outputs:
%   absorptions          - cone photon absorptions (photon counts in integrationTime)
%   absorptionsTimeAxis  - time axis for the absorptions signal
%   photocurrent         - cone photocurrent
%   photocurrentTimeAxis - time axis for photocurrent signal 
%
% NPC ISETBIO Team 2016
%

function [absorptions, absorptionsTimeAxis, varargout] = computeForOISequenceV2(obj, oiSequence, varargin)

    p = inputParser;
    p.addRequired('oiSequence', @(x)isa(x, 'oiSequence'));
    p.addParameter('emPaths', [], @isnumeric);
    p.addParameter('currentFlag', false, @islogical);
    p.addParameter('newNoise', true, @islogical);
    p.addParameter('workerID', [], @isnumeric);
    p.addParameter('workDescription', '', @ischar);
    p.parse(oiSequence, varargin{:});
    
    oiSequence = p.Results.oiSequence;
    emPaths = p.Results.emPaths;
    currentFlag = p.Results.currentFlag;
    newNoise = p.Results.newNoise;
    workerID = p.Results.workerID;
    workDescription = p.Results.workDescription;
    
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
    
    % Save default integration time
    defaultIntegrationTime = obj.integrationTime;
    
    % Compute OIrefresh
    if (numel(oiTimeAxis) == 1)
        oiRefreshInterval = defaultIntegrationTime;
    else
        oiRefreshInterval = oiTimeAxis(2)-oiTimeAxis(1);
    end
    
    % Only allocate memory for the non-null cones in a 3D matrix [instances x numel(nonNullConesIndices) x time]
    nonNullConesIndices = find(obj.pattern>1);
    absorptions = zeros(instancesNum, numel(nonNullConesIndices), numel(eyeMovementTimeAxis), 'single');
    
    
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

            if (~isempty(workerID))
                % Update progress bar
                displayProgress(workerID, workDescription, 0.5*oiIndex/oiSequence.length);
            end
            
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
                obj.absorptions = [];
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
            obj.absorptions = [];
            absorptionsDuringCurrentFrame =  obj.compute(...
                                oiSequence.frameAtIndex(oiIndex), ...
                                'emPath', emSubPath, ...
                                'newNoise', newNoise, ...
                                'currentFlag', false ...                                                    
                                );
            
            % summed absorptions
            absorptionInstances = absorptionsDuringPreviousFrame+absorptionsDuringCurrentFrame;
                
            % Reformat and insert to time series  
            insertionIndices = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1; 
            reformatAbsorptionInstancesMatrix(instancesNum, numel(insertionIndices), size(obj.pattern,1), size(obj.pattern,2));
            absorptions(1:instancesNum, :, insertionIndices) = absorptionInstances;
      
            % Full absorptions with current oi and default integration time)
            if (numel(indices)>1)    
                % Update the @coneMosaic with the default integration time
                obj.integrationTime = defaultIntegrationTime;
                % Compute absorptions for all remaining the OIs
                idx = indices(2:end);
                emSubPath = reshape(emPaths(1:instancesNum, idx,:), [instancesNum*numel(idx) 2]);
                obj.absorptions = [];
                absorptionInstances = obj.compute(...
                        oiSequence.frameAtIndex(oiIndex), ...
                        'emPath', emSubPath, ...      
                        'newNoise', newNoise, ...
                        'currentFlag', false ...      
                        );
                
                % Reformat and insert to time series     
                insertionIndices = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
                reformatAbsorptionInstancesMatrix(instancesNum, numel(insertionIndices), size(obj.pattern,1), size(obj.pattern,2));
                absorptions(1:instancesNum, :, insertionIndices) = absorptionInstances;
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
            obj.absorptions = [];
            absorptionInstances = obj.compute(...
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
                    obj.absorptions = [];
                    absorptionInstances = absorptionInstances + obj.compute(...
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
           reformatAbsorptionInstancesMatrix(instancesNum, numel(insertionIndices), size(obj.pattern,1), size(obj.pattern,2));
           absorptions(1:instancesNum, :, insertionIndices) = absorptionInstances;
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
        if (currentFlag)
            % Add one more absorption at the end
            absorptions = cat(3, absorptions, squeeze(absorptions(:,:,end)));
        end
    else
        % Reshape absorptions to correct dimensions [instances, cone_rows, cone_cols, time]
        absorptions = reshape(absorptions, [instancesNum size(obj.pattern,1) size(obj.pattern,2) numel(eyeMovementTimeAxis)]);
        
        if (currentFlag)
            % Add one more absorption at the end
            absorptions = cat(4, absorptions, absorptions(:,:,:,end));
        end
    end

    % Reload the absorptions signal from the last instance
    obj.absorptions = squeeze(absorptions(end,:,:,:));
    
    % align absorptions time axis with respect to optical image sequence time axis
    absorptionsTimeAxis = oiTimeAxis(1) + obj.absorptionsTimeAxis; 

    % Special case where we only have a time series with just 1 point (+the
    % extra time point inserted above)
    if (numel(absorptionsTimeAxis) == 2)
        % Remove the last absorption we inserted at the end
        if (isa(obj, 'coneMosaicHex')) 
            absorptions = absorptions(:,:,1:(end-1));
        else
            absorptions = absorptions(:,:,:,1:(end-1));
        end
        
        varargout{1} = [];
        varargout{2} = [];
        return;
    end
    
    if (currentFlag)
        % compute the photocurrent time axis
        dtOS = obj.os.timeStep;
        osTimeAxis = absorptionsTimeAxis(1): dtOS : absorptionsTimeAxis(end);
        
        if (isa(obj, 'coneMosaicHex'))
            photocurrents = zeros(instancesNum, numel(nonNullConesIndices), numel(osTimeAxis), 'single');
            for instanceIndex = 1:instancesNum
                
                if (~isempty(workerID))
                    displayProgress(workerID, workDescription, 0.5 + 0.5*instanceIndex/instancesNum);
                end
                
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
                
                if (~isempty(workerID))
                    displayProgres(workerID, 0.5 + 0.5*instanceIndex/instancesNum);
                end
                
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
        instancesNum = size(absorptions,1);
        timePoints = size(absorptions,4);
        if (isa(obj, 'coneMosaicHex')) 
            obj.absorptions = squeeze(absorptions(instancesNum,:,:,:));
        else
            obj.absorptions = reshape(squeeze(absorptions(instancesNum,:,:,:)), [size(obj.pattern,1) size(obj.pattern,2) timePoints]);
        end
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
 
    if (~isempty(workerID))
        displayProgress(workerID, workDescription, nan);
    end
    
    % Function to reformat absorptions
    function reformatAbsorptionInstancesMatrix(instancesNum, timePointsNum, coneRows, coneCols)
        % Will save all absorptions as singles
        absorptionInstances = single(absorptionInstances);
        % Reshape to cones x instances x timePoints. Note the 3rd dimension of
        % absorptionInstances is traditionally time, but here it is instances * time
        absorptionInstances = reshape(absorptionInstances, [coneRows*coneCols instancesNum timePointsNum]);
        % Only get the absorptions for the non-null cones - this has an
        % effect only for coneMosaicHex
        absorptionInstances = absorptionInstances(nonNullConesIndices,:,:);
        % Reshape to [instances x cones x timePoints]
        absorptionInstances = permute(absorptionInstances, [2 1 3]);
     end
        
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
                