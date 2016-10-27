% Compute the pattern of cone absorptions and typically the
% photocurrent
%    [absorptions, absorpionsTimeAxis, [current, currentTimeAxis]] = ...
%         cMosaic.compute(oiSequence, varargin);
%
% Inputs:
%   oiSequence  - an @oiSequence object
%   oiTimeAxis  - time axis for the optical image sequence
%
% Optional inputs:
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

function [absorptions, absorptionsTimeAxis, varargout] = computeForOISequence(obj, oiSequence,  varargin)

    p = inputParser;
    p.addRequired('oiSequence', @(x)isa(x, 'oiSequence'));
    p.addParameter('currentFlag', false, @islogical);
    p.addParameter('newNoise', true, @islogical);
    p.parse(oiSequence,  varargin{:});
    
    oiSequence = p.Results.oiSequence;
    currentFlag = p.Results.currentFlag;
    newNoise = p.Results.newNoise;
    oiTimeAxis = oiSequence.oiTimeAxis;
    
    if (oiSequence.length ~= numel(oiTimeAxis))
        error('oiTimeAxis and oiSequence must have equal length\n');
    end
    
    % Save a copy of the entire eye movement sequence because we will be 
    % applying different emPath segments for different optical images
    eyeMovementsForOISequence = obj.emPositions;
    eyeMovementTimeAxis = oiTimeAxis(1) + (0:1:(size(eyeMovementsForOISequence,1)-1)) * obj.integrationTime;
    
    % Compute OIrefresh
    oiRefreshInterval = oiTimeAxis(2)-oiTimeAxis(1);
    
    % Save default integration time
    defaultIntegrationTime = obj.integrationTime;

    % Initialize our time series
    absorptions = zeros(size(obj.pattern,1), size(obj.pattern,2), numel(eyeMovementTimeAxis));
    
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
            if (oiIndex > 1)
                % Update the @coneMosaic with the partial integration time
                obj.integrationTime = integrationTimeForFirstPartialAbsorption;
                % Compute partial absorptions
                absorptionsDuringPreviousFrame = obj.compute(...
                                oiSequence.frameAtIndex(oiIndex-1), ...
                                'emPath', eyeMovementsForOISequence(idx,:), ...
                                'newNoise', newNoise, ...
                                'currentFlag', false ...                                                    
                                ); 
            else
                absorptionsDuringPreviousFrame = zeros(size(obj.pattern,1), size(obj.pattern,2), 1);
            end

            % Partial absorptions (p2 in graph above) with current oi
            % Update the @coneMosaic with the partial integration time
            obj.integrationTime = integrationTimeForSecondPartialAbsorption;
            % Compute partial absorptions
            absorptionsDuringCurrentFrame =  obj.compute(...
                                oiSequence.frameAtIndex(oiIndex), ...
                                'emPath', eyeMovementsForOISequence(idx,:), ...
                                'newNoise', newNoise, ...
                                'currentFlag', false ...                                                    
                                );
            % insert the sum of the two partial absorptions in the time series  
            insertionIndex = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
            absorptions(:,:, insertionIndex) = absorptionsDuringPreviousFrame+absorptionsDuringCurrentFrame;

            % Full absorptions with current oi and default integration time)
            if (numel(indices)>1)    
                % Update the @coneMosaic with the default integration time
                obj.integrationTime = defaultIntegrationTime;
                % Compute absorptions for all remaining the OIs
                idx = indices(2:end);
                absorptionsForRemainingEyeMovements = obj.compute(...
                        oiSequence.frameAtIndex(oiIndex), ...
                        'emPath', eyeMovementsForOISequence(idx,:), ...      
                        'newNoise', newNoise, ...
                        'currentFlag', false ...      
                        );
                % insert in time series  
                insertionIndices = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
                absorptions(:,:, insertionIndices) = absorptionsForRemainingEyeMovements; 
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
%
        % Loop over the eye movements
        for emIndex = 1:numel(eyeMovementTimeAxis) 
        
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
            absorptionsAccum = obj.compute(...
                                oiSequence.frameAtIndex(idx), ...
                                'emPath', eyeMovementsForOISequence(emIndex,:), ...
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
                    absorptionsAccum = absorptionsAccum + obj.compute(...
                            oiSequence.frameAtIndex(indices(k)), ...
                            'emPath', eyeMovementsForOISequence(emIndex,:), ...      
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
           absorptions(:,:, insertionIndices) = absorptionsAccum;
        end % emIndex  
    end

    % Restore default integrationTime
    obj.integrationTime = defaultIntegrationTime;

    % Reload the full eye movement sequence
    obj.emPositions = eyeMovementsForOISequence;
    
    % Reload the full absorptions signal
    obj.absorptions = absorptions;
    
    % align absorptions time axis with respect to optical image sequence time axis
    absorptionsTimeAxis = oiTimeAxis(1) + obj.absorptionsTimeAxis; 
    
    if (currentFlag)
        % compute the photocurrent time axis
        dtOS = obj.os.timeStep;
        osTimeAxis = absorptionsTimeAxis(1): dtOS :absorptionsTimeAxis(end);
    
        % Resample absorptions to osTimeAxis timebase
        resampledAbsorptionsSequence = coneMosaic.tResample(absorptions, absorptionsTimeAxis, osTimeAxis);
    
        % Convert to photon rate in photons/sec for the osTimeStep
        pRate = resampledAbsorptionsSequence/dtOS;
        
        % Compute and return the photocurrent response
        varargout{1} = obj.os.osCompute(pRate, obj.pattern, 'append', false);
        
        % Return the photocurrent time axis 
        varargout{2} = osTimeAxis;
    end % currentFlag
end