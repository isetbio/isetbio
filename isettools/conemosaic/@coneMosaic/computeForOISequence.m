% Compute the pattern of cone absorptions and typically the
% photocurrent
%    [absorptions, absorpionsTimeAxis, [current, currentTimeAxis]] = ...
%         cMosaic.compute(oiSequence, oiTimeAxis);
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

function [absorptions, absorptionsTimeAxis, varargout] = computeForOISequence(obj, oiSequence, oiTimeAxis, varargin)

    p = inputParser;
    p.addRequired('oiSequence', @(x)isa(x, 'oiSequence'));
    p.addRequired('oiTimeAxis', @isnumeric);
    p.addParameter('currentFlag', false, @islogical);
    p.addParameter('newNoise', true, @islogical);
    p.parse(oiSequence, oiTimeAxis,varargin{:});
    
    oiSequence = p.Results.oiSequence;
    oiTimeAxis = p.Results.oiTimeAxis;
    currentFlag = p.Results.currentFlag;
    newNoise = p.Results.newNoise;

    if (oiSequence.length ~= numel(oiTimeAxis))
        error('oiTimeAxis and oiSequence must have equal length\n');
    end
    
    % Save a copy of the entire eye movement sequence because we will be 
    % applying different emPath segments for different optical images
    eyeMovementsForOISequence = obj.emPositions;
    eyeMovementTimeAxis = oiTimeAxis(1) + (0:1:(size(eyeMovementsForOISequence,1)-1)) * obj.integrationTime;
    
    timingPrecisionDigits = 7;
    % Round eyeMovementTimeAxis
    eyeMovementTimeAxis(eyeMovementTimeAxis>=0) = round(eyeMovementTimeAxis(eyeMovementTimeAxis>=0), timingPrecisionDigits);
    eyeMovementTimeAxis(eyeMovementTimeAxis<0) = -round(-eyeMovementTimeAxis(eyeMovementTimeAxis<0), timingPrecisionDigits);
    
    % Round oiTimeAxis 
    oiTimeAxis = sign(oiTimeAxis) .* round(abs(oiTimeAxis), timingPrecisionDigits);
    
    % Round oiRefresh and save a copy of the default integrationTime
    oiRefreshInterval = round(oiTimeAxis(2)-oiTimeAxis(1), timingPrecisionDigits);
    defaultIntegrationTime = round(obj.integrationTime, timingPrecisionDigits);
    
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
%                                       p1   p2     full       full
%                               partial_/     \_partial
%

        % Loop over the optical images
        for oiIndex = 1:oiSequence.length

            % Current oi time limits
            tFrameStart = oiTimeAxis(oiIndex);
            tFrameEnd   = tFrameStart + oiRefreshInterval;

            % Find eye movement indices withing the oi limits
            indices = find( (eyeMovementTimeAxis >  tFrameStart-defaultIntegrationTime) & ...      
                            (eyeMovementTimeAxis <= tFrameEnd-defaultIntegrationTime+eps) );
            
            % the first eye movement requires special treatment as it may have started before the current frame,
            % so we need to compute partial absorptions over the previous frame and over the current frame
            idx = indices(1);
            integrationTimeForFirstPartialAbsorption = tFrameStart-eyeMovementTimeAxis(idx);
            integrationTimeForSecondPartialAbsorption = eyeMovementTimeAxis(idx)+defaultIntegrationTime-tFrameStart;

            % Partial absorptions (p1 in graph above) with previous oi
            if (oiIndex > 1)
                % Update coneMosaic before compute()
                obj.integrationTime = integrationTimeForFirstPartialAbsorption;
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
            % Update coneMosaic before compute()
            obj.integrationTime = integrationTimeForSecondPartialAbsorption;
            absorptionsDuringCurrentFrame =  obj.compute(...
                                oiSequence.frameAtIndex(oiIndex), ...
                                'emPath', eyeMovementsForOISequence(idx,:), ...
                                'newNoise', newNoise, ...
                                'currentFlag', false ...                                                    
                                );
            % insert the sum in the to time series  
            insertionIndex = round((eyeMovementTimeAxis(idx)-eyeMovementTimeAxis(1))/defaultIntegrationTime)+1;
            absorptions(:,:, insertionIndex) = absorptionsDuringPreviousFrame+absorptionsDuringCurrentFrame;

            % Full absorptions with current oi and default integration time)
            if (numel(indices)>1)
                % Update coneMosaic before compute()
                idx = indices(2:end);
                obj.integrationTime = defaultIntegrationTime;
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

        % Loop over the eye movements
        for emIndex = 1:numel(eyeMovementTimeAxis) 
        
            % Current eye movement time limits
            emStart = eyeMovementTimeAxis(emIndex);
            emEnd   = emStart + defaultIntegrationTime;
            
            actualIntegrationTime = 0;
            
            % Find oi indices withing the eye movement frame time limits
            indices = find( (oiTimeAxis > emStart-oiRefreshInterval) & ...
                            (oiTimeAxis <= emEnd + eps) );
            
            % Partial absorptions during the ovelap with the oi that started before the emStart
            idx = indices(1);
            integrationTimeForFirstPartialAbsorption = oiTimeAxis(idx)+oiRefreshInterval-emStart;
            % Update coneMosaic before compute()
            obj.integrationTime = integrationTimeForFirstPartialAbsorption;
            actualIntegrationTime = actualIntegrationTime + obj.integrationTime;
            absorptionsAccum = obj.compute(...
                                oiSequence.frameAtIndex(idx), ...
                                'emPath', eyeMovementsForOISequence(emIndex,:), ...
                                'newNoise', newNoise, ...
                                'currentFlag', false ...                                                    
                                ); 
 
            % Next, compute full absorptions for the remaining OIs
            if (numel(indices)>1) 
                for k = 2:numel(indices)
                    % Update coneMosaic before compute()
                    if (k < numel(indices))
                        obj.integrationTime = oiRefreshInterval;
                    else
                        idx = indices(end);
                        integrationTimeForLastPartialAbsorption = emEnd - oiTimeAxis(idx);
                        obj.integrationTime = integrationTimeForLastPartialAbsorption;
                    end
                    actualIntegrationTime = actualIntegrationTime + obj.integrationTime;
                    absorptionsAccum = absorptionsAccum + obj.compute(...
                            oiSequence.frameAtIndex(indices(k)), ...
                            'emPath', eyeMovementsForOISequence(emIndex,:), ...      
                            'newNoise', newNoise, ...
                            'currentFlag', false ...      
                            );
                end  % for k              
            end
           
           if (abs(actualIntegrationTime-defaultIntegrationTime) > defaultIntegrationTime*0.01) && (emIndex < numel(eyeMovementTimeAxis))
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
    
        % Resample absorptions to osTimeAxis
        resampledAbsorptionsSequence = coneMosaic.tResample(absorptions, absorptionsTimeAxis, osTimeAxis);
    
        % Convert to photon rate in photons/sec for the osTimeStep
        pRate = resampledAbsorptionsSequence/dtOS;
        
        % Compute and return the photocurrent response
        varargout{1} = obj.os.osCompute(pRate, obj.pattern, 'append', false);
        
        % Return the photocurrent time axis 
        varargout{2} = osTimeAxis;
    end % currentFlag
end