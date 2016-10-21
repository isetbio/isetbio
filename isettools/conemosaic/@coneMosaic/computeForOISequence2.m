% Compute the pattern of cone absorptions and typically the
% photocurrent
%    [absorptions, absorpionsTimeAxis, [current, currentTimeAxis]] = ...
%         cMosaic.compute(oiSequence, oiTimeAxis);
%
% Inputs:
%   oiSequence  - cell array with a sequence of optical images, see oiCreate for more details
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

function [absorptions, absorptionsTimeAxis, varargout] = computeForOISequence2(obj, oiSequence, oiTimeAxis, varargin)

    p = inputParser;
    p.addRequired('oiSequence',@iscell);
    p.addRequired('oiTimeAxis',@isnumeric);
    p.addParameter('currentFlag', false, @islogical);
    p.addParameter('newNoise', true, @islogical);
    p.parse(oiSequence, oiTimeAxis,varargin{:});
    
    oiSequence = p.Results.oiSequence;
    oiTimeAxis = p.Results.oiTimeAxis;
    currentFlag = p.Results.currentFlag;
    newNoise = p.Results.newNoise;

    if (numel(oiSequence) ~= numel(oiTimeAxis))
        error('oiTimeAxis and oiSequence must have equal length\n');
    end
    
    % Save a copy of the entire eye movement sequence because we will be 
    % applying different emPath segments for different optical images
    eyeMovementsForOISequence = obj.emPositions;
    eyeMovementTimeAxis = oiTimeAxis(1) + (0:1:(size(eyeMovementsForOISequence,1)-1)) * obj.integrationTime;
    
    size(eyeMovementTimeAxis)
    
    % Initialize our time series
    absorptions = [];
    
    oiRefreshInterval = oiTimeAxis(2)-oiTimeAxis(1);

    % Case where we have an integer number of eye movements / refresh interval    
    if (mod(round(oiRefreshInterval*1000000),round(obj.integrationTime*1000000)) == 0)
        fprintf('Integer\n');
        % Loop over the optical images and compute isomerizations
        for oiIndex = 1:numel(oiSequence)
            % Retrieve eye movements for current OI
            eyeMovementIndices = find((eyeMovementTimeAxis >= oiTimeAxis(oiIndex)) & (eyeMovementTimeAxis < oiTimeAxis(oiIndex) + oiRefreshInterval));
            if (isempty(eyeMovementIndices))
                continue;
            end

            % Compute absorptions for current OI and eye movement path
            absorptionsForThisOI = obj.compute(oiSequence{oiIndex}, ...
                'emPath', eyeMovementsForOISequence(eyeMovementIndices,:), ...      % current OI eye movement path
                'newNoise', newNoise, ...
                'currentFlag', false ...                                            % current computation not computed for each oi - current will be computed for the entire sequence of absorptions
                );

            % Concatenate sequences
            absorptions = cat(3, absorptions, absorptionsForThisOI);
        end % oiIndex
    
    % Case where we have a non-integer number of eye movements / refresh interval
    else

        fprintf('Non Integer\n');
        
        % Save default integration time for later restoration
        originalIntegrationTime = obj.integrationTime;
            
        % Subcase where oiRefreshInterval > obj.integrationTime
        if (oiRefreshInterval > originalIntegrationTime)
                 
            % Loop over the optical images
            for oiIndex = 1:numel(oiSequence) 
                % Retrieve eye movements for current OI
                eyeMovementIndices = find((eyeMovementTimeAxis >= oiTimeAxis(oiIndex)) & (eyeMovementTimeAxis < oiTimeAxis(oiIndex) + oiRefreshInterval));
                if (isempty(eyeMovementIndices))
                    continue;
                end
                   
                % Compute absorptions for current OI and ALL BUT THE LAST eye movement
                obj.integrationTime = originalIntegrationTime;
                absorptionsForThisOI = obj.compute(oiSequence{oiIndex}, ...
                        'emPath', eyeMovementsForOISequence(eyeMovementIndices(1:end-1),:), ...      % do not include first and last eye movement for current OI
                        'newNoise', newNoise, ...
                        'currentFlag', false ...                                                    
                        ); 
                % Concatenate sequences 
                absorptions = cat(3, absorptions, absorptionsForThisOI);
                disp('1');
                size(absorptions)
                
                % Compute absorptions for the LAST eye movement which could be extending into the next OI
                if (eyeMovementTimeAxis(eyeMovementIndices(end))+originalIntegrationTime >= oiTimeAxis(oiIndex) + oiRefreshInterval)

                    % Set new (partial) integrationTime applicable for the last eye movement duration into the current OI
                    obj.integrationTime = (oiTimeAxis(oiIndex) + oiRefreshInterval) - eyeMovementTimeAxis(eyeMovementIndices(end));
                    
                    absorptionsForThisPartialOI = obj.compute(oiSequence{oiIndex}, ...
                            'emPath', eyeMovementsForOISequence(eyeMovementIndices(end),:), ...      % only include last eye movement for current OI
                            'newNoise', newNoise, ...
                            'currentFlag', false ...                                                    
                            );
                        
                    % Check if the current OI is the last OI
                    if (oiIndex == numel(oiSequence)) 
                        % Concatenate by adding the current OI partial absorptions
                        absorptions = cat(3, absorptions, absorptionsForThisPartialOI);
                        disp('2');
                size(absorptions)
                    else
                        % Set new (partial) integrationTime applicable to the last eye movement duration into the next OI
                        obj.integrationTime = eyeMovementTimeAxis(eyeMovementIndices(end)) + originalIntegrationTime - oiTimeAxis(oiIndex+1);
                        
                        % Compute absorptions for next OI and the last eye movement
                        absorptionsForTheNextPartialOI = obj.compute(oiSequence{oiIndex+1}, ...
                            'emPath', eyeMovementsForOISequence(eyeMovementIndices(end),:), ...      % only include last eye movement for current OI
                            'newNoise', newNoise, ...
                            'currentFlag', false ...                                                    
                            );
                        
                        % Concatenate by adding the sum of partial absorptions for the current  OI + next  OI
                        absorptions = cat(3, absorptions, absorptionsForThisPartialOI + absorptionsForTheNextPartialOI);
                        disp('2');
                size(absorptions)
                    end
                else
                    % Compute absorptions for current OI
                    absorptionsForThisOI = obj.compute(oiSequence{oiIndex}, ...
                        'emPath', eyeMovementsForOISequence(eyeMovementIndices(end),:), ...
                        'newNoise', newNoise, ...
                        'currentFlag', false ...                                                    
                        );
                    % Concatenate sequences
                    absorptions = cat(3, absorptions, absorptionsForThisOI);
                    disp('2');
                 size(absorptions)
                end
            end
            
        % Subcase where obj.integrationTime > oiRefreshInterval
        else
            
            % Loop over the eye movements
            for eyeMovementIndex = 1:numel(eyeMovementTimeAxis)
                
                % Retrieve OIs during the duration of the current eye movement
                oiIndices = find((oiTimeAxis >= eyeMovementTimeAxis(eyeMovementIndex)-oiRefreshInterval) & (oiTimeAxis < eyeMovementTimeAxis(eyeMovementIndex) + originalIntegrationTime));
                if (isempty(oiIndices))
                    continue;
                end
                
                accumAbsorptions = 0;
                
                oiIndex = 1;
                % Set new (partial) integrationTime for the first OI
                obj.integrationTime = oiTimeAxis(oiIndices(oiIndex)) + oiRefreshInterval - eyeMovementTimeAxis(eyeMovementIndex);
         %       fprintf('int. time for point 1= %f, oiRefresh = %f eye t = %f\n', obj.integrationTime, oiRefreshInterval, eyeMovementTimeAxis(eyeMovementIndex));

                absorptionsForThisOI = obj.compute(oiSequence{oiIndices(oiIndex)}, ...
                        'emPath', eyeMovementsForOISequence(eyeMovementIndex,:), ...      
                        'newNoise', newNoise, ...
                        'currentFlag', false ...                                                    
                        );
                accumAbsorptions = accumAbsorptions + absorptionsForThisOI;
                     
                % Accumulate absorptions over the course of the OIs during
                % the current eye movement (except for the last OI)
                % Set new integrationTime equal to the OI interval
                obj.integrationTime = oiRefreshInterval;
                    
                for oiIndex = 2:numel(oiIndices)-1
                    absorptionsForThisOI = obj.compute(oiSequence{oiIndices(oiIndex)}, ...
                        'emPath', eyeMovementsForOISequence(eyeMovementIndex,:), ...    
                        'newNoise', newNoise, ...
                        'currentFlag', false ...                                                    
                        );
                     accumAbsorptions = accumAbsorptions + absorptionsForThisOI;
                end % oiIndex
                
                
                % Update accumAbsorptions with the absorptions during for the LAST OI which could be extending into the next eye movement
                % Set new (partial) integrationTime 
                oiIndex = numel(oiIndices);
                obj.integrationTime = eyeMovementTimeAxis(eyeMovementIndex) + originalIntegrationTime - oiTimeAxis(oiIndices(oiIndex));
        %        fprintf('int. time for point N = %f, oiRefresh = %f eye t = %f\n', obj.integrationTime, oiRefreshInterval, eyeMovementTimeAxis(eyeMovementIndex));

                absorptionsForThisOI = obj.compute(oiSequence{oiIndices(oiIndex)}, ...
                        'emPath', eyeMovementsForOISequence(eyeMovementIndex,:), ...  
                        'newNoise', newNoise, ...
                        'currentFlag', false ...                                                    
                        );
                accumAbsorptions = accumAbsorptions + absorptionsForThisOI;

                % Concatenate sequences 
                if (eyeMovementIndex == 1)
                    absorptions = zeros([size(accumAbsorptions) numel(eyeMovementTimeAxis)]);
                end
                absorptions(:,:,eyeMovementIndex) = accumAbsorptions; 
            end
        end
        
        % Restore original integration time
        obj.integrationTime = originalIntegrationTime;
            
    end
       
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
        resampledAbsorptionsSequence = coneMosaic.resample(absorptions, absorptionsTimeAxis, osTimeAxis);
    
        % Convert to photon rate in photons/sec for the osTimeStep
        pRate = resampledAbsorptionsSequence/dtOS;
        
        % Compute and return the photocurrent response
        varargout{1} = obj.os.osCompute(pRate, obj.pattern, 'append', false);
        
        % Return the photocurrent time axis 
        varargout{2} = osTimeAxis;
    end % currentFlag
end