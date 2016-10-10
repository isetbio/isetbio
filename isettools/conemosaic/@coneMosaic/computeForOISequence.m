% Compute the pattern of cone absorptions and typically the
% photocurrent
%    [absorptions, absorpionsTimeAxis, [current, currentTimeAxis]] = ...
%         cMosaic.compute(oiSequence, oiTimeAxis);
%
% Inputs:
%   oiSequence  - cell array with a sequence of optical images, see oiCreate for more details
%   oiTimeAxis  - time axis for the optical images
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
    
    % Initialize our time series
    absorptions = [];
    
    % Loop over the optical images and compute isomerizations
    for oiIndex = 1:numel(oiSequence)
        
        % Retrieve eye movements for current OI
        if (oiIndex < numel(oiSequence))
            eyeMovementIndices = find((eyeMovementTimeAxis >= oiTimeAxis(oiIndex)) & (eyeMovementTimeAxis < oiTimeAxis(oiIndex+1)));
        else
            eyeMovementIndices = find(eyeMovementTimeAxis >= oiTimeAxis(oiIndex));
        end
        
        fprintf('Eye movements num for optical image %d: %2.2f\n', oiIndex, numel(eyeMovementIndices));
        
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
       
    % Reload the full eye movement sequence
    obj.emPositions = eyeMovementsForOISequence;
    
    % Reload the full absorptions signal
    obj.absorptions = absorptions;
    
    % align absorptions time axis with respect to optical image sequence time axis
    absorptionsTimeAxis = oiTimeAxis(1) + obj.absorptionsTimeAxis; 
    
    if (currentFlag)
        % compute the os time axis
        dtOS = obj.os.timeStep;
        osTimeAxis = absorptionsTimeAxis(1): dtOS :absorptionsTimeAxis(end);
    
        % Resample absorptions to osTimeAxis
        resampledAbsorptionsSequence = coneMosaic.resample(absorptions, absorptionsTimeAxis, osTimeAxis);
    
        % Convert to photon rate in photons/sec for the osTimeStep
        pRate = resampledAbsorptionsSequence/dtOS;
        varargout{1} = obj.os.osCompute(pRate, obj.pattern, 'append', false);
    
        % Compute photocurrent time axis (starting at the origin of the oiTimeAxis)
        varargout{2} = osTimeAxis;
    end % currentFlag

    
end
