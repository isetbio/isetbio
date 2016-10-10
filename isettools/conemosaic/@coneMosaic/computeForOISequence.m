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
    
    % Check that there is at least 1 eye movement per OI
    eyeMovementsNumPerOpticalImage = (oiTimeAxis(2)-oiTimeAxis(1))/obj.integrationTime;
    if (eyeMovementsNumPerOpticalImage < 1.000-1e-6)
        error('\nEye movements per optical image (%f) is less than 1. Integration time: %g, stimulus interval: %g.\nEither decrease the integrationTime or increase the stimulusSamplingInterval \n', eyeMovementFramesPerOpticalImage, obj.integrationTime, oiTimeAxis(2)-oiTimeAxis(1));
    end
    
    % Initialize our time series
    absorptions = []; absorptionsTimeAxis = []; eyeMovementSequence = [];  
    lastEyeMovementIndex = 0;
    
    % Loop over the optical images and compute isomerizations
    for oiIndex = 1:numel(oiSequence)
        
        % Retrieve eye movements for current OI
        firstEyeMovementIndex = lastEyeMovementIndex+1;
        lastEyeMovementIndex = round(oiIndex*eyeMovementsNumPerOpticalImage);
        eyeMovementIndices = firstEyeMovementIndex:lastEyeMovementIndex;
        % fprintf('Eye movements num for optical image %d: %2.2f\n', oiIndex, numel(eyeMovementIndices));
        
        eyeMovementPathForThisOI = eyeMovementsForOISequence(eyeMovementIndices,:);
        
        % Compute absorptions for current OI and eye movement path
        absorptionsForThisOI = obj.compute(oiSequence{oiIndex}, ...
            'emPath', eyeMovementPathForThisOI, ...      % current OI eye movement path
            'newNoise', newNoise, ...
            'currentFlag', false ...                     % no current computation for each oi -current will be computed for the entire sequence of absorptions
            );
        
        % Concatenate sequences
        if (isempty(absorptionsTimeAxis))
            absorptionsTimeAxis = obj.absorptionsTimeAxis;
            eyeMovementSequence = eyeMovementPathForThisOI;
            absorptions = absorptionsForThisOI;
        else
            absorptionsTimeAxis = cat(2, absorptionsTimeAxis, absorptionsTimeAxis(end)+obj.absorptionsTimeAxis);
            absorptions = cat(3, absorptions, absorptionsForThisOI);
            eyeMovementSequence = cat(1, eyeMovementSequence, eyeMovementPathForThisOI);
        end
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
        resampledAbsorptionsSequence = coneMosaic.resampleAbsorptionsSequence(absorptions, absorptionsTimeAxis, osTimeAxis);
    
        % Convert to photon rate in photons/sec for the osTimeStep
        pRate = resampledAbsorptionsSequence/dtOS;
        varargout{1} = obj.os.osCompute(pRate, obj.pattern, 'append', false);
    
        % Compute photocurrent time axis (starting at the origin of the oiTimeAxis)
        varargout{2} = osTimeAxis;
    end % currentFlag

    
end
