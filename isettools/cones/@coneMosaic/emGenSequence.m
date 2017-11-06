function [nTrialsPos, nTrialsPosMicrons] = emGenSequence(obj, nFrames, varargin)
%EMGENSEQUENCE  Generate sequence of eye movements
%
%   [nTrialsPos, nTrialsPosMicrons] = EMGENSEQUENCE(obj,nFrames,'nTrials',1,'em',emCreate);
%
% The eye movement samples are created at the same temporal sample rate as the
% cone integration time.  We only update the position at the beginning of each
% integration time.
%
% This has some consequences for how we store and retrieve the eye movement
% parameters. Consider the case of tremor. We need to store the amplitudes of
% the tremor in a way that permits us to do the calculation correctly no matter
% what the cone integration time is. The literature specifies that the maximum
% frequency of tremor is 100 Hz. In that case, the typical sample time between
% each tremor is 10 ms. The amplitude of the tremor at  100 Hz is about 1 cone
% width (17e-3 radians/sample).
%
% Our calculations make a best estimate of the eye position at the start of
% every cone integration period.  We assume the eye is fixed throughout the
% integration period.  This approximation makes sense for calculations done at a
% shorter cone integration time than 10 ms, which is the tremor rate.
%
% For an example, suppose the update time is 10 ms. We store the amplitude in
% terms of radians per second, which is (17e-3)/(10 ms) = (17e-1/sec).  Thus,
% when you  look at the amplitude value, it may appear frighteningly large -
% 17e-1 radians is about 100 cones.  But remember this is the amplitude per
% second. When we calculate the tremor amplitude at the proper update rate,
% which is 10 ms, the amplitude is 100 times smaller (about 1 cone).
%
% In normal programming use, such as in this routine, dealing with the units is
% simplified for you because you can ask for the amplitude in units of
% cones/sample.  This takes into account the tremor parameters as well as the
% cone integration time. See the code below.
%
% Inputs:
%     obj        - cone mosaic object
%     nFrames    - number of frames to generate
%
% Ouputs:
%     nTrialsPos        - nTrials x nFrames x 2 matrix of eye positions in units of cone positions
%     nTrialsPosMicrons - nTrials x nFrames x 2 matrix of eye positions in units of microns
%
% Optional parameter name/value pairs chosen from the following:
%    'em'              Eye movement structure, see emCreate for details
%    'rSeed'           Random seed to be used
%    'nTrials'         Multiple trial case, default = 1
%
% The eye movement model combines tremor, drift and microsaccades.
% 
%   Examples:
%  
%    To control the eye movement parameters, set the eye movement parameters
%    structure, say ..
%
%      emParameters = emCreate;
%      <Set em parameters> 
%      coneMosaic.emGenSequence(nFrames,'nTrials',1,'em',emParameters);
%
% See also EMCREATE, EMSET, EMGET
%
% HJ/BW ISETBIO Team, 2016
% 11/06/17  ncp       Added line to make drift magnitude independent of
%                     sample time.
% 11/06/17  dhb, npc  Added comments on microsaccade algorithm.

%% parse input
p = inputParser;
p.addRequired('nFrames', @isscalar);
p.addParameter('nTrials',1,@isscalar);
p.addParameter('em', emCreate, @isstruct);
p.addParameter('rSeed', [], @isscalar);

% set parameters
sampTime  = obj.integrationTime;

p.parse(nFrames, varargin{:});
em = p.Results.em;
em = emSet(em, 'sample time', sampTime);
nTrials = p.Results.nTrials;

if ~isempty(p.Results.rSeed), rng(p.Results.rSeed); end
emFlag = emGet(em, 'em flag');
pos = zeros(nFrames, 2);

nTrialsPos = zeros(nTrials,nFrames,2);

% define cone parameters needed to convert units of mm or deg to cones, and vice versa
% pattern sample size in meters - for rect mosaics this is the same as the cone size, but not so for hex-mosaics
params.w = obj.patternSampleSize(1); 

% Are we only doing micro-saccades?
microSaccadesOnlyFlag = (emFlag(3) && (~emFlag(1)) && (~emFlag(2)));

for nn=1:nTrials
    %% generate eye movement for tremor
    if emFlag(1) || microSaccadesOnlyFlag 
        % Load parameters
        amplitude  = emGet(em, 'tremor amplitude', 'cones/sample', params);
        interval   = emGet(em, 'tremor interval');
        intervalSD = emGet(em, 'tremor interval SD');
        
        % Compute time of tremor occurs
        t = interval + randn(nFrames, 1) * intervalSD;
        t(t < 0.001) = 0.001; % get rid of negative values
        tPos = cumsum(t);
        tPos = max(1, round(tPos / sampTime));
        indx = 1:find(tPos <= nFrames, 1, 'last');
        tPos = tPos(indx);
        
        % Generate random step at the selected times
        direction = rand(length(tPos),1);
        
        % Unit length direction
        pos(tPos, :) = amplitude*[direction sqrt(1-direction.^2)];
        
        pos(tPos, :) = bsxfun(@times,pos(tPos,:), t(indx(:))/sampTime);
        pos = pos .* (2*(randn(size(pos))>0)-1); % shuffle the sign
        pos = cumsum(pos, 1);
    end
     
    % Generate eye movement for drift
    if emFlag(2)  ||  microSaccadesOnlyFlag  
        % Load Parameters
        speed     = emGet(em, 'drift speed', 'cones/sample', params);
        speedSD   = emGet(em, 'drift speed SD', 'cones/sample', params);
        
        % Make drift magnitude independent of sample time.
        correctionForSampleTime = 1000*sampTime;
    
        % Generate random move at each sample time
        theta = 360 * randn + 0.1 * correctionForSampleTime*(1 : nFrames)';
        direction = [cosd(theta) sind(theta)];
        s = speed + speedSD * randn(nFrames, 1);
        pos = filter(1,[1 -1],bsxfun(@times, direction, s)) + pos;
    end
    
    if (microSaccadesOnlyFlag)
        % Save eye movements due to drift and tremor
        driftTremorEMpos = pos;
    end
    
    % generate eye movement for micro-saccade
    if emFlag(3)
        % Load microsaccade model parameters
        interval = emGet(em, 'msaccade interval');
        intervalSD = emGet(em, 'msaccade interval SD');
        dirSD = emGet(em, 'msaccade dir SD', 'deg');
        speed = emGet(em, 'msaccade speed', 'cones/sample', params);
        speedSD = emGet(em, 'msaccade speed SD', 'cones/sample', params);
        
        % For the logic of this routine to be correct, we are assuming
        % that the integration time (sampTime) is one second or less.
        % Check that condition here.
        if (sampTime > 1)
            error('Sample time (aka integration time is too long (greater than 1 sec)');
        end
        
        % Compute microsaccade occurence times.
        %
        % The use of randn means that we might end up with some
        % negative intervals, we just get rid of these.

        minInterval = 0.010;
        saccadeIntervals = interval + randn(nFrames, 1) * intervalSD;
        saccadeIntervals(saccadeIntervals < 0) = minInterval;
        % Old way of getting rid of negative intervals, makes no sense to
        % us
        % t(t < 0.3) = 0.3 + 0.1*rand;     % get rid of negative times
        
        % Convert saccade intervals to saccade times
        saccadeTimes = cumsum(saccadeIntervals);                  
        
        % Convert to integer sample indices when we'll make a saccade
        saccadeTimeIndices = round(saccadeTimes / sampTime);
        
        % We now truncate the sample time indices to the duration that
        % we are actually simulating.
        saccadeTimeIndices = saccadeTimeIndices(saccadeTimeIndices <= nFrames); 
        %saccadeTimeIndices = tPos(1:find(saccadeTimeIndices <= nFrames, 1, 'last'));
        
        % Make sure we do at least one saccade, and in this case put it in
        % the middle of the time sequence
        %tPos = max(tPos,1);   
        if (isempty(saccadeTimeIndices))
            saccadeTimeIndices = round(nFrames/2);
        end
        
        % Comptue saccades by looping over when they happen
        for ii = 1:length(saccadeTimeIndices)
            % Find out where we are at the time of a saccade
            curPos = pos(saccadeTimeIndices(ii), :);
            
            % Saccade duration and direction are computed from the current
            % position.  The code computes a unit vector headed generally
            % towards the origin, with some angular noise.  This code makes
            % our algorithm have the feature that saccade amplitude and
            % duration is determined by the average distance from fixation
            % that drift and tremor produce in a inter-saccade interval.
            % This might be good, but we can't have this feature and
            % control saccade duration and magnitude independently of drift
            % and tremor. Indeed, without drift or tremor, there will be no
            % microsaccades.
            duration = round(sqrt(curPos(1)^2 + curPos(2)^2)/speed);
            saccadeAngle = atand(curPos(2)/curPos(1)) + dirSD * randn;
            saccadeDirection = [cosd(saccadeAngle) sind(saccadeAngle)];
            saccadeDirection = abs(saccadeDirection) .* (2*(curPos < 0) - 1);
            
            % Compute the offset we need to add to the positions vector, to
            % put in this saccade.  Start with an array of zeros, and
            % compute position offsets over the duration of the saccade.
            offset = zeros(nFrames, 2);
            indx = saccadeTimeIndices(ii):min(saccadeTimeIndices(ii) + duration - 1, nFrames);
            
            % Each saccade has a slightly different speed, compute the
            % speed for this saccade.
            cursaccadeSpeed = speed + speedSD * randn;
            if cursaccadeSpeed < 0, cursaccadeSpeed = speed; end
            
            % For each frame of the saccade, add in the offset for that
            % frame.
            offset(indx, 1) = cursaccadeSpeed*saccadeDirection(1);
            offset(indx, 2) = cursaccadeSpeed*saccadeDirection(2);
            
            % Convert offsets into cumulative saccade position delta, and
            % add into the positions we already had.
            pos = pos + cumsum(offset);
        end
    end
    
    % Remove eye movements due to drift and tremor
    if (microSaccadesOnlyFlag)
        pos = pos - driftTremorEMpos;
    end
    
    nTrialsPos(nn,:,:) = pos;
    pos = zeros(nFrames, 2);
end

%% Adjustments for the return

% Round to discrete cone steps
nTrialsPos = round(nTrialsPos);

% Return in microns as well
nTrialsPosMicrons = nTrialsPos * obj.patternSampleSize(1)*1e6;

% The positions in the mosaic is always just 1 sequence of eye movements
obj.emPositions = squeeze(nTrialsPos(1,:,:)); 

end