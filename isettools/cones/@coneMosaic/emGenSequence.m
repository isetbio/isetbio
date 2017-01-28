function nTrialsPos = emGenSequence(obj, nFrames, varargin)
%EMGENSEQUENCE  Generate sequence of eye movements
%
%   nTrialsPos = EMGENSEQUENCE(obj,nFrames,'nTrials',1,'em',emCreate);
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
%     nTrialsPos - nTrials x nFrames x 2 matrix of eye positions in units of cone positions
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

for nn=1:nTrials
    %% generate eye movement for tremor
    if emFlag(1)
        % Load parameters
        amplitude  = emGet(em, 'tremor amplitude', 'cones/sample');
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
    
    % generate eye movement for drift
    if emFlag(2)
        % Load Parameters
        speed     = emGet(em, 'drift speed', 'cones/sample');
        speedSD   = emGet(em, 'drift speed SD', 'cones/sample');
        
        % Generate random move at each sample time
        theta = 360 * randn + 0.1 * (1 : nFrames)';
        direction = [cosd(theta) sind(theta)];
        s = speed + speedSD * randn(nFrames, 1);
        pos = filter(1,[1 -1],bsxfun(@times, direction, s)) + pos;
    end
    
    % generate eye movement for micro-saccade
    if emFlag(3)
        % Load microsaccade model parameters
        interval = emGet(em, 'msaccade interval');
        intervalSD = emGet(em, 'msaccade interval SD');
        dirSD = emGet(em, 'msaccade dir SD', 'deg');
        speed = emGet(em, 'msaccade speed', 'cones/sample');
        speedSD = emGet(em, 'msaccade speed SD', 'cones/sample');
        
        % Compute microsaccade occurence times
        t = interval + randn(nFrames, 1) * intervalSD;
        t(t < 0.3) = 0.3 + 0.1*rand;     % get rid of negative times
        t = cumsum(t);                   % Add them up
        
        % Convert to integer locations of the positions
        tPos = round(t / sampTime);
        
        % Finds the last nonzero element in the tPos array with a value less
        % than nFrames.  But that element has to be at least 1. BW doesn't
        % really understand the logic of the model here.
        tPos = tPos(1:find(tPos <= nFrames, 1, 'last'));
        tPos = max(tPos,1);   % HJ to check.
        
        % Compute positions
        for ii = 1 : length(tPos)
            curPos = pos(tPos(ii), :);
            duration = round(sqrt(curPos(1)^2 + curPos(2)^2)/speed);
            direction = atand(curPos(2)/curPos(1)) + dirSD * randn;
            direction = [cosd(direction) sind(direction)];
            direction = abs(direction) .* (2*(curPos < 0) - 1);
            
            offset = zeros(nFrames, 2);
            indx = tPos(ii):min(tPos(ii) + duration - 1, nFrames);
            curSpeed = speed + speedSD * randn;
            if curSpeed < 0, curSpeed = speed; end
            offset(indx, 1) = curSpeed*direction(1);
            offset(indx, 2) = curSpeed*direction(2);
            
            pos = pos + cumsum(offset);
        end
    end
    nTrialsPos(nn,:,:) = pos;
    pos = zeros(nFrames, 2);
end

%% Adjustments for the return

% If a hex mosaic, we  amplify the positions by the resampling factor.
if (isa(obj, 'coneMosaicHex'))
    nTrialsPos = nTrialsPos * obj.resamplingFactor;
end

% Round to discrete cone steps
nTrialsPos = round(nTrialsPos);

% The positions in the mosaic is always just 1 sequence of eye movements
obj.emPositions = squeeze(nTrialsPos(1,:,:)); 

end