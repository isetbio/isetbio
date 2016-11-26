function pos = emGenSequence(obj, nFrames, varargin)
% Generate eye movement path
%
% Inputs
%   obj     - cone mosaic object
%   nFrames - number of frames to generate
%
% Optional inputs (key value pairs in varargin):
%   'em'    - eye movement structure, see emCreate for details
%   'rSeed' - random seed to be used
%
% Ouputs
%   pos     - nFramesx2 matrix of eye positions in units of
%             cone positions
%
% See also:
%   emCreate

%% parse input
p = inputParser;
p.addRequired('nFrames', @isscalar);
p.addParameter('em', emCreate, @isstruct);
p.addParameter('rSeed', [], @isscalar);

% set parameters
sampTime  = obj.integrationTime;

p.parse(nFrames, varargin{:});
em = p.Results.em;
em = emSet(em, 'sample time', sampTime);

if ~isempty(p.Results.rSeed), rng(p.Results.rSeed); end
emFlag = emGet(em, 'em flag');
pos = zeros(nFrames, 2);

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

% If the mosaic is a hex mosaic, we need to amplify the positions by
% the resampling factor.
if (isa(obj, 'coneMosaicHex'))
    pos = pos * obj.resamplingFactor;
end

pos = round(pos);
obj.emPositions = pos;

end