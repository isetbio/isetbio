function sensor = emGenSequence(sensor)
%% Generate eye movement sensor sequence
%    sensor = emGenSequence([sensor]);
%
%  Inputs:
%    sensor - human sensor structure, see sensorCreate('human') for detail
%             Here, we require that all the eye-movement related parameters
%             set. For details about the parameters, see eyemovementInit
%             for detail
%
%  Outputs:
%    sensor - human sensor structure with eye movement position sequence
%             set. Position sequence could be retrievedy be
%             sensorGet(sensor, 'sensor positions');
%
%  Example:
%    sensor = sensorCreate('human');
%    sensor = eyemovementInit(sensor);
%    sensor = emGenSequence(sensor);
%
%  See also:
%    eyemovementInit
%
%  (HJ) Copyright PDCSOFT TEAM 2014

%% Init
if notDefined('sensor'), error('human sensor required'); end
rSeed = sensorGet(sensor, 'human rseed');
if isempty(rSeed)
    rSeed = rng;
    sensor = sensorSet(sensor, 'human rseed', rSeed);
else
    rng(rSeed);
end
emFlag = sensorGet(sensor, 'em flag');
if isempty(emFlag), error('eye movement type not defined'); end

% Init positions
pos = sensorGet(sensor, 'sensorpositions');
if isempty(pos), error('sensor positions length unknown'); end
pos = zeros(size(pos));

% Load general parameters
em        = sensorGet(sensor, 'eye movement');
sampTime  = emGet(em, 'sample time');
seqLen    = size(pos, 1);
% mperdeg   = vcConstants('mmperdeg') / 1000;
% coneWidth = pixelGet(sensorGet(sensor, 'pixel'), 'width');

%% Generate eye movement for tremor
if emFlag(1)
    % Load parameters
    amplitude  = emGet(em, 'tremor amplitude', 'cones/sample');
    interval   = emGet(em, 'tremor interval');
    intervalSD = emGet(em, 'tremor interval SD');
    
    % Compute time of tremor occurs
    t = interval + randn(seqLen, 1) * intervalSD;
    t(t < 0.001) = 0.001; % get rid of negative values
    t = cumsum(t);
    tPos = round(t / sampTime);
    tPos = tPos(1:find(tPos <= seqLen, 1, 'last'));
    
    % Generate random move on the selected time
    direction = rand(length(tPos),1);
    pos(tPos, :) = amplitude * [direction sqrt(1-direction.^2)];
    pos = pos .* (2*(randn(size(pos))>0)-1); % shuffle the sign
    pos = cumsum(pos, 1);
end

%% Generate eye movement for drift
if emFlag(2)
    % Load Parameters
    speed     = emGet(em, 'drift speed', 'cones/sample');
    speedSD   = emGet(em, 'drift speed SD', 'cones/sample');
    
    % Generate random move at each sample time
    theta = 360 * randn + 0.1 * (1 : seqLen)';
    direction = [cosd(theta) sind(theta)];
    s = speed + speedSD * randn(seqLen, 1);
    pos = filter(1, [1 -1], bsxfun(@times, direction, s)) + pos;
end

%% Generate eye movement for micro-saccade
if emFlag(3)
    % Load parameters
    interval   = emGet(em, 'msaccade interval');
    intervalSD = emGet(em, 'msaccade interval SD');
    dirSD      = emGet(em, 'msaccade dir SD', 'deg');
    speed      = emGet(em, 'msaccade speed', 'cones/sample');
    speedSD    = emGet(em, 'msaccade speed SD', 'cones/sample');
    
    % compute time of occurance
    t = interval + randn(seqLen, 1) * intervalSD;
    t(t < 0.3) = 0.3 + 0.1*rand; % get rid of negative values
    t = cumsum(t);
    tPos = round(t / sampTime);
    tPos = tPos(1:find(tPos <= seqLen, 1, 'last'));
    
    % Compute positions
    for ii = 1 : length(tPos)
        curPos = pos(tPos(ii), :);
        duration = round(sqrt(curPos(1)^2 + curPos(2)^2)/speed);
        direction = atand(curPos(2) / curPos(1)) + dirSD * randn;
        direction = [cosd(direction) sind(direction)];
        direction = abs(direction) .* (2*(curPos < 0) - 1);
        
        offset = zeros(seqLen, 2);
        indx = tPos(ii):min(tPos(ii) + duration - 1, seqLen);
        curSpeed = speed + speedSD * randn;
        if curSpeed < 0, curSpeed = speed; end
        offset(indx, 1) = curSpeed*direction(1);
        offset(indx, 2) = curSpeed*direction(2);
        
        pos = pos + cumsum(offset);
    end
end

%% Set sensor position back to sensor
%  pos = round(cumsum(pos, 1));
pos = round(pos);
sensor = sensorSet(sensor, 'sensorpositions', pos);
end