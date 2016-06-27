function [LMS, msk] = coneAbsorptionsLMS(sensor, oi)
% [LMS, msk] = coneAbsorptionsLMS(sensor, oi)
% 
% This function performs the first half of the code in coneAbsorptions.
% The sensor is resized to a larger image based on the bounds of the eye
% movement path.  This larger sensor is used to calculate LMS absorptions
% which will cropped accordingt to the eyemovement path in
% coneAbsorptionsApplyPath
%
% 7/27/15  xd  moved code from coneAbsorptions

%% check inputs
if notDefined('sensor'),  error('Need sensor'); end
if notDefined('oi'),      error('Need optical image'); end

% Get and verify eye movement parameters
xpos = sensorGet(sensor,'positions x');
ypos = sensorGet(sensor,'positions y');

%% Calculate volts.
% Pad to the sensor to max size
rows = [-min([ypos(:); 0]) max([ypos(:); 0])];
cols = [max([xpos(:); 0]) -min([xpos(:); 0])];

% need to pad the same to both sides so that the sensor center is not
% changed
rows = [max(rows) max(rows)];
cols = [max(cols) max(cols)];

coneType = sensorGet(sensor, 'cone type');
sensor2  = sensorHumanResize(sensor, rows, cols);

% With the noise flag set to 0, we compute only the mean.  No photon noise
% or sensor noise.

% Compute a full mosaic of L, M, and S absorptions
sz = sensorGet(sensor2, 'size');
LMS = zeros([sz 3]);
msk = zeros([size(coneType) 3]);
for ii = 2 : 4 % L, M, S
    pattern = ii * ones(sz);
    
    % There are two uses of pattern.  They could be in conflict.  We should
    % eliminate the need for two. Probably we need to preserve 'pattern'
    % because it is long-standing, and we need to deal with cone type in
    % some other way.
    sensor2 = sensorSet(sensor2,'pattern', pattern);
    sensor2 = sensorSet(sensor2,'cone type', pattern);
    sensor2 = sensorSet(sensor2,'noise flag',0);

    sensor2 = sensorCompute(sensor2, oi);
    LMS(:,:,ii-1) = sensorGet(sensor2, 'volts');
    
    % Store the positions of the cones of this type.  We use this later
    % when we fill in the voltages in the original sensor.
    msk(:,:,ii-1) = double(coneType == ii);
end

end

