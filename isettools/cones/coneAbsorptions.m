function sensor = coneAbsorptions(sensor, oi, varargin)
% Compute a temporal array of cone absorptions including eye movements
%
%   sensor = coneAbsorptions(sensor, oi);
%
% This function is an extended form of sensorCompute, which was designed
% for a single capture.  This one loops through a series of eye fixations
% and calculates the voltage image at each position, concatenating the
% voltage images to (row,col,total nSamp).  The calculation for each
% position uses sensorCompute.
%
% The computational strategy is for each fixation, first calculate the mean
% (noiseless) voltage and then add frame-by-frame noise. The reason for
% this strategy is that sometimes the eye positions are the same.  By only
% using sensorCompute() for the different positions, we save computational
% cost.
%
% To simulate eye (sensor) movements, we expand the sensor using
% sensorHumanResize to be large enough to account for all possible eye
% positions. We compute with the large image and then remove unneeded
% rows/columns from the output voltage image.
%
% Example: An eye movement example
%
%   scene = sceneCreate;
%   oi = oiCreate;
%   oi = oiCompute(scene,oi);
%   sensor = sensorCreate('human');
%   positions = randi(5,20,2);    % Units are number of cones
%   sensor = sensorSet(sensor,'sensor positions',positions);
%   sensor = coneAbsorptions(sensor, oi);
%   v = sensorGet(sensor,'volts');
%   vcNewGraphWin; implay(v/max(v(:)))
%
% See also:  s_rgcScene2Cones, sensorHumanResize, sensorComputeSamples
%
% (c) Stanford VISTA Team
% 
% 5/28/15 xd, dhb  make this routine respect sensor's noise flag


%% check inputs
if notDefined('sensor'),  error('Need sensor'); end
if notDefined('oi'),      error('Need optical image'); end

% Get and verify eye movement parameters
xpos = sensorGet(sensor,'positions x');
ypos = sensorGet(sensor,'positions y');

%% Calculate volts.


% loop across positions, suppress waitbar in sensorCompute
% wbarState = ieSessionGet('wait bar'); ieSessionSet('wait bar','off');
nPos = length(xpos);

% Pad to the sensor to max size
rows = [-min([ypos(:); 0]) max([ypos(:); 0])];
cols = [max([xpos(:); 0]) -min([xpos(:); 0])];
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

% The LMS has a large, full-stack representation of the cones.  We know at
% each moment in time where the eye is.  This next piece of code takes the
% LMS stick accounts for eye position, and places the relevant LMS data
% into the volts for the original size sensor.
sz = sensorGet(sensor, 'size');
volts = zeros([sz nPos]);   % 3D - Spatial size by number of eye positions
for p = 1:nPos
    
    % cropping
    tmp = LMS(1+rows(1)+ypos(p):end-rows(2)+ypos(p), ...
        1+cols(1)-xpos(p):end-cols(2)-xpos(p),:);
    
    % select cone type to match the mosaic
    volts(:,:,p) = sum(tmp .* msk, 3);
end

% Set the volts field
% Add photon noise - though perhaps we should be adding general noise using
% sensorAddNoise.
%
sensor = sensorSet(sensor, 'volts', volts);

% Add photon noise
% This does not distinguish between noise flag value 1 and 2, which for
% camera sensors determines whether it's just shot noise or shot noise plus
% electronics noise.
if (sensorGet(sensor, 'noise flag'))
    sensor = sensorSet(sensor, 'volts', noiseShot(sensor));
end

% Someday, consider running Rieke adapt to get the current and then add
% noise based on Rieke's data (riekeAddNoise))


end