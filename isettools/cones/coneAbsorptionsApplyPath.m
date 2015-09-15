function sensor = coneAbsorptionsApplyPath(sensor, LMS, msk, varargin)
% sensor = coneAbsorptionsApplyPath(sensor, LMS, msk)
%   Detailed explanation goes here

% Get and verify eye movement parameters
xpos = sensorGet(sensor,'positions x');
ypos = sensorGet(sensor,'positions y');
nPos = length(xpos);

% Pad to the sensor to max size
rows = [-min([ypos(:); 0]) max([ypos(:); 0])];
cols = [max([xpos(:); 0]) -min([xpos(:); 0])];
rows = [max(rows) max(rows)];
cols = [max(cols) max(cols)];

ip = inputParser;
ip.addOptional('rows', rows);
ip.addOptional('cols', cols);
ip.parse(varargin{:});

rows = ip.Results.rows;
cols = ip.Results.cols;

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

end

