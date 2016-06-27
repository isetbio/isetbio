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

[LMS, msk] = coneAbsorptionsLMS(sensor, oi);
sensor = coneAbsorptionsApplyPath(sensor, LMS, msk);

end