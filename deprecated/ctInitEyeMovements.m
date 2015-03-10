function [sensor,ic] = ctInitEyeMovements(sensor, scene, oi, nSamples, randSeed)
% THIS FUNCTION IS NOT IN USE ANYMORE. PLEASE USE EMINIT INSTEAD
% Init eye movement parameters in the sensor structure
%
%   sensor = ctInitEyeMovements(sensor, scene, oiD, nSamples, eyeMoveType, randSeed)
%
% General Process:
%   1. Check eyeMoveType and set random seed
%   2. Generate moving postion x and y
%   3. Generate frames per position according to distribution and nSamples
%   4. Generate linear eye movement for testing
%
% Input Parameters:
%   sensor       - eye sensor to be used, if empty, program will create one
%   scene        - ISET scene structure
%   oiD          - Optical image data
%   nSamples     - number of samples to be generated
%   randSeed     - seed to be used in random generator
%
% Output Parameter:
%   sensor       - sensor with eye movement related parameters set
%
% Example:
%    sensor = ctInitEyeMovements(sensor, scene, oi, 100, 1)
%
% (HJ) Copyright PDCSOFT TEAM 2013

%% Check inputs and Init
warning('This function is not in use. See emCreate instead');
if isempty(sensor)
    sensor = sensorCreate('human');
    sensor = sensorSet(sensor,'exp time',0.05);
end
if nargin < 4, nSamples = 1000; end

% Set rand seed
if nargin < 5,  rng('shuffle'); % Random according to time
else            rng(randSeed); % Set random seed
end

% Init parameters for Gaussian RV
% This has to be scaled for the sensor exposure time.  We think that these
% are the sd for 25 ms.  So the std dev for 1 ms is [0.02, 0.03]/sqrt(25)
% In general, we should get the sensor exposure time and make the std dev
% [0.02, 0.03] * sqrt(expTime)/sqrt(25)
% (BW)
center = [0 0];
%sigmaX = 0.02 / sqrt(5);
%sigmaY = 0.03 / sqrt(5); % Could be changed to correlation matrix later
sigmaX = 0;
sigmaY = 0;
% (BW)
%s = sqrt(sensorGet(sensor,'exp time','ms')/25);
%sigmaX = sigmaX*s;
%sigmaY = sigmaY*s;

% Compute sensor fov and size
fov = sensorGet(sensor, 'fov', scene, oi); 
sz = sensorGet(sensor, 'size');

%% Generate Eye
% Staring Mode - eye is fixed at a position but moves randomly
% around it, with some microsaccade distribution.
%
% Set eye move to be in range [50% 50%] with Guassian distribution
x    = randn(nSamples,1)*sigmaX + center(1);
y    = randn(nSamples,1)*sigmaY + center(2);
% Round to sz
% x    = x * fov;  y = y * fov;
xPos = round(x*sz(1)/fov);
yPos = round(y*sz(2)/fov);

% Group the same positions
[pos,~,ic]  = unique([xPos yPos],'rows');
x = pos(:,1)/2/sz(1)*fov;
y = pos(:,2)/2/sz(2)*fov;
% Compute frame per position
f    = hist(ic,unique(ic)); % frames per position
f(1) = f(1) + nSamples - sum(f); % make sure sum(f) == nSamples
% Set parameters to sensor
sensor = sensorSet(sensor,'movement positions',[x y]);
sensor = sensorSet(sensor,'frames per position',f);
% Add vcObject
vcAddAndSelectObject('scene',scene);
vcAddAndSelectObject('oi',oi);

end