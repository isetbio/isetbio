%% s_rieke_speedTest
%  This script tests the speed for rieke's cone adaptation model
%  For the detailed description of the model, please see coneAdapt,
%  riekeAdapt, riekeInit, etc.
%
%  (HJ) ISETBIO TEAM, 2014

%% Init parameters
s_initISET;
p = riekeInit;
deltaT  = 0.001; % sample every 1 ms
expTime = 0.05;  % integration time 50 ms

sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'exposure time', expTime);
sensor = sensorSet(sensor, 'time interval', deltaT);
sensor = sensorSet(sensor, 'size', [100 100]); % make it a square

%% Small test case
%  This should finish within a couple of seconds
duration = 1; % 1 second
fov = 1; % 1 deg

sensor = sensorSetSizeToFOV(sensor, fov);
sz = sensorGet(sensor, 'size');

% create random input
% the random isomerization rate is controlled within 1000 to 100000
pRate = 99000 * rand([sz duration/deltaT]) + 1000;
sensor = sensorSet(sensor, 'photon rate', pRate);

% compute adapted current
% this computation uses a default background (initial) isomerizaiton rate
% equals to the mean of the stimulus
tic;
[~, cur_small] = coneAdapt(sensor, 'rieke');
toc;

%% Large test cases
%  This could takes you 10 minutes if there is no 'not enough memory' error
% duration = 1; % 1 second
% fov = 10; % 10 x 10 deg, more than 3M cones
% 
% sensor = sensorSetSizeToFOV(sensor, fov);
% sz = sensorGet(sensor, 'size');
% 
% % create random input
% % the random isomerization rate is controlled within 1000 to 100000
% pRate = 99000 * rand([sz duration/deltaT]) + 1000;
% sensor = sensorSet(sensor, 'photon rate', pRate);
% 
% % compute adapted current
% % this computation uses a default background (initial) isomerizaiton rate
% % equals to the median of the stimulus
% tic;
% [~, cur_large] = coneAdapt(sensor, 'rieke');
% toc;
