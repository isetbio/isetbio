% funcftion [scene, sceneRGB, oi, sensor] = movieCreate(varargin)


% (HJ) ISETBIO TEAM, 2014
% (JRG) modified 10/2015

%% Init
clear
ieInit;

%% Compute scene, oi, sensor over time

% Set up Gabor stimulus using sceneCreate('harmonic',params)
fov = 0.6;
params.freq = 6; params.contrast = 1;
params.ph  = 0;  params.ang = 0;
params.row = 64; params.col = 64;
params.GaborFlag = 0.2; % standard deviation of the Gaussian window
params.nSteps = 20;
params.expTime = 0.01;
params.timeInterval = 0.01;

% Set up scene, oi and sensor
scene = sceneCreate('harmonic', params);
scene = sceneSet(scene, 'h fov', fov);
oi  = oiCreate('wvf human');
sensor = sensorCreate('human');
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);
sensor = sensorSet(sensor, 'exp time', params.expTime); 
sensor = sensorSet(sensor, 'time interval', params.timeInterval); 

% Compute sceneRGB and cone absorptions 
sceneRGB = zeros([sceneGet(scene, 'size') params.nSteps 3]); % 3 is for R, G, B
sensorPhotons = zeros([sensorGet(sensor, 'size') params.nSteps]);
stimulus = zeros(1, params.nSteps);
fprintf('Computing cone isomerization:    ');

% Loop through frames to build movie
for t = 1 : params.nSteps
    fprintf('\b\b\b%02d%%', round(100*t/params.nSteps));
    
    params.ph = 2*pi*(t-1)/params.nSteps; % one period over nSteps
    scene = sceneCreate('harmonic', params);
    scene = sceneSet(scene, 'h fov', fov);
    
    % Get scene RGB data    
    sceneRGB(:,:,t,:) = sceneGet(scene,'rgb');    
    % Compute optical image
    oi = oiCompute(scene, oi);    
    % Compute absorptions
    sensor = sensorCompute(sensor, oi);
    sensorPhotons(:,:,t) = sensorGet(sensor, 'photons');
    
    vcAddObject(scene); sceneWindow
end
fprintf('\n');

% Set the stimuls into the sensor object
sensor = sensorSet(sensor, 'photons', sensorPhotons);

%% Build outer segment

% Input = cone current
os = osCreate('linear');
os = osCompute(os, sensor);

% % Input = RGB
% os = osCreate('identity');
% os = osSet(os, 'rgbData', sceneRGB);

%% Build rgc

rgc1 = rgcCreate('linear', scene, sensor, os, 'right', 3.0, 180);

rgc1 = rgcCompute(rgc1, os);

rgcPlot(rgc1, 'mosaic');
rgcPlot(rgc1, 'linearResponse');
rgcPlot(rgc1, 'spikeResponse');
%% Build rgc response movie
 
% rgcMovie(rgc1, os);