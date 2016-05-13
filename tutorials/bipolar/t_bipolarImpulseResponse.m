% t_bipolar
% 
% Implement full retinal pathway with sequential computations of the cone,
% bipolar and RGC responses.
% 
% The stimulus is an impulse. The goal is to examine the impulse repsonse
% of the bipolar object with the differentiator 
% 
% 5/2016 JRG (c) isetbio team

%% Initialize
clear
ieInit

%% Load image sequence
% Fred's code for examining impulse response

% Set up parameters for stimulus.
nSamples = 4000;        % 2000 samples
timeStep = 1e-4;        % time step
flashIntens = 50000;    % flash intensity in R*/cone/sec (maintained for 1 bin only)

% Create human sensor.
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', timeStep);

% Create stimulus.
stimulus = zeros(nSamples, 1);
stimulus(1) = flashIntens;
stimulus = reshape(stimulus, [1 1 nSamples]);

% Set photon rates. This is a kluge that appeared
% just for this test, and that should probably go
% away again. This is an artifact of directly specifying the stimulus
% in the sensor, and will not be an issue when the sensor
% is the result of a sensorCompute command on a scene and oi.
sensor = sensorSet(sensor, 'photon rate', stimulus);

% Create outersegment object and get the adapted response.
noiseFlag = 0;
% os = osBioPhys();
os = osLinear()
os = osSet(os, 'noiseFlag', noiseFlag);
os = osCompute(os, sensor);

% Set size of retinal patch
patchSize = sensorGet(sensor,'width','m');
os = osSet(os, 'patch size', patchSize);

% % Set time step of simulation equal to absorptions
% timeStep = sensorGet(absorptions,'time interval','sec');
% os = osSet(os, 'time step', timeStep);
        
% figure; plot(squeeze(stimulus))
% xlabel('Time (msec)','fontsize',14); ylabel('Stimulus Intensity','fontsize',14)
osPlot(os,sensor);

%% Find bipolar responses

bp = bipolar(os);

% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

bp = bipolarCompute(bp, os);

% bipolarPlot(bp);

bpResponse = bipolarGet(bp,'responseCenter');
figure; plot(.1:.1:.1*length(bpResponse),-squeeze(bpResponse - bpResponse(end))./max(abs(squeeze(bpResponse - bpResponse(end)))));

load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/OnParasolExcFilters.mat')
% figure; 
hold on;
% plot(tme,max(abs(1))*mean(LinFilt)./max(abs(mean(LinFilt))));

meanLinFilt = mean(LinFilt);
% Downsample to appropriate time course
dsrate = length(tme)/max(tme);
rgcFilt = meanLinFilt(1:2*dsrate:end);
plot(rgcFilt);
xlabel('Time (msec)','fontsize',16);
ylabel('Normalized Amplitude','fontsize',16);
title('On Parsol Synaptic Input Impulse Response','fontsize',16);
legend('Bipolar IR','RGC Synaptic IR');
set(gca,'fontsize',16);
grid on;