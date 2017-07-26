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
% ieInit

%% Load image sequence
% Fred's code for examining impulse response
osModelType = 'linear';
% osModelType = 'biophys';

% Set up parameters for stimulus.
timeStep = (1/150);        % time step
nSamples = .4/timeStep;        % 2000 samples

flashIntens = 5000;    % flash intensity in R*/cone/sec (maintained for 1 bin only)

% Create human sensor.
sensor = sensorCreate('human');
% sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'size', [64 64]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', timeStep);


% Create stimulus.
stimulus = zeros(nSamples, 1);
if strcmp(osModelType,'linear')
    stimulus = zeros(64,64,nSamples);
%     stimulus(31,31,round(23*.001/timeStep)) = flashIntens;
    stimulus(31,31,round(23)) = flashIntens;
else
    stimulus = zeros(64,64,nSamples);
    stimulus(round(24*.001/timeStep)) = flashIntens;
end
% stimulus = reshape(stimulus, [1 1 nSamples]);

% Set photon rates. This is a kluge that appeared
% just for this test, and that should probably go
% away again. This is an artifact of directly specifying the stimulus
% in the sensor, and will not be an issue when the sensor
% is the result of a sensorCompute command on a scene and oi.
sensor = sensorSet(sensor, 'photon rate', stimulus);


% Create outersegment object and get the adapted response.
noiseFlag = 0;
if strcmp(osModelType,'linear')
    os = osLinear(sensor); 
    paramsOS.convolutionType = 1; 
else
    os = osBioPhys(sensor); 
    paramsOS.bgVolts = 10*mean(vectorize(sensorGet(sensor,'volts')));
end

os = osSet(os, 'noiseFlag', noiseFlag);
os = osCompute(os, sensor, paramsOS);

% Set size of retinal patch
patchSize = sensorGet(sensor,'width','m');
os = osSet(os, 'patch size', patchSize);

osPlot(os,sensor);

%%

%% Find bipolar responses


bpParams.filterType = 4;
bpParams.cellLocation = 58;
bp = bipolar(os,bpParams);
bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);
% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

bp = bipolarCompute(bp, os);

% bipolarPlot(bp,'response');

%% Build rgc

clear params
params.name      = 'Macaque inner retina 1'; % This instance

params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaSU = irCreate(bp, params);

% innerRetinaSU.mosaicCreate('model','glm','type','off midget');
innerRetinaSU.mosaicCreate('model','glm','type','on midget');

innerRetinaSU = irCompute(innerRetinaSU, bp);
irPlot(innerRetinaSU,'linear');


%% Measure the temporal impulse response for the original GLM model
% This response can be compared to the above response 
% Image parameters
params.image_size = 63;
params.meanLuminance = 1;
params.nsteps = 100;
params.fov = 0.8;

% Set impulse stimulus by only turning on a signal pixel at t = 1
sceneRGB = zeros(params.image_size, params.image_size, params.nsteps, 3);
sceneRGB(16,16,23,:) = 100*[10 10 10];
os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep', 1/150);

% Attach the movie to the object
os1 = osSet(os1, 'rgbData', sceneRGB);

%% Generate RGC object for simulated GLM prediction of response
% Set the parameters for the inner retina RGC mosaic. For the inner retina
% type irPhys, the values for eyeSide, eyeRadius and eyeAngle have no
% effect, because those are dependent on the properties of the retinal
% piece used in the Chichilnisky Lab experiment.
clear params
% Set parameters
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = 12; 
params.eyeAngle = 0; ntrials = 0;

% Determined at beginning to allow looping
params.experimentID = '2013-08-19-6'; % Experimental dataset
params.stimulusTest = 'WN'; % WN or NSEM
params.cellType = 'on parasol';         % ON or OFF Parasol

% Create object
innerRetina = irPhys(os1, params);
nTrials = 57; innerRetina = irSet(innerRetina,'numberTrials',nTrials);

innerRetina = irCompute(innerRetina,os1);
irPlot(innerRetina,'linear')

%%

% % % Run the t_bipolarImpulseResponse to get this ir
mpt = 40;
filter1 = squeeze(innerRetina.mosaic{1}.responseLinear(bpParams.cellLocation,1,23:end));
figure; plot((1/150)*(-0+[0:length(filter1)-1]),(filter1-filter1(mpt))./max(abs(filter1-filter1(mpt))),'linewidth',4)
filter2 = -squeeze(innerRetinaSU.mosaic{1}.responseLinear(3,3,23:end));
% figure;
hold on; 
plot((1/150)*(-0+[0:length(filter2)-1]),(filter2-filter2(mpt))./max(abs(filter2-filter2(mpt))),':r','linewidth',3);


grid on;

legend('Black Box','Cascade');
% title('Outer segment biophysical impulse response');
% title('RGC impulse response, averaged');
title('Validation of impulse response, On Parasol');
grid on;
axis([0 0.4 -1 1]);
xlabel('Time (sec)'); ylabel('Normalized Response');
set(gca,'fontsize',14);

axis([0 0.25 -.4 1]);