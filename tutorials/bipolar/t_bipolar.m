% t_bipolar
% 
% Implement full retinal pathway with sequential computations of the cone,
% bipolar and RGC responses. 
% 
% The stimulus is a contrast reversing grating with high spatial frequency.
% This implements the classic Hochstein & Shapley (1976) experiment which
% reveals the nonlinear spatial responses of RGCs. If RGCs had purely
% linear spatial processing, they would not fire in response to
% high-frequency gratings, because the linear response over the receptive
% field would average out to zero. However, if there were rectifying
% spatial subunits with the receptive field, this would lead to RGC firing,
% as the Hochstein & Shapley paper definitively demonstrated.
% 
% A contrast-reversing grating is either generated or loaded from the RDT.
% The scene, the oi, the sensor and the cone outer segment responses are
% computed. The cone outer segment responses are passed to the isetbio
% bipolar object, which can have nonlinear spatial subunits and has a
% temporal response that reproduces the appropriate RGC IR downstream. The
% bipolar responses are fed into the RGC object, which only carries out
% spatial processing at this point.
% 
% This is the first implementation of the bioplar object.
% 
% 5/2016 JRG (c) isetbio team

%% Initialize
clear;
ieInit;

%% Load contrast-reversing grating movie
% Grating subunit stimulus

clear params
stimP.fov      = 1; % degrees
stimP.barWidth = 2;
stimP.nSteps   = 400;
stimP.expTime = 0.001; % sec
stimP.timeInterval = 0.001; % sec
% % Select stimulus generation method
% % Uncomment here to change parameters. Warning: slow!
% iStim = ieStimulusGratingSubunit(stimP);

% % Uncomment here to grab a precomputed scene and sensor for the grating
% % subunit from the RDT.
rdt = RdtClient('isetbio');
% rdt.credentialsDialog();
rdt.crp('/resources/data/istim');
data = rdt.readArtifact('iStim_subunitGrating', 'type', 'mat');
iStim = data.iStim;

absorptions = iStim.absorptions; % cone isomerizations

% Show movie. Hold mouse over a bar to check that contrast reverses.
ieMovie(iStim.sceneRGB);

%% Plot one pixel from stimulus movie
vcNewGraphWin([],'upperleftbig'); 
plot(squeeze(iStim.sceneRGB(40,40,:,1)))
xlabel('Time (msec)','fontsize',14); 
ylabel('Stimulus level','fontsize',14)
title('Stimulus level for reversing grating at 5 Hz');
set(gca,'fontsize',16); grid on

%% Outer segment calculation - linear model
% The iStim structure generates the movie, the scene, the oi and the
% cone absorptions. The next step is to get the outer segment current. The
% linear outer segment model is employed here.

% Initialize
osL = osCreate('linear');

% Set size of retinal patch based on absorptions sensor object
patchSize = sensorGet(absorptions,'width','m');
osL = osSet(osL, 'patch size', patchSize);

% Set time step of simulation equal to absorptions sensor object
timeStep = sensorGet(absorptions,'time interval','sec');
osL = osSet(osL, 'time step', timeStep);

% Compute the outer segment response to the absorptions with the linear
% model.
osL = osCompute(osL,absorptions);

% % Plot the photocurrent for a pixel.
osPlot(osL,absorptions);

%% Outer segment calculation - biophysical model
% As an alternative to the linear model, we can employ Fred's biophysical
% model that captures properties of temporal adaptation.
osBp = osCreate('BioPhys');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','m');
osBp = osSet(osBp, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
osBp = osSet(osBp, 'time step', timeStep);

% Compute the outer segment response to the absorptions with the
% biophysical model.
osBp = osCompute(osBp,absorptions);

% % Plot the photocurrent for a pixel.
osPlot(osBp,absorptions);

%% Find bipolar responses
% The bipolar object takes as input the outer segment current. Bipolar
% processing consists of a spatial convolution and a temporal
% differentiator that matches the IR of the measured RGC output.

% Choose which outer segment to use for input.
os = osL;
% os = osBp;

% Create the bipolar object
bp = bipolar(os);

% Compute the bipolar response
bp = bipolarCompute(bp, os);

% Plot the response of every bipolar cell.
bipolarPlot(bp,'response');

%% Find RGC responses
% Build and IR object that takes as input the bipolar mosaic.

% Initialize.
clear params innerRetinaBpSu
params.name      = 'Bipolar with nonlinear subunits'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 9;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaBpSu = irCreate(bp, params);

% Create a subunit model for the on midget ganglion cell parameters
innerRetinaBpSu.mosaicCreate('model','Subunit','type','on midget');

% % Uncomment to get rid of spatial nonlinearity
% newRectifyFunction = @(x) x;
% innerRetinaBpSu.mosaic{1}.mosaicSet('rectifyFunction',newRectifyFunction);

% irPlot(innerRetinaBpSu,'mosaic');

% Compute RGC mosaic responses
innerRetinaBpSu = irCompute(innerRetinaBpSu, bp);
irPlot(innerRetinaBpSu, 'psth');