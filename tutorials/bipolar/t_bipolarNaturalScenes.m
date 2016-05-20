% t_bipolar
% 
% Implement full retinal pathway with sequential computations of the cone,
% bipolar and RGC responses.
% 
% This is the first implementation of the bioplar object.
% 
% 5/2016 JRG (c) isetbio team

%% Initialize
clear
ieInit

%% Load image sequence
tic
nSteps = 3600;%;

% Natural scene with eye movements stimulus
rdt = RdtClient('isetbio');
rdt.crp('resources/data/rgc');
data = rdt.readArtifact('testmovie_schemeA_8pix_Identity_8pix', 'type', 'mat');
testmovie = data.testmovie.matrix;

% load('/Users/james/Documents/MATLAB/akheitman/NSEM_mapPRJ/Stimuli/NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat')

paramsMovie.nSteps = nSteps;
paramsMovie.timeInterval = .00025;
paramsMovie.expTime = .00025;
paramsMovie.fov = 1;
iStim = ieStimulusMovie(testmovie.matrix(:,:,0+[1:3600]),paramsMovie);

% load('/Users/james/Documents/MATLAB/isetbio misc/iStim_NSEM_800fr.mat')
% load('/Users/james/Documents/MATLAB/isetbio misc/iStim_NSEM_3600fr.mat')
% load('/Users/james/Documents/MATLAB/isetbio misc/iStim_NSEM_9600fr.mat')
ieMovie(iStim.sceneRGB(:,:,1:nSteps))
absorptions = iStim.sensor; % cone isomerizations

figure; plot(squeeze(iStim.sceneRGB(4,10,:,1)))
% figure; plot(reshape(iStim.sceneRGB(:,:,:,1),size(iStim.sceneRGB,1)*size(iStim.sceneRGB,2),size(iStim.sceneRGB,3))')
xlabel('Time (msec)','fontsize',14); ylabel('Stimulus Intensity','fontsize',14)

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

% Set circular convolution, only steady state
paramsOSL.convolutionType = 1; 

% Compute the outer segment response to the absorptions with the linear
% model.
osL = osCompute(osL,absorptions,paramsOSL);

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

absorptionsVolts = sensorGet(absorptions,'volts');
paramsOS.bgVolts = 10*mean(absorptionsVolts(:));
osBp = osCompute(osBp,absorptions,paramsOS);

% % Plot the photocurrent for a pixel.
osPlot(osBp,absorptions);
%% Find bipolar responses
% The bipolar object takes as input the outer segment current. Bipolar
% processing consists of a spatial convolution and a temporal
% differentiator that matches the IR of the measured RGC output.

% Choose which outer segment to use for input.
% os = osL;
os = osBp;

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
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaBpSu = irCreate(bp, params);

% Create a subunit model for the on midget ganglion cell parameters
innerRetinaBpSu.mosaicCreate('model','Subunit','type','off parasol');
innerRetinaBpSu.mosaic{1}.mosaicSet('numberTrials',1);
% % Uncomment to get rid of spatial nonlinearity
newRectifyFunction = @(x) x;
innerRetinaBpSu.mosaic{1}.mosaicSet('rectifyFunction',newRectifyFunction);

% irPlot(innerRetinaBpSu,'mosaic');

% Compute RGC mosaic responses
innerRetinaBpSu = irCompute(innerRetinaBpSu, bp);
% irPlot(innerRetinaBpSu, 'psth');
% irPlot(innerRetinaBpSu, 'raster','cell',[2 2]);

irPlot(innerRetinaBpSu, 'raster','cell',[5 1]);
irPlot(innerRetinaBpSu, 'psth','cell',[5 1]);


toc