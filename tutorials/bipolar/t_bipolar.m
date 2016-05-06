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

stimulusSelect = 1;
nSteps = 400;
barWidth = 2;

switch stimulusSelect
    case 1
        % Grating subunit stimulus
        clear params
        stimP.fov      = 1; % degrees
        stimP.barWidth = barWidth;
        stimP.nSteps   = nSteps;
        stimP.expTime = 0.001; % sec
        stimP.timeInterval = 0.001; % sec
        iStim = ieStimulusGratingSubunit(stimP);
        absorptions = iStim.absorptions; % cone isomerizations
        ieMovie(iStim.sceneRGB);
    case 2
        % Natural scene with eye movements stimulus
        rdt = RdtClient('isetbio');
        rdt.crp('resources/data/rgc');
        data = rdt.readArtifact('testmovie_schemeA_8pix_Identity_8pix', 'type', 'mat');
        testmovie = data.testmovie;
        absorptions = ieStimulusMovie(testmovie(:,:,1:nSteps)); 
end

%% Outer segment calculation - linear model

osL = osCreate('linear');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','um');
osL = osSet(osL, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
osL = osSet(osL, 'time step', timeStep);

% Set osI data to raw pixel intensities of stimulus
% osI = osSet(osI, 'rgbData', 0.6-iStim.sceneRGB);
osL = osCompute(osL,absorptions);

% % Plot the photocurrent for a pixel
osPlot(osL,absorptions);

%% Outer segment calculation - biophysical model

% osBp = osCreate('BioPhys');
% 
% % Set size of retinal patch
% patchSize = sensorGet(absorptions,'width','um');
% osBp = osSet(osBp, 'patch size', patchSize);
% 
% % Set time step of simulation equal to absorptions
% timeStep = sensorGet(absorptions,'time interval','sec');
% osBp = osSet(osBp, 'time step', timeStep);
% 
% % Set osI data to raw pixel intensities of stimulus
% % osI = osSet(osI, 'rgbData', 0.6-iStim.sceneRGB);
% osBp = osCompute(osBp,absorptions);
% 
% % % Plot the photocurrent for a pixel
% osPlot(osBp,absorptions);

%% Find bipolar responses

bp = bipolar(osL);

% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

bp = bipolarCompute(bp, osL);

bipolarPlot(bp);

%% Find RGC responses

clear params innerRetina0
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 7;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina0 = irCreate(bp, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetina0.mosaicCreate('model','Subunit','type','on midget');

irPlot(innerRetina0,'mosaic');
% Compute RGC mosaic responses

innerRetina0 = irCompute(innerRetina0, bp);
% irPlot(innerRetina0, 'psth');
irPlot(innerRetina0, 'linear');
% irPlot(innerRetina0, 'raster');

%% Bipolar input with linear subunits

clear params innerRetina0
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 7;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina0 = irCreate(bp, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetina0.mosaicCreate('model','LNP','type','on midget');

irPlot(innerRetina0,'mosaic');
% Compute RGC mosaic responses

innerRetina0 = irCompute(innerRetina0, bp);
% irPlot(innerRetina0, 'psth');
irPlot(innerRetina0, 'linear');
% irPlot(innerRetina0, 'raster');

%% Compute RGC without subunits
% Input = RGB
osI = osCreate('displayRGB');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','um');
osI = osSet(osI, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
osI = osSet(osI, 'time step', timeStep);

osI = osSet(osI, 'rgbData', 2*(iStim.sceneRGB-0.5));


clear params
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 7;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina1 = irCreate(osI, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetina1.mosaicCreate('model','lnp','type','on midget');
irPlot(innerRetina1,'mosaic');

innerRetina1 = irCompute(innerRetina1, osI);
irPlot(innerRetina1, 'linear');
% irPlot(innerRetina1, 'psth');

%% Compute RGC with subunits
% Input = RGB
osI = osCreate('displayRGB');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','um');
osI = osSet(osI, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
osI = osSet(osI, 'time step', timeStep);

osI = osSet(osI, 'rgbData', 2*(iStim.sceneRGB-0.5));


clear params
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 7;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina2 = irCreate(osI, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetina2.mosaicCreate('model','subunit','type','on midget');
% irPlot(innerRetina2,'mosaic');

innerRetina2 = irCompute(innerRetina2, osI);
irPlot(innerRetina2, 'linear');
% irPlot(innerRetina1, 'psth');
