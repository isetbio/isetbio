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
nSteps = 500;
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
%         iStim = ieStimulusGratingSubunit(stimP);

        load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/iStim_subunitGrating.mat');
%        load('/Users/james/Documents/MATLAB/isetbio misc/iStim_subunitGratingLong.mat')
 
        absorptions = iStim.absorptions; % cone isomerizations
%         ieMovie(iStim.sceneRGB);
    case 2
        % Natural scene with eye movements stimulus
        rdt = RdtClient('isetbio');
        rdt.crp('resources/data/rgc');
        data = rdt.readArtifact('testmovie_schemeA_8pix_Identity_8pix', 'type', 'mat');
        testmovie = data.testmovie;
        absorptions = ieStimulusMovie(testmovie(:,:,1:nSteps)); 
end

figure; plot(squeeze(iStim.sceneRGB(40,40,:,1)))
xlabel('Time (msec)','fontsize',14); ylabel('Stimulus Intensity','fontsize',14)

%% Outer segment calculation - linear model

osL = osCreate('linear');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','m');
osL = osSet(osL, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
osL = osSet(osL, 'time step', timeStep);

% Set osI data to raw pixel intensities of stimulus
% osI = osSet(osI, 'rgbData', 0.6-iStim.sceneRGB);
osL = osCompute(osL,absorptions);

% % Plot the photocurrent for a pixel
osPlot(osL,absorptions);
% os= osL;
%% Outer segment calculation - biophysical model

osBp = osCreate('BioPhys');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','m');
osBp = osSet(osBp, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
osBp = osSet(osBp, 'time step', timeStep);

% Set osI data to raw pixel intensities of stimulus
% osI = osSet(osI, 'rgbData', 0.6-iStim.sceneRGB);
osBp = osCompute(osBp,absorptions);

% % Plot the photocurrent for a pixel
osPlot(osBp,absorptions);
% os = osBp;
%% Find bipolar responses
os = osL;
bp = bipolar(os);

% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

bp = bipolarCompute(bp, os);

bipolarPlot(bp);

%% Find RGC responses

clear params innerRetinaBpSu
params.name      = 'Bipolar with nonlinear subunits'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 6;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaBpSu = irCreate(bp, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetinaBpSu.mosaicCreate('model','Subunit','type','on midget');

tCenterOrig = innerRetinaBpSu.mosaic{1}.mosaicGet('tCenter');

tCenterNew{1} = zeros(size(tCenterOrig)); 
tCenterNew{2} = zeros(size(tCenterOrig)); 
tCenterNew{3} = zeros(size(tCenterOrig));

tImpulse = 1000;
tCenterNew{1}(1) = tImpulse; tCenterNew{2}(1) = tImpulse; tCenterNew{3}(1) = tImpulse;
innerRetinaBpSu.mosaic{1}.mosaicSet('tCenter',tCenterNew);

tSurroundNew{1} = zeros(size(tCenterOrig)); 
tSurroundNew{2} = zeros(size(tCenterOrig)); 
tSurroundNew{3} = zeros(size(tCenterOrig));
tSurroundNew{1} = -tImpulse; tSurroundNew{2} = -tImpulse; tSurroundNew{3} = -tImpulse;
innerRetinaBpSu.mosaic{1}.mosaicSet('tSurround',tSurroundNew);

% irPlot(innerRetinaBpSu,'mosaic');
% Compute RGC mosaic responses

innerRetinaBpSu = irCompute(innerRetinaBpSu, bp);
irPlot(innerRetinaBpSu, 'psth');%,'cell',[4 4]);
% irPlot(innerRetinaBpSu, 'linear');

% irPlot(innerRetinaBpSu, 'linear','cell',[4 3]);
% irPlot(innerRetinaBpSu, 'raster','cell',[4 4]);

%% Bipolar object with no convolutional subunits, just impulse function
bp2 = bipolar(osL);

% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

bipolarRF = bipolarGet(bp2,'sRFcenter');

bipolarRFnew = zeros(size(bipolarRF));
bipolarRFnew(4,4) = 1;

bp2 = bipolarSet(bp2,'sRFcenter',bipolarRFnew);

bp2 = bipolarSet(bp2,'sRFsurround',bipolarRFnew);

bp2 = bipolarCompute(bp2, osL);

bipolarPlot(bp2);

%% Bipolar input with linear subunits

clear params innerRetinaBpImpSu
params.name      = 'Bipolar impulse with no subunits'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 7;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaBpImpNoSu = irCreate(bp2, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetinaBpImpNoSu.mosaicCreate('model','LNP','type','on midget');

irPlot(innerRetinaBpImpNoSu,'mosaic');
% Compute RGC mosaic responses

innerRetinaBpImpNoSu = irCompute(innerRetinaBpImpNoSu, bp);
% irPlot(innerRetinaBpImpNoSu, 'psth');
irPlot(innerRetinaBpImpNoSu, 'linear');
% irPlot(innerRetinaBpImpNoSu, 'raster');

%%
% % % % % % % % % % % % 
% % % % % % % % % % % % 

%% Compute RGC without subunits
% Input = RGB
osI = osCreate('displayRGB');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','m');
osI = osSet(osI, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
osI = osSet(osI, 'time step', timeStep);

osI = osSet(osI, 'rgbData', 2*(iStim.sceneRGB-0.5));


clear params innerRetinaNoBpNoSU
params.name      = 'No SU no bipolar'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 7;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaNoBpNoSU = irCreate(osI, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetinaNoBpNoSU.mosaicCreate('model','lnp','type','on midget');
irPlot(innerRetinaNoBpNoSU,'mosaic');

innerRetinaNoBpNoSU = irCompute(innerRetinaNoBpNoSU, osI);
irPlot(innerRetinaNoBpNoSU, 'linear');
% irPlot(innerRetinaNoBpNoSU, 'psth');

%% Compute RGC with subunits
% Input = RGB
osI = osCreate('displayRGB');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','m');
osI = osSet(osI, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
osI = osSet(osI, 'time step', timeStep);

osI = osSet(osI, 'rgbData', 2*(iStim.sceneRGB-0.5));


clear params innerRetinaNoBPSU
params.name      = 'Nonlinear SU no bipolar'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 7;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaNoBPSU = irCreate(osI, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetinaNoBPSU.mosaicCreate('model','subunit','type','on midget');
% irPlot(innerRetinaNoBPSU,'mosaic');

innerRetinaNoBPSU = irCompute(innerRetinaNoBPSU, osI);
irPlot(innerRetinaNoBPSU, 'linear');
irPlot(innerRetinaNoBPSU, 'psth');%,'cell',[3 3]);
