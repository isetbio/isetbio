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
        iStim = ieStimulusGratingSubunit(stimP);

        % load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/iStim_subunitGrating.mat');
%         load('/Users/james/Documents/MATLAB/isetbio misc/iStim_subunitGratingLong.mat')
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

%% Find bipolar responses

bp = bipolar(osL);

% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

bp = bipolarCompute(bp, osL);

bipolarPlot(bp);

%% Find RGC responses

clear params innerRetinaBpSu
params.name      = 'Bipolar with nonlinear subunits'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 7;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaBpSu = irCreate(bp, params);

% Create a coupled GLM model for the off parasol ganglion cell parameters
innerRetinaBpSu.mosaicCreate('model','Subunit','type','on midget');

irPlot(innerRetinaBpSu,'mosaic');

% Compute RGC mosaic responses
innerRetinaBpSu = irCompute(innerRetinaBpSu, bp);

irPlot(innerRetinaBpSu,'linear')


%% Bipolar object with no convolutional subunits, just impulse function
bp2 = bipolar(osL);

% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

bipolarRF = bipolarGet(bp2,'sRF');

bipolarRFnew = zeros(size(bipolarRF));
bipolarRFnew(4,4) = 1;

bp2 = bipolarSet(bp2,'sRF',bipolarRFnew);

bp2 = bipolarCompute(bp2, osL);

bipolarPlot(bp2);

%% Bipolar input with linear subunits

clear params innerRetinaBpImpSu
params.name      = 'Bipolar impulse with no subunits'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 7;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaBpImpNoSu = irCreate(bp2, params);

% Create a coupled GLM model for the off parasol ganglion cell parameters
innerRetinaBpImpNoSu.mosaicCreate('model','LNP','type','off parasol');

irPlot(innerRetinaBpImpNoSu,'mosaic');
% Compute RGC mosaic responses

innerRetinaBpImpNoSu = irCompute(innerRetinaBpImpNoSu, bp);
% irPlot(innerRetinaBpImpNoSu, 'psth');
irPlot(innerRetinaBpImpNoSu, 'linear');
% irPlot(innerRetinaBpImpNoSu, 'raster');
