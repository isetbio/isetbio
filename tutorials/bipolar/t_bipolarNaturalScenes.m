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

stimulusSelect = 2;
nSteps = 200;
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
        testmovie = data.testmovie.matrix;
        paramsMovie.nSteps = nSteps;
        paramsMovie.timeInterval = .0005;
        paramsMovie.expTime = .0005;
        iStim = ieStimulusMovie(testmovie(:,:,1:nSteps),paramsMovie); 
        
        absorptions = iStim.sensor; % cone isomerizations
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
bipolarPlot(bp);

%% Find RGC responses
% Build an IR object that takes as input the bipolar mosaic.

clear params innerRetinaBpSu
params.name      = 'Bipolar with nonlinear subunits'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 9;        % Radius in mm
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
