%% In which our heroes lay out the basic architecture of the ir class
%
% Aspirational:
%    function [scene, sceneRGB, oi, sensor] = movieCreate(varargin)
%
% JG/BW ISETBIO Team, Copyright 2015
% (HJ) ISETBIO TEAM, 2014
% (JRG) modified 10/2015
%

%%
ieInit

%% Movie of the cone absorptions 

rd = RdtClient('isetbio');
rd.crp('/resources/data/istim');
a = rd.listArtifacts;

%%
whichA =1 ;
data = rd.readArtifact(a(whichA).artifactId);
iStim = data.iStim;
absorptions = iStim.absorptions;

% coneImageActivity(absorptions,'step',1,'dFlag',true);

%% Compute the outer segment response

% In this case we use a linear model.  Below we use a more complex model
osL = osCreate('linear');

% Set up the 
patchSize = sensorGet(absorptions,'width','um');
osL = osSet(osL, 'patch size', patchSize);

timeStep = sensorGet(absorptions,'time interval','sec');
osL = osSet(osL, 'time step', timeStep);

osL = osCompute(osL, absorptions);


%%

osD = osCreate('displayRGB');

% Set up the 
patchSize = sensorGet(absorptions,'width','um');
osD = osSet(osD, 'patch size', patchSize);

timeStep = sensorGet(absorptions,'time interval','sec');
osD = osSet(osD, 'time step', timeStep);

osD = osSet(osD,'rgbData',iStim.sceneRGB);
%% Build the inner retina object

clear params innerRetina
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina = irCreate(osL, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetina.mosaicCreate('model','glm','type','on midget');
innerRetina.mosaic{1}.mosaicSet('numberTrials',120);
%% Here is the layout of the RGC receptive fields on retinal surface

irPlot(innerRetina, 'mosaic');

%% Compute RGC mosaic responses

innerRetina = irCompute(innerRetina, osD);

%% Look at the linear inputs to the cells before spiking

irPlot(innerRetina, 'response linear');

%% Show me the raster plots for all the cells in the mosaic

irPlot(innerRetina, 'raster response');

%% Show me the PSTH for one particular cell

irPlot(innerRetina, 'psth response');

%%