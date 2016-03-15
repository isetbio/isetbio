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

%% Build the inner retina object

clear params
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina1 = irCreate(osL, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetina1.mosaicCreate('model','glm','type','on midget');


%% Compute RGC mosaic responses

innerRetina1 = irCompute(innerRetina1, osL);
irPlot(innerRetina1, 'psth response');

%% Show me the PSTH for one particular cell

irPlot(innerRetina1, 'psth response','cell',[2 2]);
title('OS Linear and Coupled GLM');

%% Compute the outer segment response

% In this case we use a linear model.  Below we use a more complex model
osB = osCreate('bioPhys');

% Set up the 
patchSize = sensorGet(absorptions,'width','um');
osB = osSet(osB, 'patch size', patchSize);

timeStep = sensorGet(absorptions,'time interval','sec');
osB = osSet(osB, 'time step', timeStep);

osB = osCompute(osB, absorptions);

%% Compute RGC mosaic responses

innerRetina2 = irCreate(osB, params);
innerRetina2.mosaicCreate('model','glm','type','on midget');

innerRetina2 = irCompute(innerRetina2, osB);
irPlot(innerRetina2, 'psth response');

%% Show me the PSTH for one particular cell

irPlot(innerRetina2, 'psth response','cell',[2 2]);
title('OS Biophys and Coupled GLM');

%%
