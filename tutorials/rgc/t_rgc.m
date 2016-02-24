%% A simple implementation and tutorial of the inner retina mosaics
%
% A sensor/cone mosaic and outersegment are loaded from the Remote Data
% Toolbox. The stimulus consists of a bar sweeping from left to right. The
% spatial nature of the stimulus allows the RGC response to be visualized
% easily. An rgcLayer object is created, and the rgc mosaic responses are
% computed. Several plotting operations are demonstrated.
% 
% 02/2016 JRG (c) isetbio team

%%
ieInit

%% Load sensor/cone mosaic and os

% The stimulus is a dynamic scene that consists of a bar sweeping from left
% to right. The scene, oi, sensor/cone mosaic and outer segment have been
% precomputed and are loaded using the Remote Data Tooblox (RDT)

% Initialize RDT
rdt = RdtClient('isetbio');
rdt.crp('resources/data/rgc');

% Read the artifact containing a description of the cone mosaic
% (coneMosaic) and the biophysical properties of the outer segment (os).
% Because these are stored in a matlab file, the data are returned.
data       = rdt.readArtifact('t_rgcData', 'type', 'mat');
coneMosaic = data.coneMosaic;
os         = data.os;

%% Visualize the cone mosaic absoprtion pattern
%
% cMosaic25 = coneMosaic;
% cMosaic25.data.volts = coneMosaic.data.volts(:,:,25);
% vcAddObject(cMosaic25);sensorWindow('scale',true);
%
% Or make a color movie
%
% coneImageActivity(coneMosaic,'step',1,'dFlag',true);
% 

%% Build the inner retina object
%
% The inner retina holds the ganglion cell mosaics. 
%
% The user names the object, specifies the type of model (linear, LNP,
% GLM, etc.) and the position of the retinal patch (which eye, radius from
% fovea and polar angle). The RGC determines the size of spatial receptive
% fields based on the temporal equivalent eccentricity calculated from th
% patch location, and build RGCs with spatial RFs over the cone mosaic and
% temporal impulse responses with the appropriate sampling rate.

% Create the inner retina that holds the rgc mosaics
innerRetina = irCreate(os, 'name','Macaque inner retina');

%% Build RGC mosaics 
% 
% This is now handled internally in rgcCompute with the same call to
% rgcMosaicCreate. This ensures the mosaics have the correct properties
% according to those set in the rgc parent object.
% 
% The mosaic property of the RGC layer object stores the mosaics of the
% different types of RGCs. The code currently supports the five most common
% types: ON parasol, OFF parasol, ON midget, OFF midget and small
% bistratified. Mosaics can be added individually or all five may be added
% automatically using the loop structure below. If no mosaics are added
% manually, the single default mosaic added is on parasol.

% Alternative syntax for creating single layers at a time
% innerRetina = rgcMosaicCreate(innerRetina,'model','glm','type','on parasol');

% innerRetina = rgcMosaicCreate(innerRetina,'model','lnp','type','on midget');

% innerRetina.mosaicCreate('model','glm','type','on parasol');
innerRetina.mosaicCreate('model','glm','type','on midget');
innerRetina.mosaicCreate('model','glm','type','on parasol');
% innerRetina.mosaicCreate('model','linear','type','on midget');
%% Compute the RGC responses

% Compute linear and nonlinear responses
% innerRetina = irCompute(innerRetina, os);
innerRetina.compute(os);

%%
% innerRetina = innerRetina.computeContinuous(os);
% 
% % Compute spiking response
% numberTrials = 5;
% for repititions = 1:numberTrials
% %     innerRetina = irSpikeCompute(innerRetina);
%     innerRetina.computeSpikes();
% end

%% Plot various aspects of the RGC response
% irPlot(innerRetina, 'mosaic');
irPlot(innerRetina, 'raster');
% irPlot(innerRetina, 'psthResponse');

% irPlot(innerRetina,'psth','type','onParasol');
% irPlot(innerRetina,'psth','cell',[1 1]);
% irPlot(innerRetina,'psth','type','onParasol','cell',[1 1]);
% 
% Create a movie of the response
% rgcMovie(rgc1, os);
% rgcMovieWave;

%%
