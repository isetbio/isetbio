% t_rgc
% 
% A tutorial that describes a simple implementation of the rgcLayer object
% in isetbio. A sensor/cone mosaic and outersegment are loaded from the
% Remote Data Toolbox. The stimulus consists of a bar sweeping from left to
% right. The spatial nature of the stimulus allows the RGC response to be
% visualized easily. An rgcLayer object is created, and the rgc mosaic
% responses are computed. Several plotting operations are demonstrated.
% 
% 02/2016 JRG (c) isetbio team

%% Load sensor/cone mosaic and os
% The stimulus is a dynamic scene that consists of a bar sweeping from left
% to right. The scene, oi, sensor/cone mosaic and outer segment have been
% precomputed and are loaded using the Remote Data Tooblox.

% Initialize RDT
rdt = RdtClient('isetbio');
client.credentialsDialog();
rdt.crp('resources/data/rgc');
% Load data
data = rdt.readArtifact('t_rgcData', 'type', 'mat');
% Separate cone mosaic/sensor and outer segment.
coneMosaic = data.coneMosaic;
os = data.os;

%% Build RGC
% The RGC layer object is constructed by specifying input parameters. The
% user sets the name of the object, the type of model (linear, LNP, GLM,
% etc.) and the position of the retinal patch (which eye, radius from
% fovea and polar angle). The RGC determines the size of spatial receptive
% fields based on the temporal equivalent eccentricity calculated from th
% patch location, and build RGCs with spatial RFs over the cone mosaic and
% temporal impulse responses with the appropriate sampling rate.

clear params
% Set the parameter values
params.name    = 'Macaque inner retina 1'; % This instance
params.model   = 'GLM';    % Computational model
% Determined by user
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees
% Inherited from cone mosaic; may be incorporated into just sensor
params.row     = sensorGet(coneMosaic,'row');  % N row samples
params.col     = sensorGet(coneMosaic,'col');  % N col samples
params.spacing = sensorGet(coneMosaic,'width','um'); % Cone width
params.timing  = sensorGet(coneMosaic,'time interval','sec'); % Temporal sampling

% Create the rgc layer object
rgc1 = rgcCreate(params);

%% Build RGC mosaics
% The mosiac property of the RGC layer object stores the mosaics of the
% different types of RGCs. The code currently supports the five most common
% types: ON parasol, OFF parasol, ON midget, OFF midget and small
% bistratified. Mosaics can be added individually or all five may be added
% automatically using the loop structure below.

numberMosaics = 5; % setting to 5 generates all five of the most common types.
for cellTypeInd = 1:numberMosaics
    rgc1 = rgcMosaicCreate(rgc1);
end

% Alternative syntax for creating single layers at a time
% rgc1 = rgcMosaicCreate(rgc1,'mosaicType','onMidget');
         
%% Compute the RGC responses
rgc1 = rgcCompute(rgc1, os);

% Plot various aspects of the RGC response
% rgcPlot(rgc1, 'mosaic');
% rgcPlot(rgc1, 'rasterResponse');
rgcPlot(rgc1, 'psthResponse');

% Create a movie of the response
% rgcMovie(rgc1, os);
