%% s_LayersFace.m
%
% Testing the new layer architecture with face image of LC
%
% BW ISETBIO Team, 2017

%%
ieInit;

%% Get the data from the RDT site

rdt = RdtClient('isetbio');
rdt.crp('/resources/data/cmosaics');
% rdt.listArtifacts('type','mat','print',true);

% coneMosaicDataFace
% For some reason, this opens the cMosaic.window.  I don't understand why.
% This happens when we load the data that is downloaded in
% rdtLoadWellKnownFileTypes.
data      = rdt.readArtifact('coneMosaicDataFace', 'type', 'mat');
cMosaic   = data.cMosaic;

% cMosaic.window;

%% Create a set of bipolar cell types in the bipolar mosaic

cellType = {'ondiffuse','offdiffuse','onmidget','offmidget','onSBC'};

clear bpL
clear bpMosaicParams

bpL = bipolarLayer(cMosaic);

bpMosaicParams.parent = bpL;
bpMosaicParams.rectifyType = 1;  % Experiment with this
bpMosaic = cell(5,1);
for ii=1:length(cellType)
    % Maybe this could be
    %   bpMosaicParams.cellType = cellType{ii};
    %   bpL.mosaicCreate(ii,bpMosaicParams);
    %   bpL.mosaic{ii}.compute;     % Knows about cMosaic input
    %
    
    bpMosaicParams.cellType = cellType{ii};
    bpMosaic{ii} = bipolarMosaic(cMosaic, bpMosaicParams);
    
    % Set up the bipolar RF 
    bpMosaic{ii}.set('sRFcenter',3);
    bpMosaic{ii}.set('sRFsurround',5);
    
    bpMosaic{ii}.compute(cMosaic);
    bpL.mosaic{ii} = bpMosaic{ii};
end
bpL.window;

%% Retinal ganlion cell model

clear rgcLayer

% Create retina ganglion cell layer object
rgcL = rgcLayer(bpL);

% There are various parameters you could set.  We will write a script
% illustrating these later.  We need a description.
mosaicParams.centerNoise = 0;
mosaicParams.ellipseParams = [1 1 0];  % Principle, minor and theta
% mosaicParams.axisVariance = .1;
mosaicParams.type  = cellType;
mosaicParams.model = 'GLM';

% diameters = [15 15 7 7 20];  % In microns.
diameters = 2*[15 15 7 7 20];  % In microns.

cellType = {'on parasol','off parasol','on midget','off midget','smallbistratified'};
for ii = 1:length(cellType)
    mosaicParams.rfDiameter = diameters(ii);
    mosaicParams.type = cellType{ii};
    mosaicParams.inMosaic = 1;   % Could switch up and match inputs to outputs
    rgcL.mosaicCreate(mosaicParams);
end
% rgcL.mosaic{1}.window
nTrials = 1; rgcL.set('numberTrials',nTrials);

%% Compute the inner retina response and visualize

% Number of trials refers to number of repeats of the same stimulus
disp('Computing rgc responses');
[rgcL, nTrialsSpikes] = rgcL.compute(bpMosaic,...
    'bipolarScale',50,...
    'bipolarContrast',0.2);

%% Show the layer

rgcL.window;
