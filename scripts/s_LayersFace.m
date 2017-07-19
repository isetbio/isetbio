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

% Something wrong with onparasol offparasol that is OK with ondiffuse and
% offdifuse.  Look into it.
cellType = {'ondiffuse','offdiffuse','onmidget','offmidget','onsbc'};

clear bpL
clear bpMosaicParams

bpL = bipolarLayer(cMosaic);

bpMosaicParams.parent = bpL;
bpMosaicParams.rectifyType = 1;  % Experiment with this
bpMosaic = cell(5,1);
for ii=1:length(cellType)
    % Maybe this could be
    bpMosaicParams.cellType = cellType{ii};
    bpL.mosaicCreate(cellType{ii},bpMosaicParams);
    
    % Should be samples on the input layer
    support = 11;
    bpL.mosaic{ii}.set('sRFCenter',  fspecial('gaussian',[support support],1));
    bpL.mosaic{ii}.set('sRFSurround',fspecial('gaussian',[support support],2));
    
    bpL.mosaic{ii}.compute(cMosaic);     % Knows about cMosaic input
end
bpL.window;

%% Retinal ganlion cell model

clear rgcLayer

% Create retina ganglion cell layer object
rgcL = rgcLayer(bpL);
cellType = {'onparasol','offparasol','onmidget','offmidget','onsbc'};

% There are various parameters you could set.  We will write a script
% illustrating these later.  We need a description.
mosaicParams.centerNoise = 0;
mosaicParams.ellipseParams = [1 1 0];  % Principle, minor and theta
% mosaicParams.axisVariance = .1;
mosaicParams.type  = cellType;
mosaicParams.model = 'LNP';
mosaicParams.coupling = false;

diameters = round([10 10 5 5 20]);  % Should be samples on the input layer
for ii = 1:length(cellType)
    mosaicParams.rfDiameter = diameters(ii);
    mosaicParams.type = cellType{ii};
    mosaicParams.inMosaic = ii;   % Could switch up and match inputs to outputs
    rgcL.mosaicCreate(mosaicParams);
end
% rgcL.mosaic{1}.window
nTrials = 1; rgcL.set('nTrials',nTrials);

%% Compute the inner retina response and visualize

% Number of trials refers to number of repeats of the same stimulus
disp('Computing rgc responses');
[rgcL, nTrialsSpikes] = rgcL.compute(bpL.mosaic, ...
    'bipolarScale',50,...
    'bipolarContrast',0.2);

%% Show the layer

rgcL.window;
