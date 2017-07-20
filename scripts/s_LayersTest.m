%% s_LayersTest.m
%
% Testing the new layer architecture.
%
% BW ISETBIO Team, 2017

%%
ieInit;

%% Get the data from the RDT site

rdt = RdtClient('isetbio');
rdt.crp('/resources/data/cmosaics');
% rdt.listArtifacts('type','mat','print',true);

% coneMosaicDataFixedEye;
% vaConeMosaic
% coneMosaicDataScones
% coneMosaicDataF
data      = rdt.readArtifact('vaConeMosaic', 'type', 'mat');
cMosaic   = data.cMosaic;
alignedC  = data.alignedC;

% cMosaic.window;

%% Create a set of bipolar cell types in the bipolar mosaic

bpL = bipolarLayer(cMosaic);

% Make each type of bipolar mosaic
cellType = {'on diffuse','off diffuse','on midget','off midget','on sbc'};

clear bpMosaicParams
bpMosaicParams.rectifyType = 1;  % Experiment with this
bpMosaicParams.spread  = 1;  % RF diameter w.r.t. input samples
bpMosaicParams.stride  = 1;  % RF diameter w.r.t. input samples

bpMosaic  = cell(1,length(cellType));
bpNTrials = cell(1,length(cellType));
for ii = 1:length(cellType)   
    bpMosaicParams.cellType = cellType{ii};
    % We want this to be mosaicCreate(bpMosaicParams) in the end.
    bpL.mosaic{ii} = bipolarMosaic(cMosaic,bpMosaicParams);
    bpL.mosaic{ii}.compute(cMosaic);   
end

bpL.window;

%% Retinal ganlion cell model

% Tried only for GLM and LNP.  Ask JRG if other types are still supported.
clear rgcL rgcParams

% Create retina ganglion cell layer object
rgcL = rgcLayer(bpL);

% There are various parameters you could set.  We will write a script
% illustrating these later.  We need a description.
rgcParams.centerNoise = 0;
rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta
rgcParams.model = 'GLM';
% mosaicParams.axisVariance = .1;

diameters = [5 5 3 3 10];  % In microns.

cellType = {'on parasol','off parasol','on midget','off midget','onsbc'};
for ii = 1:length(cellType)
    rgcParams.rfDiameter = diameters(ii);
    % Input bipolar mosaic and specific parameters
    rgcL.mosaicCreate(bpL.mosaic{ii},cellType{ii},rgcParams);
end

nTrials = 1; rgcL.set('numberTrials',nTrials);

%% Compute the inner retina response and visualize

% Every mosaic has its input and properties assigned so we should be able
% to just run through all of them.
disp('Computing rgc responses');
rgcL = rgcL.compute(...
    'bipolarScale',50,...
    'bipolarContrast',1);

%%
rgcL.window;

%%