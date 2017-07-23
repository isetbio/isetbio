%% s_bpRFSize.m
%
% Show how to set the spatial receptive field of the bipolar
%
% BW ISETBIO Team, 2017

%%
ieInit;

%% Get data from the RDT site
rdt = RdtClient('isetbio');
rdt.crp('/resources/data/cmosaics');
% rdt.listArtifacts('print',true);

data      = rdt.readArtifact('coneMosaicDataFixedEye', 'type', 'mat');
cMosaic   = data.cMosaic;

% To view the mosaic
% cMosaic.window;

%% Create a set of bipolar cell types in the bipolar mosaic

clear bpL
clear bpMosaicParams

bpL = bipolarLayer(cMosaic);

% Set the smaller size
bpMosaicParams.stride = 1;
bpMosaicParams.spread = 1;

bpL.mosaic{1} = bipolarMosaic(cMosaic,'onmidget',bpMosaicParams);

bpL.mosaic{1}.compute;     % Knows about cMosaic input

%% Now make a larger bipolar

% Set the larger size
bpMosaicParams.spread = 3;
bpMosaicParams.stride = 3;

bpL.mosaic{2} = bipolarMosaic(cMosaic,'onmidget',bpMosaicParams);

bpL.mosaic{2}.compute();     % Knows about cMosaic input

bpL.window;

%%
