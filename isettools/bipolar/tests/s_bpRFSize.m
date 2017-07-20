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
bpMosaicParams.parent = bpL;
bpMosaicParams.stride = 3;
bpL.mosaicCreate('onmidget',bpMosaicParams);

% Set the smaller size
support = 3;
bpL.mosaic{1}.set('sRFCenter',  fspecial('gaussian',[support support],1));
bpL.mosaic{1}.set('sRFSurround',fspecial('gaussian',[support support],2));
bpL.mosaic{1}.compute(cMosaic);     % Knows about cMosaic input

%% Now make a larger bipolar

bpL.mosaicCreate('onmidget',bpMosaicParams);

% Set the larger size
support = 7;
bpL.mosaic{2}.set('sRFCenter',  fspecial('gaussian',[support support],3));
bpL.mosaic{2}.set('sRFSurround',fspecial('gaussian',[support support],5));
bpL.mosaic{2}.compute(cMosaic);     % Knows about cMosaic input

bpL.window;

%%
