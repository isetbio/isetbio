%% s_bpGrid
%
% Create a scene of grid lines and check that they render on the bipolar
% cell mosaic as a grid.
%
% BW, ISETBIO Team, 2017

%%
ieInit

%% Scene and oi grid

imSize = 128;
lineSpacing = 48;
fov = 1; % deg
scene = sceneCreate('grid lines',imSize,lineSpacing);
scene = sceneSet(scene,'fov',fov);

oi = oiCreate;    % Standard human optics
oi = oiCompute(oi,scene);

%%  Make the mosaic

nMovements = 50;
cMosaic = coneMosaic;
cMosaic.setSizeToFOV(fov);
cMosaic.emGenSequence(nMovements);

cMosaic.compute(oi);
cMosaic.computeCurrent;

% cMosaic.window;

%% Make the biplar layer with just one mosaic 

clear bpL bpMosaicParams

bpL = bipolarLayer(cMosaic);

% Spread and stride are not working

ii = 1;
bpMosaicParams.spread  = 2;  % RF diameter w.r.t. input samples
bpMosaicParams.stride  = 2;  % RF diameter w.r.t. input samples
bpL.mosaic{ii} = bipolarMosaic(cMosaic,'on midget',bpMosaicParams);

% Maybe we should have a bpL.compute
bpL.mosaic{ii}.compute;
bpL.window;

%%
clear rgcL rgcParams
rgcL = rgcLayer(bpL);

% Spread and stride are not working
rgcParams.rfDiameter = 1.5;

ii = 1;
% rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1},'on midget');
rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1},'on midget',rgcParams);
% Maybe we should have a bpL.compute
rgcL.compute;
rgcL.window;
