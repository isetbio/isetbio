% s_rgcGrid
%
% Description:
%    RGC grid test
%
%    Check that the grid appears as a grid on the RGC response.
%

% History:
%    XX/XX/17  BW   ISETBIO Team, 2017
%    06/11/19  JNM  Documentation pass

%%
ieInit

%% Scene, oi, grid, cone mosaic
imSize = 128;
lineSpacing = 48;
fov = 2; % deg
scene = sceneCreate('grid lines', imSize, lineSpacing);
scene = sceneSet(scene, 'fov', fov);
oi = oiCreate;  % Standard human optics
oi = oiCompute(oi, scene);
ieAddObject(oi);
oiWindow;

%%
nMovements = 25;
cMosaic = coneMosaic;
cMosaic.setSizeToFOV(fov);
cMosaic.emGenSequence(nMovements);
cMosaic.compute(oi);
cMosaic.computeCurrent;

%% Make the biplar layer with just one mosaic
clear bpL bpMosaicParams
bpL = bipolarLayer(cMosaic);

ii = 1;
bpMosaicParams.spread  = 2;  % RF diameter w.r.t. input samples
bpMosaicParams.stride  = 2;  % RF diameter w.r.t. input samples
bpL.mosaic{ii} = bipolarMosaic(cMosaic, 'on midget', bpMosaicParams);
bpL.mosaic{ii}.compute;

% bpL.window;
%% Make the RGC layer and show it
clear rgcL rgcParams
rgcL = rgcLayer(bpL);

% Spread and stride are not working
rgcParams.rfDiameter = 2;

% rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1}, 'on midget');
rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1}, 'on midget', rgcParams);
rgcL.compute;
rgcL.window;

%%