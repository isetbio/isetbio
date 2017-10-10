%% t_bipolar
% Illustrate bipolar layer and mosaic calculations.  Show some plots.
% 
% TODO:  Notice the ringing at the end of the bipolar temporal impulse response!  
% That's no good.

ieInit
%% Test scene
% Resolution target spanning.
%%
scene = sceneCreate('rings rays');
scene = sceneSet(scene,'fov',2);   % deg horizontal
oi = oiCreate; 
oi = oiCompute(oi,scene);

% To see the retinal irradiance use
%  ieAddObject(oi); oiWindow
%% Random noise is the default
%%
cMosaic = coneMosaic;
cMosaic.emGenSequence(50);
cMosaic.setSizeToFOV(0.8*sceneGet(scene,'fov'));
cMosaic.compute(oi);
cMosaic.computeCurrent;

% To see the cone mosaic use
%  cMosaic.window;
%%
clear bpL bpMosaicParams
bpL = bipolarLayer(cMosaic);

ii = 1;
bpMosaicParams.spread  = 2;  % RF diameter w.r.t. input samples
bpMosaicParams.stride  = 2;  % RF diameter w.r.t. input samples
bpL.mosaic{ii} = bipolarMosaic(cMosaic,'on midget',bpMosaicParams);
bpL.mosaic{ii}.compute;
%% Open the bipolar layer window
%%
bpL.window;
%% The bipolar spatial receptive field
%%
vcNewGraphWin;
bpL.mosaic{1}.plot('spatial rf');
%% The  bipolar temporal filter
% The bipolarFilter routine tries to create a filter so that os convolved with 
% bipolar matches the Pillow filter.  
%%
bpFilter = bipolarFilter(bpL.mosaic{1},cMosaic,...
    'graph',true);
%% 
%