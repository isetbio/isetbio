%% Illustrate bipolar calculations and plots

%%
ieInit

%%
scene = sceneCreate('rings rays');
scene = sceneSet(scene,'fov',3);
oi = oiCreate; oi = oiCompute(oi,scene);
ieAddObject(oi);

%% Random noise is the default
cMosaic = coneMosaic;
cMosaic.emGenSequence(50);
cMosaic.setSizeToFOV(0.8*sceneGet(scene,'fov'));
cMosaic.compute(oi);
cMosaic.computeCurrent;

%%
bp = bipolar(cMosaic);
bp.compute(cMosaic);

%%
bp.plot('movie response');

%% 
bp.plot('spatial rf');

%% This is the default bipolar temporal filter

% The bipolarFilter routine tries to create a filter so that os convolved
% with bipolar matches the Pillow filter.  Maybe this should be a
% bipolar.get function.
bpFilter = bipolarFilter(bp,cMosaic,'graph',true);

%%