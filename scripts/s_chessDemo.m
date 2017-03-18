%% s_chessDemo
%
% Try for RGC calculation on chess set
% 4 deg FOV, so kind of big
%

%% 
ieInit

%%
load(fullfile(isetbioRootPath,'local','chessView1'),'scene');
scene1 = scene;
scene1 = sceneSet(scene1,'fov',2);
vcAddObject(scene1); sceneWindow;

oi1 = oiCreate;
oi1 = oiCompute(oi1,scene1);
vcAddObject(oi1); oiWindow;

%%
load(fullfile(isetbioRootPath,'local','chessView2'),'scene');
scene2 = scene;
scene2 = sceneSet(scene2,'fov',2);
vcAddObject(scene2); sceneWindow;

oi2 = oiCreate;
oi2 = oiCompute(oi2,scene2);
vcAddObject(oi2); oiWindow;

%%

cMosaic = coneMosaic;
cMosaic.emGenSequence(50);

cMosaic.setSizeToFOV(2);

cMosaic.compute(oi1);
cMosaic.computeCurrent;

cMosaic.window;
