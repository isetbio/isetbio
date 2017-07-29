%% Lens transmittance 
%
% Illustrate the effect of changing the lens transmittance
%
% BW

%%
ieInit

%% Create a colorful image
scene = sceneCreate;
scene = sceneSet(scene,'fov',1);
vcAddObject(scene); sceneWindow;

%%  Make the standard WVF optics model

oi = oiCreate('wvf human');
oi = oiCompute(oi,scene);
vcAddObject(oi); oiWindow;
oiPlot(oi,'lens transmittance');

%%
oi = oiSet(oi,'lens density',0.2);
oi = oiCompute(oi,scene);
vcAddObject(oi); oiWindow;
oiPlot(oi,'lens transmittance');

%%
oi = oiSet(oi,'lens density',0);
oi = oiCompute(oi,scene);
vcAddObject(oi); oiWindow;
oiPlot(oi,'lens transmittance');
