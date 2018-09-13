% Lens transmittance 
%
% Description:
%    Illustrate the effect of changing the lens transmittance.
%

% History:
%    xx/xx/xx  BW   Created it
%    09/11/18  jnm  Formatting

%% Initialization
ieInit

%% Create a colorful image
scene = sceneCreate;
scene = sceneSet(scene, 'fov', 1);
vcAddObject(scene);
sceneWindow;

%%  Make the standard WVF optics model
oi = oiCreate('wvf human');
oi = oiCompute(oi, scene);
vcAddObject(oi);
oiWindow;
oiPlot(oi, 'lens transmittance');

%% Modify and plot the lens density at 0.2
oi = oiSet(oi, 'lens density', 0.2);
oi = oiCompute(oi, scene);
vcAddObject(oi);
oiWindow;
oiPlot(oi, 'lens transmittance');

%% Modify and plot the lens density at 0
oi = oiSet(oi, 'lens density', 0);
oi = oiCompute(oi, scene);
vcAddObject(oi);
oiWindow;
oiPlot(oi, 'lens transmittance');
