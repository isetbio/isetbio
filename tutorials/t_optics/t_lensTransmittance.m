% Lens transmittance and a little macular pigment, too.
%
% Description:
%    Illustrate the effect of changing the lens transmittance.  Also
%    illustrates the macular.
%
% BW
%
% See also
%

%% Initialization
ieInit

%% Create a colorful image
scene = sceneCreate;
scene = sceneSet(scene, 'fov', 1);
ieAddObject(scene); sceneWindow;

%% The macular pigment is normally attached to the cones
%
% Here we illustrate the effect of that pigment by including the
% typical lens absoprtance and then multiplying by the macular
% absorptance.  We stick the whole thing into the OI and image.

% Create the OI
oi = oiCreate();

% Make sure we are set correctly
oi = oiSet(oi,'lens density',1);
oi = oiCompute(oi, scene);

% Get the macular pigment
Mac = Macular;
photons = oiGet(oi,'photons');
[photons, row, col] = RGB2XWFormat(photons);
photons = photons*diag(Mac.transmittance);

% Put it back
photons = XW2RGBFormat(photons,row,col);
oi = oiSet(oi,'photons',photons);
oi = oiSet(oi,'name','Lens and macular');
ieAddObject(oi); oiWindow;

%%  Make the standard WVF optics model
oi = oiCreate('wvf human');
oi = oiCompute(oi, scene);
oi = oiSet(oi,'name','Lens 1');
ieAddObject(oi); oiWindow;
oiPlot(oi, 'lens transmittance');

%% Modify and plot the lens density at 0.3
oi = oiSet(oi, 'lens density', 0.3);
oi = oiCompute(oi, scene);
oi = oiSet(oi,'name','Lens 0.3');

ieAddObject(oi); oiWindow;
oiPlot(oi, 'lens transmittance');

%% Modify and plot the lens density at 0
oi = oiSet(oi, 'lens density', 0);
oi = oiCompute(oi, scene);
oi = oiSet(oi,'name','Lens 0');
ieAddObject(oi); oiWindow;
oiPlot(oi, 'lens transmittance');


%%

