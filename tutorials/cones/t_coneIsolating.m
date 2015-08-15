%% t_coneIsolating
%
% Illustrate how to create cone-isolating color patches.  Illustrate
% limits of the cone isolation for patterns
%
% HJ/BW Copyright ISETBIO Team, 2015

% Initialize a new ISETBIO session
ieInit;

%% Choose a display for simulation
% We start with an Sony OLED display that we once calibrated.
dpi = 110;
d = displayCreate('OLED-Sony','dpi',dpi);

viewDist = 0.3;  % Set the subject's viewing distance in meters
d = displaySet(d, 'viewing distance', viewDist);


%% Create the RGB images for testing
% The RGB images here are actually DAC, not linear RGB. They will be
% converted into linear RGB by display gamma in vcReadImage line 180.
% Though inputs should be DAC values, they can still be in range 0~1. They
% will be converted to quantized levels in vcReadImage line 160~176
p.pattern = 0.5*ones(1,65); p.pattern(23:43) = 1;

% This converts the linear RGB values into LMS.
% This maps [r,g,b]* rgb2lms = [L,M,S]
% rgb2lms = displayGet(d,'rgb2lms');
sConeIsolating = unitLength([0 0 1]*displayGet(d,'lms2rgb'));

img = image1d(p.pattern,'rgb',sConeIsolating,'mean',0.5);

vcNewGraphWin(); 
imshow(img); axis image; 

%% We fine the S-cone isolating direction and use the RGB in that direction

% We are calculating as if we have a linear gamma table, for simplicity
dLinear = displaySet(d,'gtable','linear');
scene = sceneFromFile(img, 'rgb', [], dLinear); % mis-aligned
scenePlot(scene,'radiance hline',[1,round(sz(1)/2)]);
ieAddObject(scene); sceneWindow;

%% Now render the scene through the human optics onto a cone mosaic

oi = oiCreate('wvf human');
oi = oiCompute(oi,scene);
ieAddObject(oi); oiWindow;

cones = sensorCreate('human');
cones = sensorSetSizeToFOV(cones,sceneGet(scene,'fov'),scene,oi);
cones = sensorSet(cones,'noiseflag',0);
cones = sensorCompute(cones,oi);
ieAddObject(cones); sensorWindow('scale',1);

row = round(sensorGet(cones,'rows')/2);
sensorPlot(cones,'photons hline',[row,1]);

% Notice the little tick in the L and M cones at the edge.  Do you know why
% this arises?

%% END
