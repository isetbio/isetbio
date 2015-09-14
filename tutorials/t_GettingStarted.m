%% t_GettingStarted
%
% Welcome to ISETBIO!!!
%
% This tutorial is a very high level description of the kinds of
% calculations you can make using this tool.  There are many other
% tutorials (which we reference below) that go into more detail about
% how to change parameters and do other types of calculations.
%
% HJ/BW ISETBIO Team, Copyright 2015

%% First, we initialize the isetbio environment

% This sets up some general variables that isetbio functions use.
% Whenever you run this command, it resets the isetbio compute environment
ieInit

%% Let's a simple scene

% A scene describes the full spectral radiance of the image, the distance
% to each point in the image, and the spectral power distribution of the
% illuminant.
s = sceneCreate('radial lines');
s = sceneSet(s,'fov',1);
vcNewGraphWin;
sceneShowImage(s);

% There are many other types of scenes, and in fact, you can produce a
% scene from an rgb image and the specification of a display.  To learn
% more Please see 
%
%   t_sceneIntroduction 
%

%% Now let's create some human optics

% This is equivalent to oiCreate('human').
oi = oiCreate;

% To see more about controlling the properties of the optics have a look at
%
%    t_oiIntroduction
%

%% Now we compute the retinal irradiance

% The compute function combines the scene with the optics and stores the
% results in the optical image (oi) data structure
oi = oiCompute(oi,s);

% Let's have a look
vcNewGraphWin;
oiShowImage(oi);

%% Now, let's compute the cone absorptions for a sample human retina

% This is a default spatial array of human cones on a rectangular lattice.
cones = sensorCreate('human');

% We combine the oi and the cones to calculate the cone absorption rates
cones = sensorCompute(cones,oi);

% The lightness of each point indicates the cone absorption rate
vcNewGraphWin;
img = coneImageActivity(cones);
imagesc(img); axis off; axis image

%% Quantitative plots

% This being a computer, and us being fanatical about units, means that we
% can tell you stuff you might need to know about light levels, optics, and
% cone absorption rates.  For example

% The mean luminance
sceneGet(s,'mean luminance')

% Properties of the optics, such as the point spread function at 550nm
oiPlot(oi,'psf 550')

% The photon absorptions across a line of the cone mosaic.  The line is
% about 1/4 of the way down the rows.  The different panels show the L,M
% and S sensors.
sz = sensorGet(cones,'size');
sensorPlot(cones,'photons hline',[round(sz(1))*.25,1]);


%% END

