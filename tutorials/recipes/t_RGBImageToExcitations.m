%% t_RGBImageToExcitations
%
% Description:
%    Create a scene from an RGB image and a "typical" display,
%    compute retinal image from the scene, and the photopigment
%    excitations of a patch of retinal cone mosaic.
%    
%    This is an overview tutorial that illustrates a core ISETBio
%    processing pipeline, from stimulus to the initial daylight visual
%    encoding. This tutorial is meant to convey how this core pipeline
%    looks and to expose a few of its parameters, as a point of departure.
%
%    In other tutorials we will unpack in more detail the processing
%    elements that make up this chain.
%
%    This tutorial is available in the github repository 
%      https://github.com/ISET/ISETTutorials.git
%    It is in folder
%      ISETBio/RGBImageToExcitations/t_RGBImageToExcitations.m
%
%    You will need ISETCam and ISETBio installed, and you should also put
%    the ISETTutorials repository on your path.
%       https://github.com/ISET/ISETCam.git
%       https://github.com/ISETBio/ISETBio.git
%    Go to the wiki pages for these two repositories for detailed
%    instructions on how to install these two packages - it's easy!
%
%    If you use the ToolboxToolbox (see https://github.com/ToolboxHub/ToolboxToolbox.git), you
%    can configure in one step with the Matlab command
%      tbUse('ISETTutorials');

%% Initialize and clear
clear; close all; ieInit;

%% Set sample wavelengths
wave = 400:5:700;

%% Create scene object from image file
% 
% Here will will build the scene from a stock RGB image, and
% a "typical" display.  We need to assume some display in order
% to interpret the RGB values in physical units.  ISETCam/Bio
% has data on a number of displays. Here we use measurements we
% made at some point of an Apple LCD display.  If you were modeling
% an experiment you did on a calibrated display, you would specify
% a file describing your display.  A later tutorial will discuss
% the ISETCam/Bio display object in more detail.
%
% Although sceneFromFile will read the scene directly and create
% the image, we sometimes want a higher pixel resolution so we read and
% upsample in advance of calling sceneFromFile. Setting upsampleFactor
% to 1 leaves the pixel resolution alone.
%
% We also explicitly set the width of the field of view in degrees,
% and the mean scene luminance, since neither of these comes in 
% automatically with the RGB image.  Luminance is in cd/m2. The
% fact that we set this illustrates the fact that we always do our
% best to work in physical units in ISETBio.
%
% Available images that will work in this section include:
%   'macbeth.tif', 'bears.png'. 'eagle.jpg', 'hats.jpg', 'lion-cardinal.png'
%
% Set parameters
imgFileName = 'zebra.jpg';              
dispLCDFile = 'LCD-Apple.mat';
upsampleFactor = 1;
viewingDistanceMeters = 1;
fieldOfViewWidthDegs = 2;
sceneLuminanceCdM2 = 200;

% Read and resize the RGB image.  You don't need to resize this,
% just doing it here to illustrate the idea that the conversion
% to scene works with however may pixels you have.
im = imread(imgFileName);
im = imresize(im,upsampleFactor);

% Routine sceneFromFile does the work.  It undoes the
% monitor's gamma correction to linearize the image data,
% and then converts to spectra using the monitor's channel
% spectra.
sc = sceneFromFile(im,'rgb',[],dispLCDFile,wave);

% Then set the scene parameters we want to adjust.  These
% parameters aren't associated with a typical RGB image,
% so we move ourselves into physical units by setting them
% right after we create the scene from the RGB image.
%
% The set on 'fov' determines the width and takes the distance into
% account. The height is determined by the aspect ratio of the underlying
% image.
sc = sceneSet(sc,'name','Scene on LCD');
sc = sceneSet(sc,'distance',viewingDistanceMeters);
sc = sceneSet(sc,'fov',fieldOfViewWidthDegs);
sc = sceneSet(sc,'mean luminance',sceneLuminanceCdM2);
scSizeDegs = [sceneGet(sc,'fov horizontal') sceneGet(sc,'fov vertical')];
fprintf('Retinal image size is %0.1f by %0.1f degrees\n',scSizeDegs(1),scSizeDegs(2));

% The ieAddObject followed by sceneWindow bring up a simple
% gui that allows visualization of the scene.
%
% Note that once we create the scene, the display information
% is not associated with it.  Rather the scene stores the image
% data as the spectra at each pixel.
%
% A later tutorial will unpack ISETCam/Bio scenes in more detail.
sceneWindow(sc);

%% Set up oi object for computing the retinal image
%
% The term oi stands for optical image, and here that
% means the retinal image.  
%
% We get some reasonable human optics by default.  These
% are based on published human wavefront optics measurements
% for foveal viewing. We can dial in other optics from various
% databases that are included in ISETBio.  We will discuss how
% to do this, and the various parameters of the oi in a later
% tutorial.
%
% Note that we set the wavelengthe support of the oi to match
% that we used for the scene.  
oi = oiCreate('human');
oi = oiSet(oi,'wave',wave);

% Plot the default foveal polychromatic point spread function at
% some sample wavelengths
%
% Probably we have some way to do this in one simple call, but if
% so I couldn't easily find it.  So just make the plot we want 
% by hand
plotWls = [400 500 550 600 700];
plotLimMinutes = 100;
figure; set(gcf,'Position',[10 10 length(plotWls)*300 300]); clf;
for i = 1:length(plotWls)
    subplot(1,length(plotWls),i);
    uData = oiPlot(oi,'psf','wave',plotWls(i),'nofigure');
    mesh(uData.x,uData.y,uData.psf);
    xlim([-plotLimMinutes plotLimMinutes]);
    ylim([-plotLimMinutes plotLimMinutes]);
    xlabel('Position (min)');
    ylabel('Position (min)');
    title(sprintf('PSF at %d nm',plotWls(i)));
end

%% Computes the retinal image from the scene, using the oi
%
% The retinal image is computed at each wavelength with a polychromatic
% point spread function that takes typical human axial chromatic aberration
% into account, as shown above.  The spectral absorption of light by the
% lens is also taken into account.
%
% At the edge of the scene, we need to make some assumption about what is
% outside the scene. Here we tell the compute method to pad outside the
% scene image by its mean value at each wavelength.  The padding makes the
% retinal image larger than the scene.
oi = oiCompute(oi, sc,'pad value','mean', 'crop', true);
oiSizeDegs = [oiGet(oi,'w angular') oiGet(oi,'h angular')];
fprintf('Scene size is %0.1f by %0.1f degrees\n',oiSizeDegs(1),oiSizeDegs(2));

% As with the scene, we can have a look through a gui.  It looks yellower
% than the scene because the lens aborbs more light at shorter wavelengths.
oiWindow(oi);

%% Build a default cone mosaic and compute excitations
%
% Get the default set of cone mosaic parameters and set the ones we want.
% We match the aspect ratio of the mosaic to the oi, which itself is a bit
% larger than the scene.
%
% We place the mosaic at the fovea, but you can move it around
% with the eccentricity parameters.  The default integration time
% is 5 msec, which is good for dynamic calculations but not so
% good for single images, so we change to 100 msec.
%
% If you do move to some other eccentricity, you would want to select
% optics matched for that eccentricity.  We will unpack how to do this in a
% later tutorial.
cmParams = cMosaicParams;
cmParams.integrationTime = 0.100;
cmParams.eccentricityDegs = [0 0];
cmParams.sizeDegs = oiSizeDegs;
cmParams.micronsPerDegree = oiGet(oi,'distance per degree','um');

% Create the mosaic
cm = cMosaic(cmParams);

% Visualize how the cones are packed
cm.visualize;

%% Compute cone photopigment excitations
%
% The noise free responses are the mean excitations for each cone.  The noisy
% responses at Poisson-distributed noise to the mean excitations.
[noiseFreeExcitationList, noisyExcitationList] = cm.compute(oi);

%% Visualize excitations
%
% The plot command commented out here is simple but
% does not give much flexibility
cm.plot('excitations',noiseFreeExcitationList,'labelcones',true);

% The visualize command has many options but is more
% complex to use.
vParams = cm.visualize('params');
vParams.activation = noisyExcitationList;
vParams.activationColorMap = gray(512);
vParams.verticalActivationColorBar = true;
vParams.activationRange = [0 max(noisyExcitationList(:))];
vParams.labelcones = false;
cm.visualize(vParams);

% That image looks low contrast because the S cones occupy
% a different intensity range than the L and M cones.  We
% can adjust this by setting the range to the central portion
% of the cone excitation response range as shown here.
vParams.activationRange = prctile(noisyExcitationList(:),[5 90]);
cm.visualize(vParams);

% Here are some slices through the acivations
cm.plot('excitations horizontal line',noiseFreeExcitationList, 'y deg',0,'thickness',0.05,'conetype','l');
cm.plot('excitations horizontal line',noiseFreeExcitationList, 'y deg',0,'thickness',0.05,'conetype','m');
cm.plot('excitations horizontal line',noiseFreeExcitationList, 'y deg',0,'thickness',0.05,'conetype','s');

%% Unpack what's in the mosaic a little
%
% You might want a list of the cone positions
% A little sleuthing will be required to decide
% which is horizontal and which is vertical.
%
% The actual excitations in noiseFree and noisy are in a list ordered as the two lists below.
% The first coordinate is the x (horizontal) position of the cone, the
% second is the y.  A later tutorial will unpack the coordinate system used
% for retinal position. As the variable names suggests, these are in
% microns.
conePositionsMicrons = cm.coneRFpositionsMicrons;

% Cone types
%
% LMS coded as 1, 2, and 3 respectively.
coneTypes = cm.coneTypes;

% The continuous LMS images.  As they sit in the
% cMosaic property, they are up down flipped to
% match what happens with real retinal images.  Useful
% flip back.
8;

%% END
