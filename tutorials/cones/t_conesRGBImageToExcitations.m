%% Further introduction to the cone mosaic (cMosaic) object.
%
% Description:
%    Create a cone mosaic object and compute cone isomerizations
%    to an RGB image.
%    No eye movements.
%
%    Visualize the results with a plot.
%
% See also
%   t_cones*

%% Initialize and clear
clear; close all; ieInit;

%% Set sample wavelength spacing
wave = 400:5:700;

%% Create scene from image file
% 
% Here will will build the scene from a stock RGB image, and
% a "typical" display.
%
% Although sceneFromFile will read the scene directly and create
% the image, we sometimes want a higher pixel resolution so we read and
% upsample in advance of calling sceneFromFile. Setting upsampleFactor
% to 1 leaves the pixel resolution alone.
%
% We also explicitly set the width of the field of view in degrees,
% and the mean scene luminance, since neither of these comes in 
% automatically with the RGB image.  Luminance is in cd/m2.
%
% Set parameters
imgFileName = 'zebra.jpg';              % Available images include: 'macbeth.tif', 'bears.png'. 'eagle.jpg', 'hats.jpg', 'lion-cardinal.png'
dispLCDFile = 'LCD-Apple.mat';
upsampleFactor = 1;
fieldOfViewDegs = 1;
sceneLuminance = 200;

% Read and resize the RGB image.
im = imread(imgFileName);
im = imresize(im,upsampleFactor);

% Routine sceneFromFile does the work.  It undoes the
% monitor's gamma correction to linearize the image data,
% and then converts to spectra using the monitor's channel
% spectra.
sc = sceneFromFile(im,'rgb',[],dispLCDFile,wave);

% Then set the scene parameters we want to adjust
sc = sceneSet(sc,'name','Scene on LCD');
sc = sceneSet(sc,'fov',fieldOfViewDegs');
sc = sceneSet(sc,'mean luminance',sceneLuminance);
ieAddObject(sc); sceneWindow;

%% Compute and show the oi
%
% We get some reasonable human optics by default.  We
% can dial in other optics from various databases that
% are included in ISETBio.
oi = oiCreate('human');
oi = oiSet(oi,'wave',wave);
oi = oiCompute(oi, sc,'pad value','mean');
oiSizeDegs = [oiGet(oi,'w angular') oiGet(oi,'h angular')];
ieAddObject(oi); oiWindow;

%% Plot the polychromatic point spread function at sample wavelengths
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

%% Build a default cone mosaic and compute isomerizatoins
%
% Get the default set of cone mosaic parameters and
% set the ones we want. We match the aspect ratio
% of the mosaic to the oi, which itself is a bit
% larger than the scene.
%
% We place the mosaic at the fovea, but you can move it around
% with the eccentricity parameters.  The default integration time
% is 5 msec, which is good for dynamic calculations but not so
% good for single images, so we change to 100 msec.
cmParams = cMosaicParams;
cmParams.integrationTime = 100;
cmParams.eccentricityDegs = [0 0];
cmParams.sizeDegs = oiSizeDegs;
cmParams.micronsPerDegree = oiGet(oi,'distance per degree','um');

% Create the mosaic
cm = cMosaic(cmParams);

% Visualize how the cones are packed
cm.visualize;

%% Compute cone photopigment excitations
[noiseFree, noisy] = cm.compute(oi);

%% Visualize excitations
%
% The plot command commented out here is simple but
% does not give much flexibility
cm.plot('excitations',noiseFree,'labelcones',true);

% The visualize command has many options but is more
% complex to use.
vParams = cm.visualize('params');
vParams.activation = noisy;
vParams.activationColorMap = gray(512);
vParams.verticalActivationColorBar = true;
vParams.activationRange = [0 max(noisy(:))];
vParams.labelcones = false;
cm.visualize(vParams);

% Here are some slices through the acivations
cm.plot('excitations horizontal line',noiseFree, 'y deg',0,'thickness',0.05,'conetype','l');
cm.plot('excitations horizontal line',noiseFree, 'y deg',0,'thickness',0.05,'conetype','m');
cm.plot('excitations horizontal line',noiseFree, 'y deg',0,'thickness',0.05,'conetype','s');

%% Unpack what's in the mosaic a little
%
% You might want a list of the cone positions
% A little sleuthing will be required to decide
% which is horizontal and which is vertical.
%
% The actual excitations noiseFree, noisy are in a list ordered as the two lists below.
conePositionsMicrons = cm.coneRFpositionsMicrons;

% Cone types
%
% LMS coded as 1, 2, and 3
coneTypes = cm.coneTypes;

% The continuous LMS images.  As they sit in the
% cMosaic property, they are up down flipped to
% match what happens with real retinal images.  Useful
% flip back.
LMSImages = cm.absorptionsDensityFullMap;
LMSImages = LMSImages(end:-1:1,:,:);
figure; clf; imshow(LMSImages/max(LMSImages(:)));

%% END
