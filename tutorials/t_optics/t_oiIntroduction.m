%T_OIINTRODUCTION Introduction to ISET objects and functions: The optical image (OI)
%
% Description:
%    ISET is organized around several key objects (scene, optics, sensor
%    and processor) that reflect the important parts of the image
%    acquisition pipeline. The ISET methods generally rely on calculating
%    with these objects. By using these objects properly, the code and
%    analysis are much easier to understand.
%
%    This tutorial introduces the ISET object "optical image" or OI. The OI
%    describes the irradiance field at the sensor surface. It is computed
%    from the scene (radiance), and the optics. The optics object is stored
%    in a slot within the optical image.
%
% Notes:
%    * Use t_"TAB KEY" to see the list of tutorials 
%

% History:
%    xx/xx/12       Copyright ImagEval Consultants, LLC, 2012.
%    09/10/18  jnm  Formatting

%% Initialize ISET 
ieInit;

%% Create a point array scene
% This and other simple synthetic scenes are useful for examining simple
% optical properties, this will be accomplished by creating an array of
% points, and then setting the FoV to 1 degree.
scene = sceneCreate('pointarray');
scene = sceneSet(scene, 'hfov', 1);
sceneWindow(scene);

%% Create and show an optical image
% ISET has several optics models that you can experiment with. These
% include shift-invariant optics, in which there is a different
% shift-invariant pointspread function for each wavelength, and a ray-trace
% method, in which we read in data from Zemax and create a shift-variant
% set of pointspread functions along with a geometric distortion function.
% The default optical image has diffraction limited optics. We compute the
% irradiance field from the scene.
oi = oiCreate('diffraction');
oi = oiCompute(scene, oi);
oi = oiSet(oi, 'name', 'Small f number');
oiWindow(oi);

%% Increase the f-number of the optics and compute again. 
% Increasing the f-number of the optics will blur the irradiance image. The
% optics structure is attached to the optical image. To increase the
% f-number we get the optics object, set its' f number, and then reattach
% it to the optical image.
fnSmall = oiGet(oi, 'optics fNumber');
fnBig = 3 * fnSmall;
oiBigF = oiSet(oi, 'optics fNumber', fnBig);
oiBigF = oiCompute(scene, oiBigF);

oiBigF = oiSet(oiBigF, 'name', 'Big f number');  % Name for the GUI
oiWindow(oiBigF);

%% Use oiPlot to compare the two different optics 
% For this example, plot the PSF at 600nm. First for the small f/#.
thisWave = 600;

% Most ISET plotting routines can return the plotted data into a variable
pData = oiPlot(oi, 'psf', [], thisWave);
set(gca, 'xlim', [-20 20], 'ylim', [-20 20])
colormap(0.5 * gray + 0.5 * ones(size(gray)))
title(sprintf('F-number = %d', fnSmall))

% Here are the values returned by oiPlot
pData

%% Show the effect of the larger f/# on the PSF
% This is for the larger f/#. We just show a 20 um square region.
oiPlot(oiBigF, 'psf', [], thisWave);
set(gca, 'xlim', [-20 20], 'ylim', [-20 20])
colormap(0.5 * copper + 0.5 * ones(size(copper)))
title(sprintf('F-number = %d', fnBig))

% We can also get the plotted values from the figure, shown here.
pData = get(gcf, 'userdata')

%% Use oiPlot to quantify the optical image 
% This example plots an RGB representation of the irradiance field in the
% optical image.

% For this plot, the xy variable is irrelevant, therefore it is entered
% below as []. The grid spacing, however, is set at 17 microns.
gSpacing = 17;
oiPlot(oi, 'irradiance image with grid', [], gSpacing);
title(sprintf('F-number = %d', fnSmall))

%% Plot an RGB representation of the new irradiance field
% The xy variable is irrelevant, and so entered below as [], however, the
% grid spacing is crucial to the representation, which is recorded here as
% 17 microns.
gSpacing = 17;
oiPlot(oiBigF, 'irradiance image with grid', [], gSpacing);
title(sprintf('Zoomed view: F-number = %d', fnBig))

%% Now, run using parameters from human optics
% The defcous of the human optics is quite surprising. In particular, the
% blurring of the short wavelength light is very extreme. That is
% illustrated in this example.
oiHuman = oiCreate('human');
oiHuman = oiCompute(scene, oiHuman);
oiHuman = oiSet(oiHuman, 'name', 'human optics');
oiWindow(oiHuman);

%% End
