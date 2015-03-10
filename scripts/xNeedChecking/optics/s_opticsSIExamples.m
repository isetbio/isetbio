% s_opticsSIExamples.m
%
%   Shift-invariant optics examples
%
% This script illustrates how to create an optics structure used for the
% shift-invariant optics calculations. It also illustrates how to place
% that optics structure into the Optics window (optical image window) and
% run the calculation.
%
% Three examples are below.  One creates a simple pillbox (averaging) PSF.
% The second creates a simple sharpening filter.  The third illustrates the
% built-in ability to create a series of Gaussian blur filters with a
% spread that varies with wavelength and with different x- and y-dimension
% spread ratios.
%
% You can use this script as a model for creating the optics structure of a
% shift-invariant lens.  After creating the psf data, we suggest that you
% save the data in a file called SI-yourLensName within the directory
% data\optics.  The file format for the data, illustrated below, is
% very simple.  The save procedure is described in the routine
% ieSaveSIDataFile (see below).  We suggest you use ieSaveSIDataFile to save
% your shift-invariation point spread data. 
% 
% You can load those data into the optical image window directly (Optics |
% Load SI Data). Or, you can use this script to create the optics and save
% those from the window using Optics | Export, or you can simply save the
% optics using vcSaveObject.
%
% Copyright ImagEval Consultants, LLC, 2008


%% Example 1: Create a shift-invariant optics structure with custom PSF data.
%
% The structure is created using a set of point spread functions.  The
% point spread functions are a set of matrices defined on a grid.  THere is
% one point spread for each wavelength.  In this example, there are 31
% wavelengths (400:10:700).  The point spread functions are simply random
% numbers.  The grid is 128 x 128 with samples spaced every 0.25 microns.
% Hence, the total grid size is 32 microns on a side.
%

%%
s_initISET

%% Create a simple test scene

% Let's work with a small checkerboard scene
pixPerCheck = 16;
nChecks = 6; 
scene = sceneCreate('checkerboard',pixPerCheck,nChecks);
wave  = sceneGet(scene,'wave');
scene = sceneSet(scene,'fov',3);

% Replace the optical image into your ISET window
vcAddAndSelectObject(scene);
sceneWindow

%% Create Shift-invariant data

% Now, write out a file containing the relevant point spread function
% data, along with related variables.
umPerSample = [0.25,0.25];                % Sample spacing

% Point spread is a little square in the middle of the image
h = zeros(128,128); h(48:79,48:79) = 1; h = h/sum(h(:));
for ii=1:length(wave), psf(:,:,ii) = h; end     % PSF data

% Save the data
ieSaveSIDataFile(psf,wave,umPerSample,'SI-pillBox');

%% Read the custom data and put it into an optics structure.
oi = oiCreate;
optics = siSynthetic('custom',oi,'SI-pillBox',[]);

% Make sure the program knows you want to use shift invariant
optics = opticsSet(optics,'model','shiftInvariant');

% Attach the optics structure to the optical image structure
oi = oiSet(oi,'optics',optics);

% You can now compute using your current scene.
oi = oiCompute(scene,oi);
oi = oiSet(oi,'name','Pillbox');

% Show the OI window
vcAddAndSelectObject(oi);
oiWindow;

% Use Analyze | Optics | XXX to plot various functions in the optics
% (optical image) window.


%% Example 2: Create a a slight sharpening filter. 
h1 = fspecial('gaussian', 128, 5);
h2 = fspecial('gaussian', 128, 10);
h = h1 - 0.5*h2;
% figure(1); mesh(h)

psf = zeros(128,128,length(wave));
for ii=1:length(wave), psf(:,:,ii) = h; end     % PSF data

% Save the data and all the rest, in compact form 
ieSaveSIDataFile(psf,wave,umPerSample,'customFile');

optics = siSynthetic('custom',oi,'customFile','deleteMe');
optics = opticsSet(optics,'model','shiftInvariant');
oi     = oiSet(oi,'optics',optics);

oi = oiCompute(scene,oi);
oi = oiSet(oi,'name','Sharpened');
vcAddAndSelectObject(oi);
oiWindow;

%% Example 3: Create a wavelength-varying shift-invariant gaussian PSF

wave    = oiGet(oi,'wave');
psfType = 'gaussian';

% The special set up for the wavelength-dependent Gaussian spread functions
% is illustrated here.

%We set the spread as a function of wavelength.  This function is
%arbitrary - it makes the long wavelength much worse than the short.
waveSpread = (wave/wave(1)).^3;

% Make point spreads with a symmetric Gaussian
xyRatio = ones(1,length(wave));

% Now call the routine with these parameters
optics  = siSynthetic(psfType,oi,waveSpread,xyRatio,'gaussPSF');

% Here is the rest of the computation, as above
optics  = opticsSet(optics,'model','shiftInvariant');
oi      = oiSet(oi,'optics',optics);
scene   = vcGetObject('scene');
oi      = oiCompute(scene,oi);

oi = oiSet(oi,'name','Chromatic Gaussian');
vcAddAndSelectObject(oi);
oiWindow;

% We saved these shift invariant optics using:
%
% fullName = fullfile(isetRootPath,'data','optics','si2x1GaussianWaveVarying');
% (or) fullName = fullfile(isetRootPath,'data','optics','siGaussianWaveVarying');
% vcExportObject(optics,fullName)
% You can import these optics using:  Optics | Import Optics

% Use Analyze | Optics | PSF Mesh to see the point spread a different
% wavelengths

% Use Analyze | Line x Wave Plot to see the loss of contrast a long
% wavelengths compared to short wavelengths

%% Loading (2x1) asymmetric optics data from a file

fullName = fullfile(isetRootPath,'data','optics','si2x1GaussianWaveVarying.mat');
newVal   = vcImportObject('OPTICS',fullName);
oi       = vcGetObject('oi');
oi       = oiCompute(scene,oi);

oi = oiSet(oi,'name','Asymmetric Gaussian');
vcAddAndSelectObject(oi);
oiWindow;

%%
imageMultiview('oi',1:4,1)

%% End
