% Script s_opticsRTPSFandFigs
%
% This script demonstrates the complete set of calculations in the ray trace
% code, including 
%
%    * geometric distortion
%    * relative illumination
%    * field-height dependent PSF blurring.  
% 
% This script takes a while to compute
%
% The script brings up an array of figures that show different aspects of
% the calculation.  These are useful for teaching and explaining the
% calculations.
%
% The scene is transformed to an optical image using ray trace methods
% based on the aspherical, 2mm lens computed in Zemax.
%
% The scene is also transformed using diffraction limited methods
% (shift-invariant).  The f# and focal length of the diffraction model are
% set equal to those of the ray trace lens.
%
% The illuminance computed the two ways is then compared.
%
% Copyright ImagEval, LLC, 2005


%%
s_initISET

%% Set up a wide angle scene for the wide angle lens below
scene = sceneCreate('gridlines',[384 384], 48);
scene = sceneInterpolateW(scene,(550:100:650));  % Small wavelength sample
scene = sceneSet(scene,'hfov',45);
scene = sceneSet(scene,'name','rtDemo-Large-grid');

% Show the grid line scene
vcAddAndSelectObject(scene); sceneWindow;

%% Import wide angle lens optics.

% Set up default optical image
oi = oiCreate;

% Load in the wide angle lens optics file created by Zemax (zm)
opticsFileName = fullfile(isetRootPath,'data','optics','zmWideAngle.mat');
tmp = load(opticsFileName);

% Set the oi with the optics loaded from the file
oi = oiSet(oi,'optics',tmp.optics);

% Retrieve it and print its name to verify and inform user
optics = oiGet(oi,'optics');
fprintf('Ray trace optics: %s\n',opticsGet(optics,'lensFile'));

%% Set up diffraction limited parameters to match the ray trace numbers

% Now, match the scene properties
oi = oiSet(oi,'wangular',sceneGet(scene,'wangular'));
oi = oiSet(oi,'spectrum',sceneGet(scene,'spectrum'));

% Match the scene distance and the rt distance.  They are both essentially
% infinite.
scene   = sceneSet(scene,'distance',2);  % Two meters - essentially inf
optics  = opticsSet(optics,'rtObjectDistance',sceneGet(scene,'distance','mm'));   

%% Compute the ray trace distortion and show it in the OI

% We calculate in the order of (a) Distortion, (b) Relative
% illumination, and then (c) OTF blurring The function rtGeometry
% calculates the distortion and relative illumination at the same time.
irradiance = rtGeometry(scene,optics);

% Copy the resulting data into the optical image structure
oi = oiSet(oi,'cphotons',irradiance);
oi = oiSet(oi,'name','Geometry only');
vcAddAndSelectObject('oi',oi); oiWindow;

%% Precompute the PSF and then apply
%
angStep = 20;   % Very coarse for speed
svPSF    = rtPrecomputePSF(oi,angStep);
oi = oiSet(oi,'psfStruct',svPSF);

% Apply
outIrrad = rtPrecomputePSFApply(oi,angStep);
oi = oiSet(oi,'cphotons',outIrrad);
oi = oiSet(oi,'name','Stepwise-RT');
vcAddAndSelectObject('oi',oi); oiWindow;

%% We choose ray trace by setting the optics method
%
oi = oiSet(oi,'optics model','ray trace');
oi = oiCompute(scene,oi);
oi = oiSet(oi,'name','Automated ray trace');

% Have a look - barrell distortion and all
vcAddAndSelectObject('oi',oi); oiWindow;

% Here is a horizontal line of illuminance
rtData = plotOI(oi,'illuminance hline',[1,64]);

%% Compute using the diffraction-limited method 
%
oiDL = oiSet(oi,'optics model','diffraction limited');
optics = oiGet(oiDL,'optics');

% Set the diffraction limited f# from the ray trace values
fNumber = opticsGet(optics,'rt fnumber');
optics = opticsSet(optics,'fnumber',fNumber);
oiDL = oiSet(oiDL,'optics',optics);

% Now set the method to diffraction limited and compute
oiDL = oiSet(oiDL,'name','DL method');
oiDL = oiCompute(scene,oiDL);

% No barrel distortion, less blurring
vcAddAndSelectObject('oi',oiDL); oiWindow;

% Here is a horizontal line of illuminance
dlData = plotOI(oiDL,'illuminance hline',[1,64]);

%% Make the FOV smaller to see the ray trace blurring
%
% The first calculation was spatially coarse, and inappropriate for a
% geometric calculation such as barrel distortion. 
%
% Here, we make the scene smaller and recalculate. With this field of view
% there is no noticeable distortion, but the sample spacing is much finer
% so we can see the various point spread functions
%
% At this resolution, the calculation takes a little while.

sceneSmall = sceneSet(scene,'name','rt-Small-Grid');
sceneSmall = sceneSet(sceneSmall,'fov',20);

% Ray trace calculation with distortion and shift-variant blurring
oi = oiCompute(sceneSmall,oi);
oi = oiSet(oi,'name','rt-Small-RT');
vcAddAndSelectObject('oi',oi); oiWindow;

%% Equivalent diffraction limited
% There is no distortion computed and the scene is small
% This calculation is pretty quick.
oiDL = oiCompute(sceneSmall,oiDL);
oiDL = oiSet(oiDL,'name','rt-Small-DL');
vcAddAndSelectObject('oi',oiDL); oiWindow;

%% End

