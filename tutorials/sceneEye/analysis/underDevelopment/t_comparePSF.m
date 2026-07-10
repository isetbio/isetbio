% SkipFile
%% t_comparePSF.m
%
% Visually compare the PSF generated from the 3D eye modeling and the PSF
% generated from ISETBIO's wavefront tools. Let's also compare how the
% PSF's change with defocus/accommodation.
%
% Depends on: iset3d, isetbio, Docker, RemoteDataToolbox
%
% TL ISETBIO Team, 2017

%% Initialize ISETBIO
ieInit;

%% Generate a PSF using ISETBIO
% To parallel the way we generate the PSF in the next section, we will not
% generate the PSF from the OTF, but instead pass ISETBIO a scene
% consisting of a single lit pixel in an image and calculate the response
% through ISETBIO optics.

% Create an equal energy display
d = displayCreate('equal energy');

pFile = fullfile(piRootPath,'data','imageTextures','pointTest.png');
scene = sceneFromFile(pFile, 'rgb', [], d);

% Set equal energy illumination
wave = sceneGet(scene, 'wave');
onesPhotons = ones(size(wave)) * 1e+15;
equalPhotonsEnergy = Quanta2Energy(wave, onesPhotons);
scene = sceneAdjustIlluminant(scene, equalPhotonsEnergy);

% Set the FOV
horFieldofView = 1;
scene = sceneSet(scene,'fov',horFieldofView);

% Place scene at infinity
scene = sceneSet(scene, 'distance', 100000001);

% Get size and parameters of the lit pixel
pixelWidth_m = sceneGet(scene,'sample size');
sceneDistance_m = sceneGet(scene,'distance');
sceneFOV = sceneGet(scene,'fov');
sceneSize_m = sceneGet(scene,'height and width');

% Show the scene
ieAddObject(scene); sceneWindow;

% Create human optics
oi_2d = oiCreate; % By default, ISETBIO uses Marimont and Wandell

% ---- WVF -----
maxUM = 30;    
measPupilMM = 4.5;    % This selects which Thibos data set to load
calcPupilMM = 4.0;    % Calculate for this pupil size

[sample_mean ~] = wvfLoadThibosVirtualEyes(measPupilMM);

% Allocate space and fill in the lower order Zernicke coefficients
z = zeros(65,1);
z(1:13) = sample_mean(1:13);

% Create the example subject
thisGuy = wvfCreate;                                  % Initialize
thisGuy = wvfSet(thisGuy,'zcoeffs',z);                % Zernike
thisGuy = wvfSet(thisGuy,'measured pupil',measPupilMM);   % Data
thisGuy = wvfSet(thisGuy,'calculated pupil',calcPupilMM); % What we calculate
thisGuy = wvfSet(thisGuy,'measured wavelength',550);
thisGuy = wvfSet(thisGuy,'calc wave',[450:100:750]');     % Must be a column vector 
thisGuy = wvfComputePSF(thisGuy);

oi_2d = wvf2oi(thisGuy);

% ---------------

% Remove lens transmittance
wave = oi_2d.optics.lens.get('wave');
oi_2d.optics.lens.set('unitdensity',ones(size(wave)));

oi_2d = oiSet(oi_2d,'name','psf_2d');
oi_2d = oiCompute(oi_2d,scene);

% Get the oi parameters we need
oiRes = oiGet(oi_2d,'rows');
oiFOV = oiGet(oi_2d,'fov');
oiSS = oiGet(oi_2d,'sample size');

% Show the optical image
ieAddObject(oi_2d);
oiWindow;

% Plot the PSF
psfData = oiPlot(oi_2d,'psf 550');

% figure(1);
% x = psfData.x(1,:);
% psfXLine = psfData.psf(30,:);
% plot(x,psfXLine);
% xlabel('Position (um)');
% grid on;

% Plot oi cross section
% irradData = oiPlot(oi_2d,'illuminancehline','roiLocs',[320 320]);
%
% figure(2);
% x = irradData.pos;
% irradXLine = irradData.data(16,:); % 550 nm
% plot(x,irradXLine);
% xlabel('Position (um)');
% grid on;


%% Generate a PSF at focus using the eye model

% Load scene
% A black disk located 100 meters away (on the y-axis) and spanning 120 deg FOV.
scene3d = sceneEye('blackBackdrop');

% Calculate the FOV of the lit pixel
pixelFOV = 2*atand((pixelWidth_m/2)/sceneDistance_m);

% Calculate size of sphere given an arbitrary sphere distance. For some
% reason, using to large a pixel distance breaks things. Maybe overflow
% in the C++ code?
pointDistance = 100; % In meters
pointWidth = 2*pointDistance*tand(pixelFOV/2);

% Add equal energy sphere on the optical axis. We will make it the same size as the
% pixel above.
scene3d.recipe = piAddSphere(scene3d.recipe,...
    'spectrum',[400 1 800 1],...
    'radius',pointWidth/2,...
    'location',[0 pointDistance 0]);

% Set rendering parameters
% Drop FOV to have higher chance of hitting the point
fovScale = 0.1;
scene3d.fov = oiFOV*fovScale;

% Accommodate to the point
scene3d.accommodation = 1/pointDistance;

% Calculate new resolution
scene3d.resolution = round(scene3d.width/(oiSS*10^3));
scene3d.numRays = 2^14; % This needs to be high, since the point is small and hard to hit!  
% scene3d.numRays = 256;

% Add chromatic aberration
scene3d.numCABands = 16;

% Loop over pupil diameters
for pd = [4] % 3]
    
    % Change pupil diameter
    scene3d.pupilDiameter = pd;
    
    % Render through eye
    scene3d.name = sprintf('psf_3deye_%imm',pd);
    [oi_3d, result] = scene3d.render;
    
    % Show the optical image
    ieAddObject(oi_3d);
    oiWindow;
    
    % Save it!
    saveName = sprintf('%s.mat',scene3d.name);
    oi = oi_3d; % Rename for easier access
    save(saveName,'oi','scene3d');
    
end

% Do some calculations
sensorPointWidth_mm = 2*tand(pixelFOV/2)*16.32;
fprintf('ANALYSIS: \n')
fprintf('Pixel FOV = %f deg \n',pixelFOV)
fprintf('For distance %f m, the point width is %f mm. \n',pointDistance,pointWidth*10^3);
fprintf(['Assume a pinhole system, the sphere size on the sensor '...
    'should be around: %0.2f um \n'],sensorPointWidth_mm*10^3);

sensorPixelWidth = oiGet(oi_3d,'sample size');
sensorPixelWidth = sensorPixelWidth(1);
fprintf('Each pixel on the sensor has width: %0.2f um \n',sensorPixelWidth*10^6)

% Crop the oi_2d so we can make a comparison
cropSize = oiGet(oi_3d,'rows');
cropSizeHalf = round(cropSize/2);
oiCenter = oiGet(oi_2d,'rows')/2;
oi_2d_crop = oiCrop(oi_2d,[oiCenter-cropSizeHalf oiCenter-cropSizeHalf cropSize cropSize]);
ieAddObject(oi_2d_crop);
oiWindow;
