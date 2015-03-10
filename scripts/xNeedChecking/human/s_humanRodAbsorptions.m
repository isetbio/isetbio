%% Calculate the number of rod absorptions to a uniform scene 
%
% Redraft - not checked yet - Wandell
%
% Hiroshi Horiguchi, 2012.

%%
s_initISET

%% Set variables

% Photopigment optical density
rodpod = 0.05;

RodinnerSegmentDiameter = 2.22; % 15 deg ecc. Curio 1993
meanluminance  = 2060; % cd/m2
pupilsize = 2; % mm 

rodPeakAbsorptance = 0.66; % from Rodieck

%%
scene = sceneCreate('uniform ee');
wave  = sceneGet(scene,'wave');

% Create a file with your primaries here.
% fullpathname = ieSaveSpectralFile(wavelength,data,comment,[fullpathname]);
% Or just load primaries.
% Note that it requires a path to colorTime in vistaproj
d = displayCreate('lcdExample',wave);
primaries = displayGet(d,'spd');

% multiply your primaries by illEnergy
illEnergy = primaries * ones(size(primaries,2),1);

% apply illuminant energy to scene
scene = sceneAdjustIlluminant(scene,illEnergy);
% sceneGet(scene,'mean luminance') % you'll probably get 100 Cd/m2.

% set luminance you desire
scene = sceneSet(scene,'mean luminance', meanluminance);   % Cd/m2
vcAddAndSelectObject(scene);sceneWindow(scene);

%% create an optical image of human eye
oi = oiCreate('human');
optics = opticsCreate('human', pupilsize / 2 / 1000);
oi = oiSet(oi,'optics',optics);


% Calc rod responses
% Note that it requires a path to colorTime in vistaproj
rodabsorbance = ieReadSpectra('rods',wave);
r = coneCreate;
r = coneSet(r,'absorbance',rodabsorbance);
r = coneSet(r,'name','rod');
r = coneSet(r,'peak efficiency',rodPeakAbsorptance);
% r = coneSet(r,'spatial density',[0 1]);

% rods = cm_variableLMSI_PODandLambda(rodabsorbance, rodpod, [], LensTransmittance(wave));
% rods = rods * rodPeakAbsorptance;

% or
% rods = ieReadSpectra('scotopicLuminosity.mat',wave);
% vcNewGraphWin; plot(wave,rods)


%% open an optical image window
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi);
oiWindow;

%%  Now set up the rod sensor parameters a little better

RodArea  = (RodinnerSegmentDiameter./2)^2 * pi();
Rodpixels = sqrt(RodArea);

% Peak sensitivity - includes lens and rod pigment
pixSize = Rodpixels*1e-6;
sensor = sensorCreateIdeal('monochrome');
sensor = sensorSet(sensor,'pixel size keep fill factor',pixSize);

sensor = sensorSet(sensor,'pixel voltageSwing',        300); % just for visualization

sensor = sensorSet(sensor,'autoexposure',0);
sensor = sensorSet(sensor,'exposureTime',1);

sensor = sensorSet(sensor,'filter spectra',coneGet(r,'absorbance'));
sensor = sensorSet(sensor,'filter names',{'wrod'});

sensor = sensorCompute(sensor,oi);
vcAddAndSelectObject(sensor); sensorImageWindow
%% Calculate number of absorptions (electrons) per rod

roi    = sensorROI(sensor,'center');
sensor = sensorSet(sensor,'roi',roi);
elROI  = sensorGet(sensor,'roi electrons');

% mean of electron
fprintf('Mean isomerizations per second per rod %f\n',mean(elROI));

%% End