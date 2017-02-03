%% Store iStim data in the RDT
%
% These stimulus files are used for testing RGC responses. The purpose is to
% illustrate the firing of the RGCs in response to stimuli with obvious
% spatial patterns. The movies generated from the RGC linear and PSTH
% responses ought to closely resemble the spatial activity of these
% stimuli.
% 
% Here, we build the stimulus, optical image and cone mosaic for the
% test stimuli. We precompute them here and store them with the RDT to save
% time.
% 
% The user must set publishFlag to 1 in order to publish the data to the
% RDT. Default is do not publish.
% 
% Stimuli:
%   * Moving bar
%   * Drifting Gabor
%   * Rings and rays (with eye movements)
%   * Vernier (with eye movements)
%   % Letter (with eye movements)
% 
% See t_rgc*.m for the individual test scripts.
%
% JRG/BW ISETBIO Team, 2016

%% Open the RDT
rd = RdtClient('isetbio');
rd.credentialsDialog;
rd.crp('/resources/data/istim');

publishFlag = 1;

%% Create the moving bar iStim for RDT upload

clear params iStim
params.radius = 1e-3;
params.barWidth = 10; 
params.fov      = 0.6;
params.startFrames = 50;
params.endFrames = 70;
params.integrationTime = 0.001;
params.os = 'linear';
% params.os = 'hex';
iStim = ieStimulusBar(params);  % Full params are returned in iStim
%%
fname = fullfile(isetbioRootPath,'local','barMovie_osLinear.mat');
save(fname,'iStim');

% Publish to RDT
if publishFlag; rd.publishArtifact(fname); end;

%% Create the file with the Gabor movie iStim in it

clear params iStim
params.radius = 1e-3;
params.startFrames = 50;
params.endFrames = 70;
params.freq = 3;
params.nSteps = 350;
params.GaborFlag = 0.3;

params.os = 'biophys';
params.integrationTime = 0.001;
iStim = ieStimulusGabor(params); %#ok<NASGU>

fname = fullfile(isetbioRootPath,'local','gaborMovie_osBioPhys.mat');
save(fname,'iStim');

% Publish to RDT
if publishFlag; rd.publishArtifact(fname); end;

%% Create the rings-rays stimulus with eye movements

% Initialize parameters of simulated retinal patch
ecc = [0,0]*1e-3;   % Cone mosaic eccentricity in meters from fovea
fov = 2;            % Scene Field of view in degrees
emLength = 250;     % Eye movement frames

sceneType = 'rings rays';
cellType = 'on parasol';

osFlag  = 0;        % 0 for osLinear, 1 for osBioPhys

% Build a scene and oi for computing

s = sceneCreate(sceneType);
s = sceneSet(s,'fov',fov);
s = sceneAdjustLuminance(s,10);
vcAddObject(s);

oi = oiCreate;
oi = oiCompute(oi,s);
vcAddObject(oi); % oiWindow;

% Build a default cone mosaic and compute the OI

if osFlag % osBioPhys
    osCM = osBioPhys(); 
    cMosaic = coneMosaic('center',[0 0]*1e-3,'os',osCM);  % Create the object
else      % osLinear
    cMosaic = coneMosaic('center',[0 0]*1e-3);  % Create the object
end

%  Set cone mosaic size
% cMosaic.rows = 100; cMosaic.cols = 120;
cMosaic.rows = 144; cMosaic.cols = 176;
cMosaic.emGenSequence(emLength);

cMosaic.compute(oi,'currentFlag',true);


fname = fullfile(isetbioRootPath,'local','ringsRays_cMosaic.mat');
save(fname,'cMosaic');

% Publish to RDT
if publishFlag; rd.publishArtifact(fname); end;

%% Create the Vernier stimulus cMosaic

% Initialize parameters of simulated retinal patch
emLength = 250;     % Eye movement frames

% Create the display

% Create a display with a linear gamma table, though
% in general it could be the default or anything.
dpi = 300; d = displayCreate('LCD-Apple','dpi',dpi);
viewDist = 1; % viewing distance in meters
d = displaySet(d, 'viewing distance', viewDist);
d = displaySet(d, 'gamma', 'linear');

clear p;
p.display = d;
p.sceneSz = [128,129]; p.barWidth = 2; p.offset = 3; p.lineSpace = 2;
p.meanLum = 100;p.barColor = [.9 0.9 0.9]; p.bgColor = .3;
s = sceneCreate('vernier','display',p);
% vcAddObject(s); sceneWindow;

oi = oiCreate;
oi = oiCompute(oi,s);
% vcAddObject(oi); oiWindow;

% Build a default cone mosaic and compute the OI

cMosaic = coneMosaic;  % Create the object
cMosaic.emGenSequence(emLength);
cMosaic.setSizeToFOV(sceneGet(s,'fov'),...
    'sceneDist',sceneGet(s,'distance'),...
    'focallength',oiGet(oi,'optics focal length'));
cMosaic.compute(oi);
cMosaic.computeCurrent();

fname = fullfile(isetbioRootPath,'local','vernier_cMosaic.mat');
save(fname,'cMosaic');

% Publish to RDT
if publishFlag; rd.publishArtifact(fname); end;

%% Create the isetbio letter stimulus with eye movements

% Initialize parameters of simulated retinal patch
emLength = 250;     % Eye movement frames

cellType = 'on parasol';

% Create the display

% Create a display with a linear gamma table, though
% in general it could be the default or anything.
dpi = 500; d = displayCreate('LCD-Apple','dpi',dpi);
viewDist = 2; % viewing distance in meters
d = displaySet(d, 'viewing distance', viewDist);
d = displaySet(d, 'gamma', 'linear');

font = fontCreate; 
s = sceneCreate('letter', font, d);
% vcAddObject(s); sceneWindow;

oi = oiCreate;
oi = oiCompute(oi,s);
% vcAddObject(oi); oiWindow;

% Build a default cone mosaic and compute the OI
cMosaic = coneMosaic;  % Create the object
cMosaic.emGenSequence(emLength);
cMosaic.setSizeToFOV(sceneGet(s,'fov'),...
    'sceneDist',sceneGet(s,'distance'),...
    'focallength',oiGet(oi,'optics focal length'));
cMosaic.compute(oi);
cMosaic.computeCurrent();

fname = fullfile(isetbioRootPath,'local','letter_cMosaic.mat');
save(fname,'cMosaic');

% Publish to RDT
if publishFlag; rd.publishArtifact(fname); end

%%
