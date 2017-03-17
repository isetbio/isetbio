% t_blenderSceneResponse
% 
% Take a hyperspectral image generated from a Blender 3D model and find the
% response in isetbio.
% 
% 3/2017 JRG/HB/BW (c) isetbio team

%% Initialize parameters

% Initialize parameters of simulated retinal patch
ecc = [0,0]*1e-3;   % Cone mosaic eccentricity in meters from fovea
fov = 2;            % Scene Field of view in degrees
emLength = 40;     % Eye movement frames

cellType = 'on parasol';
% cellType = 'off midget';

osFlag  = 0;        % 0 for osLinear, 1 for osBioPhys
%% Save hyperspectral scene to RDT

% rd = RdtClient('isetbio');
% rd.credentialsDialog;
% rd.crp('/resources/data/istim');
% fname = [isetbioRootPath 'local/City_pinhole_19.mat'];
% % % rd.publishArtifact(fname);

%% Load scene from the RDT

rdt = RdtClient('isetbio');
rdt.crp('/resources/data/istim');
data = rdt.readArtifact('City_pinhole_7', 'type', 'mat');
% data = rdt.readArtifact('City_pinhole_19', 'type', 'mat');

%% Build a scene and oi for computing

s = sceneCreate('multispectral');

s = sceneSet(s,'photons',data.multispectralImage);
s = sceneSet(s,'illuminant',[400:10:700]);
s = sceneSet(s,'fov',fov);
% s = sceneAdjustLuminance(s,10);
vcAddObject(s);

oi = oiCreate;
oi = oiCompute(oi,s);
vcAddObject(oi); % oiWindow;

%% Build a default cone mosaic and compute the OI

if osFlag % osBioPhys
    osCM = osBioPhys(); 
    cMosaic = coneMosaic('center',[0 0]*1e-3,'os',osCM);  % Create the object
else      % osLinear
    cMosaic = coneMosaic('center',[0 0]*1e-3);  % Create the object
end


% % Set cone mosaic size

% cMosaic.rows = 50; cMosaic.cols = 60;
cMosaic.rows = 100; cMosaic.cols = 120;
% cMosaic.rows = 144; cMosaic.cols = 176;
cMosaic.emGenSequence(emLength);

cMosaic.compute(oi);
cMosaic.computeCurrent();

% Show the window
% cMosaic.window;

% Examine the outer segment current
% cMosaic.plot('movie absorptions','vname','deleteme.avi','step',5);

%% Compute the bipolar response

bp = bipolar(cMosaic);
bp.set('sRFcenter',10);
bp.set('sRFsurround',0);
bp.compute(cMosaic);

%% Set other RGC mosaic parameters

clear params innerRetinaSU
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = sqrt(sum(ecc.^2)); 
% params.fov = fov;
params.eyeAngle = 0; ntrials = 0;

% Create RGC object
innerRetina = ir(bp, params);
innerRetina.mosaicCreate('type',cellType,'model','GLM');

nTrials = 1; innerRetina = irSet(innerRetina,'numberTrials',nTrials);

%% Compute the inner retina response

innerRetina = irCompute(innerRetina, bp); 

%% Make the PSTH movie
innerRetina.mosaic{1}.set('dt',1);
lastTime = innerRetina.mosaic{1}.get('last spike time');
psth = innerRetina.mosaic{1}.get('psth');

clear params

% params.vname = fullfile(isetbioRootPath,'local','vernier.avi'); 
param.FrameRate = 5; params.step = 2; params.show = true;

% % View movie of RGC linear response
%  vcNewGraphWin; ieMovie(innerRetinaSU.mosaic{1}.responseLinear);

% View movie of PSTH for mosaic
steadyStateFrame = 40;
% vcNewGraphWin; ieMovie(psth(:,:,steadyStateFrame:end),params);

% % View average of PSTH movie
vcNewGraphWin; 
subplot(121);
oiShowImage(oi);
subplot(122);
imagesc(mean(psth,3)); axis image