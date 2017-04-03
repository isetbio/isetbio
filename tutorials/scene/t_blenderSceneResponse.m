% t_blenderSceneResponse
% 
% Take a hyperspectral image generated from a Blender 3D model and find the
% response in isetbio.
% 
% 3/2017 JRG/HB/BW (c) isetbio team

%% Initialize parameters

% Initialize parameters of simulated retinal patch
ecc = [0,0]*1e-3;    % Cone mosaic eccentricity in meters from fovea
fov = 2;             % Scene Field of view in degrees
emLength = 20;      % Eye movement frames

cellType = 'on parasol';
% cellType = 'off midget';

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

s = sceneCreate();

% Hack to make reflectance about the right range
sFactor = 2.5;
ill = illuminantCreate('d65');
mx = max(data.multispectralImage(:));
data.multispectralImage = (data.multispectralImage/mx)*max(ill.data.photons(:))*sFactor;

s = sceneSet(s,'photons',data.multispectralImage);
s = sceneSet(s,'fov',fov);
s = sceneSet(s,'illuminant',ill);

s = sceneAdjustLuminance(s,100);
vcAddObject(s); sceneWindow;

%%
oi = oiCreate;
oi = oiCompute(oi,s);
vcAddObject(oi); % oiWindow;

%% Build a default cone mosaic and compute the OI
%
% if osFlag % osBioPhys
%     osCM = osBioPhys();
%     cMosaic = coneMosaic('center',[0 0]*1e-3,'os',osCM);  % Create the object
% else      % osLinear
%         cMosaic = coneMosaic('center',[0 0]*1e-3);  % Create the object
% end

cMosaic = coneMosaic('center',[0 0]*1e-3);  % Create the object
cMosaic.integrationTime = 0.005;            % Five ms integration time
cMosaic.setSizeToFOV(fov*1.1);              % A little bigger for eye move
cMosaic.emGenSequence(emLength);

% Compute and show
cMosaic.compute(oi);
%
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
bp.window;


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
innerRetina.mosaic{1}.window;

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

%%