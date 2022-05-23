%% Use the coneMosaic object to simulate responses of foveal RGC mosaic
%
% This tutorial generates RGC responses to the isetbio vernier
% stimulus with eye movements.
%
%   * Create a scene
%   * Create an optical image
%   * Calculate a cone mosaic of the scene with eye movements 
%   * Calculate bipolar
%   * Calculate RGC for on parasol 
%
% Based on t_coneMosaic and t_VernierClassifier.
% 
% 7/2016 JRG HJ BW (c) isetbio team

%% Initialize parameters

clx; ieInit;

% Initialize parameters of simulated retinal patch
ecc = [0,0]*1e-3;   % Cone mosaic eccentricity in meters from fovea
fov = 1;            % Scene Field of view in degrees
emLength = 500;     % Eye movement frames

cellType = 'on parasol';

%% Get iStim structure for vernier movie from RDT
% rdt = RdtClient('isetbio');
% rdt.crp('/resources/data/istim');
% 
% switch osFlag
%     case 0 % osLinear
%         data = rdt.readArtifact('vernier_cMosaic', 'type', 'mat');
%     case 1 % osBioPhys
%         data = rdt.readArtifact('vernier_cMosaic_osBioPhys', 'type', 'mat');
% end
% 
% % iStim = data.iStim; clear data;
% cMosaic = data.cMosaic;

%% Create the display

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

%%

oi = oiCreate;
oi = oiCompute(oi,s);
% vcAddObject(oi); oiWindow;

%% Build a default cone mosaic and compute the OI

cMosaic = coneMosaic;  % Create the object
% cMosaic = coneMosaic('center',[0 0]*1e-3);  % Create the object
% cMosaic = coneMosaicHex(5,true);  % Create the object
cMosaic.emGenSequence(emLength);
cMosaic.setSizeToFOV(sceneGet(s,'fov'),...
    'sceneDist',sceneGet(s,'distance'),...
    'focallength',oiGet(oi,'optics focal length'));
cMosaic.compute(oi);
cMosaic.computeCurrent;

% Show the window
% cMosaic.window;

% Examine the outer segment current
% cMosaic.plot('movie absorptions','vname','deleteme.avi','step',5);

%% Compute the bipolar response

bp = bipolar(cMosaic);
bp.set('sRFcenter',1);
bp.set('sRFsurround',0);
bp.compute(cMosaic);
% bp.plot('movie response')

%% Set other RGC mosaic parameters

clear params innerRetinaSU
params.name = 'vernier test';
params.eyeSide = 'left'; 
params.eyeRadius = sqrt(sum(ecc.^2)); 
% params.fov = fov;
params.eyeAngle = 0; ntrials = 0;

% Create RGC object
innerRetinaSU = ir(bp, params);
innerRetinaSU.mosaicCreate('type',cellType,'model','GLM');

nTrials = 1; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);

%% Compute the inner retina response

innerRetinaSU = irCompute(innerRetinaSU, bp); 
lastTime = innerRetinaSU.mosaic{1}.get('last spike time');

%% Make the PSTH movie
innerRetinaSU.mosaic{1}.set('dt',1);
psth = innerRetinaSU.mosaic{1}.get('psth');

clear params
% params.vname = fullfile(isetbioRootPath,'local','vernier.avi'); 
param.FrameRate = 3; params.step = 1; params.show = true;

% % View movie of RGC linear response
%  vcNewGraphWin; ieMovie(innerRetinaSU.mosaic{1}.responseLinear);

% Show the movie after the signal has reached steady state
steadyStateFrame = 30;  % In ms if dt is 1

% View movie of PSTH for mosaic
vcNewGraphWin; ieMovie(psth(:,:,steadyStateFrame:end),params);

% % View average of PSTH movie
vcNewGraphWin; 
subplot(121);
oiShowImage(oi);
subplot(122);
imagesc(mean(psth,3)); axis image

% % Plots of RGC linear response and OS current
% vcNewGraphWin; plot(RGB2XWFormat(innerRetinaSU.mosaic{1}.responseLinear)')
% vcNewGraphWin; plot(RGB2XWFormat(iStim.cMosaic.current)')

%% Make GIF
params.vname = [isetbioRootPath '/local/vernierMovieTest.gif'];
% ieGIF(psth(:,:,steadyStateFrame:end),params);