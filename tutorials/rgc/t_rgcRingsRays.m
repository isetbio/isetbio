%% Use the coneMosaic object to simulate responses of foveal RGC mosaic
%
% This tutorial generates RGC responses to the isetbio rings and rays
% stimulus with eye movements.
%
%   * Create a scene
%   * Create an optical image
%   * Calculate a cone mosaic of the fixed scene with eye movements 
%   * Calculate bipolar
%   * Calculate RGC for on parasol 
%
% Based on t_coneMosaic and t_rgcAverageFull.
% 
% 7/2016 JRG HJ BW (c) isetbio team

%% Initialize parameters

% clx; ieInit;

% Initialize parameters of simulated retinal patch
ecc = [0,0]*1e-3;   % Cone mosaic eccentricity in meters from fovea
fov = 2;            % Scene Field of view in degrees
emLength = 250;     % Eye movement frames

sceneType = 'rings rays';
cellType = 'on parasol';

osFlag  = 1; % 0 for osLinear, 1 for osBioPhys

%% Build a scene and oi for computing

s = sceneCreate(sceneType);
s = sceneSet(s,'fov',fov);
s = sceneAdjustLuminance(s,10);
vcAddObject(s);

%%
oi = oiCreate;
oi = oiCompute(oi,s);
vcAddObject(oi); % oiWindow;

%% Build a default cone mosaic and compute the OI

if osFlag
    osCM = osBioPhys();
    cMosaic = coneMosaic('center',[0 0]*1e-3,'os',osCM);  % Create the object
else
    cMosaic = coneMosaic('center',[0 0]*1e-3);  % Create the object
end


% cMosaic.rows = 100; cMosaic.cols = 120;
cMosaic.rows = 144; cMosaic.cols = 176;
cMosaic.emGenSequence(emLength);

cMosaic.compute(oi,'currentFlag',true);

% Show the window
% cMosaic.window;

% Examine the outer segment current
% cMosaic.plot('movie absorptions','vname','deleteme.avi','step',5);

%% Compute the bipolar response

bp = bipolar(cMosaic.os);
bp.set('sRFcenter',1);
bp.set('sRFsurround',1);
bp.compute(cMosaic.os);

%% Set other RGC mosaic parameters

clear params innerRetinaSU
params.name = 'macaque phys';
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

%%
innerRetinaSU.mosaic{1}.set('dt',1);
psth = innerRetinaSU.mosaic{1}.get('psth');

clear params
params.vname = fullfile(isetbioRootPath,'local','vernier.avi'); 
param.FrameRate = 5; params.step = 2; params.show = false;
%  vcNewGraphWin; ieMovie(innerRetinaSU.mosaic{1}.responseLinear);
steadyStateFrame = 40;
vcNewGraphWin; ieMovie(psth(:,:,steadyStateFrame:end),params);

vcNewGraphWin; imagesc(mean(psth,3))