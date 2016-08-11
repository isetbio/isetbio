%% Use the coneMosaic object to simulate responses of foveal RGC mosaic
%
% This tutorial generates RGC responses to a drifting Gabor. This tutorial
% also includes
%
%   * Create a scene, oi and cMosaic
%   * Calculate bipolar
%   * Calculate RGC for on parasol 
%
% Based on t_coneMosaic.
% 
% 7/2016 JRG HJ BW (c) isetbio team

%% Initialize parameters

% clx; ieInit;
clear;
% Initialize parameters of simulated retinal patch
ecc = [0,0]*1e-3;   % Cone mosaic eccentricity in meters from fovea
fov = 2.8;            % Scene Field of view in degrees

params.nSteps = 150;

rdtFlag = 1; % If set, get cone mosaic from RDT; otherwise compute locally
osFlag  = 0; % If set, use osBioPhys; otherwise use osLinear
%%
switch rdtFlag
    case 0
%% Generate iStim structure locally

if osFlag; params.os = 'bioPhys'; end; 
        
params.freq = 3;
params.GaborFlag = 0.3;
iStim = ieStimulusGabor(params);

cMosaic = iStim.cMosaic;

%%
    case 1
%% Get iStim structure for barMovie from RDT
rdt = RdtClient('isetbio');
rdt.crp('/resources/data/istim');

switch osFlag
    case 0
data = rdt.readArtifact('gaborDrifting_cMosaic', 'type', 'mat');
    case 1
data = rdt.readArtifact('gaborDrifting_cMosaic_osBioPhys', 'type', 'mat');
end

iStim = data.iStim; clear data;
cMosaic = iStim.cMosaic;
%%
end
%% Compute the bipolar response

bp = bipolar(cMosaic.os);
bp.set('sRFcenter',1);
bp.set('sRFsurround',1);
bp.compute(cMosaic.os);

%% Set other RGC mosaic parameters

clear params innerRetinaSU
cellType = 'onParasol';
% cellType = 'offParasol';
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
param.FrameRate = 5; params.step = 2; params.show = true;
% vcNewGraphWin; ieMovie(innerRetinaSU.mosaic{1}.responseLinear(:,:,1:120));
vcNewGraphWin; ieMovie(psth(:,:,50:end),params);

% vcNewGraphWin; imagesc(mean(psth,3))

% vcNewGraphWin; plot(RGB2XWFormat(innerRetinaSU.mosaic{1}.responseLinear)')
% vcNewGraphWin; plot(RGB2XWFormat(iStim.cMosaic.current)')
