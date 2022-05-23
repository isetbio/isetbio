%% Use the coneMosaic object to simulate responses of foveal RGC mosaic
%
% This tutorial generates RGC responses to a drifting Gabor. This tutorial
% also includes
%
%   * Get precomputed cone mosaic response from RDT
%   * Calculate bipolar
%   * Calculate RGC for on parasol 
%
% Based on t_coneMosaic.
% 
% 7/2016 JRG HJ BW (c) isetbio team

%% Initialize parameters

clx; ieInit;

% Initialize parameters of simulated retinal patch
ecc = [0,0]*1e-3;   % Cone mosaic eccentricity in meters from fovea
fov = 2.8;          % Scene Field of view in degrees
nSteps = 150;       % Number of temporal frames

rdtFlag = 1;        % 0 = compute locally, 1 = pull cMosaic from RDT
osFlag  = 0;        % 0 = osLinear, 1 = osBioPhys

%%
%% Get iStim structure for barMovie from RDT
rdt = RdtClient('isetbio');
rdt.crp('/resources/data/istim');

switch osFlag
    case 0 % osLinear
%         data = rdt.readArtifact('gaborDrifting_cMosaic', 'type', 'mat');        
        data = rdt.readArtifact('gaborMovie', 'type', 'mat');
    case 1 % osBioPhys
%         data = rdt.readArtifact('gaborDrifting_cMosaic_osBioPhys', 'type', 'mat');
        data = rdt.readArtifact('gaborMovie_osBioPhys', 'type', 'mat');
end

iStim = data.iStim; clear data;
cMosaic = iStim.cMosaic;
cMosaic.computeCurrent;
%% Compute the bipolar response

bp = bipolar(cMosaic);
bp.set('sRFcenter',1);
bp.set('sRFsurround',0);
bp.compute(cMosaic);
% bp.plot('movie response')

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

%% Make the PSTH movie
innerRetinaSU.mosaic{1}.set('dt',1);
psth = innerRetinaSU.mosaic{1}.get('psth');

clear params 
param.FrameRate = 3; params.step = 1; params.show = true;

% % View movie of RGC linear response
% vcNewGraphWin; ieMovie(innerRetinaSU.mosaic{1}.responseLinear(:,:,1:120));

% View movie of PSTH for mosaic
steadyStateFrame = 50; % Get rid of transient spiking
vcNewGraphWin; ieMovie(psth(:,:,steadyStateFrame:end),params);

% % View average of PSTH movie
% vcNewGraphWin; imagesc(mean(psth,3))

% % Plots of RGC linear response and OS current
% vcNewGraphWin; plot(RGB2XWFormat(innerRetinaSU.mosaic{1}.responseLinear)')
% vcNewGraphWin; plot(RGB2XWFormat(iStim.cMosaic.current)')

%% Make GIF
params.vname = [isetbioRootPath '/local/gaborMovieTest.gif'];
% ieGIF(psth(:,:,steadyStateFrame:end),params);