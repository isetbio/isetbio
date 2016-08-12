%% Use the coneMosaic object to simulate responses of foveal RGC mosaic
%
% This tutorial generates RGC responses to a drifting Gabor. This tutorial
% also includes
%
%   * Create a scene, oi and cMosaic
%   * Run computation locally or pull from RDT
%   * Select osLinear or osBioPhys for cone current
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
switch rdtFlag
    case 0 % Compute locally
        
        %% Generate iStim structure locally
        
        if osFlag; params.os = 'bioPhys'; end;
        
        params.freq = 3;
        params.nSteps = nSteps;
        params.GaborFlag = 0.3;
        
        % This function generates a movie of a Gabor drifing in phase; the
        % movie is converted into an isetbio scene, oi and cone mosaic
        % within the function, and these are returned within the iStim
        % structure.

        iStim = ieStimulusGabor(params);
        
        cMosaic = iStim.cMosaic;
        
        %%
    case 1 % Use RDT
        
        %% Get iStim structure for barMovie from RDT
        rdt = RdtClient('isetbio');
        rdt.crp('/resources/data/istim');
        
        switch osFlag
            case 0 % osLinear
                data = rdt.readArtifact('gaborDrifting_cMosaic', 'type', 'mat');
            case 1 % osBioPhys
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
