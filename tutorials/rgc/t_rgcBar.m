%% Use the coneMosaic object to simulate responses of foveal RGC mosaic
%
% This tutorial generates RGC responses to a moving bar. This tutorial also
% includes
%
%   * Create a scene, oi and cMosaic
%   * Run computation locally or pull from RDT
%   * Select osLinear or osBioPhys for cone current
%   * Calculate bipolar
%   * Calculate RGC for on parasol 
%   * Estimate position of bar from RGC firing
%
% Based on t_coneMosaic.
% 
% 7/2016 JRG HJ BW (c) isetbio team

%% Initialize parameters

clx; ieInit;

% Initialize parameters of simulated retinal patch
ecc = [0,0]*1e-3;     % Cone mosaic eccentricity in meters from fovea
fov = 2.8;            % Scene Field of view in degrees

rdtFlag = 1;          % 0 = compute locally, 1 = pull cMosaic from RDT
osFlag  = 0;          % 0 = osLinear, 1 = osBioPhys

%% Local or RDT computation

switch rdtFlag
    case 0 % Compute locally
        
        %% Generate iStim structure locally
        
        if osFlag; params.os = 'bioPhys'; end;
        
        % This function generates a movie of a bar sweeping from left to right; the
        % movie is converted into an isetbio scene, oi and cone mosaic within the
        % function, and these are returned within the iStim structure.
        iStim = ieStimulusBar(params);
        
        cMosaic = iStim.cMosaic;
        
        %%
    case 1 % Use RDT
        
        %% Get iStim structure for barMovie from RDT
        rdt = RdtClient('isetbio');
        rdt.crp('/resources/data/istim');
        
        switch osFlag
            case 0 % osLinear
                data = rdt.readArtifact('barMovie_cMosaic', 'type', 'mat');
            case 1 % osBioPhys
                data = rdt.readArtifact('barMovie_cMosaic_osBioPhys', 'type', 'mat');
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

%% Make the PSTH movie
innerRetinaSU.mosaic{1}.set('dt',1);
psth = innerRetinaSU.mosaic{1}.get('psth');

clear params 
param.FrameRate = 5; params.step = 2; params.show = true;

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

%% Estimate position of bar

clear colLocation trueLocation

szPsth = size(psth);
sizeRGB = size(iStim.sceneRGB);
sizeMosaic = innerRetinaSU.mosaic{1}.get('mosaicSize');

% Plot estimated position of bar versus true position of bar
frStart = steadyStateFrame; frEnd = 15;
for fr = 1:szPsth(3)-frStart-frEnd
    % Find indices of spiking cells for each frame
    [spikingCells1 spikingCells2] = find(psth(:,:,fr+frStart)==1);
    % Find average column position of spiking RGC for that frame
    colLocation(fr) = mean(spikingCells2);
    
    % Find indicies of true bar pixels for each frame
    [barPx1 barPx2] = find(iStim.sceneRGB(:,:,fr+frStart)>.75);
    % Find average column position
    trueLocation(fr) = mean(barPx2)*sizeMosaic(2)./sizeRGB(2);
end


vcNewGraphWin; scatter(1:szPsth(3)-frStart-frEnd,colLocation-0*min(colLocation),'x');
hold on; plot(1:szPsth(3)-frStart-frEnd,trueLocation);
xlabel('Time (msec)');
ylabel('Position estimate (mm)');
legend('Estimated Position','True Position','location','NW');
