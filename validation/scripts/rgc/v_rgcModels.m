% v_rgcModels
% 
% Validation script for RGC spiking. 
% 
% Test the RGC firing in response to a simple stimlus, e.g. a moving bar,
% and compare to validation data.
% 
% The RGC response is computed for different types of computational models:
% 
%   1. OS = linear,  bipolar = linear,  RGC = LNP
%   2. OS = linear,  bipolar = linear,  RGC = GLM (coupled)
%   3. OS = linear,  bipolar = subunit, RGC = LNP
%   4. OS = linear,  bipolar = subunit, RGC = GLM (coupled)
%   5. OS = biophys, bipolar = linear,  RGC = LNP
%   6. OS = biophys, bipolar = linear,  RGC = GLM (coupled)
%   7. OS = biophys, bipolar = subunit, RGC = LNP
%   8. OS = biophys, bipolar = subunit, RGC = GLM (coupled)
% 
% Movies of the PSTH are shown for each simulation. Note that there are
% differences between the spiking activity for the different combinations
% of OS, bipolar and RGC. This script compares the responses of each
% simulation from a certain isetbio github commit and checks it with the
% same simulation in the current commit. If there have been no changes that
% affect the value of the PSTH for each simulation, the validationRun
% variable will have entries close to zero (< 0.01).
% 
% 8/2016 JRG (c) isetbio team

%%
ecc = [0,0]*1e-3;   % Cone mosaic eccentricity in meters from fovea
fov = 2.8;          % Scene Field of view in degrees

showMovieFlag = 1;

osFlag      = [0 0 0 0 1 1 1 1]; % 0 = osLinear,        1 = osBioPhys
bipolarFlag = [0 0 1 1 0 0 1 1]; % 0 = bipolar linear,  1 = bipolar rectify   
rgcFlag     = [0 1 0 1 0 1 0 1]; % 0 = rgc uncoupled,   1 = rgc coupled

for flagInd = 1:length(osFlag)
tic
flagInd
%% RDT computation
rdt = RdtClient('isetbio');
rdt.crp('/resources/data/istim');

switch osFlag(flagInd)
    case 0 % osLinear
        data = rdt.readArtifact('barMovie_cMosaic', 'type', 'mat');
    case 1 % osBioPhys
        data = rdt.readArtifact('barMovie_cMosaic_osBioPhys', 'type', 'mat');
end

% We are only using the cMosaic
% sceneRGB = data.iStim.sceneRGB;
cMosaic = data.iStim.cMosaic;
clear data;

%% Compute the bipolar response

switch bipolarFlag(flagInd)
    case 0 % linear bipolar
        bp = bipolar(cMosaic.os);
    case 1 % subunit bipolar
        bp = bipolar(cMosaic.os,'rectifyType',2);
end
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

nTrials = 4; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);

%% Compute the inner retina response

switch rgcFlag(flagInd)
    case 0 % no coupling
        innerRetinaSU = irCompute(innerRetinaSU, bp);
    case 1 % coupling
        innerRetinaSU = irCompute(innerRetinaSU, bp, 'coupling', true);
end

toc

lastTime = innerRetinaSU.mosaic{1}.get('last spike time');

%% Make the PSTH movie
innerRetinaSU.mosaic{1}.set('dt',1);
psthTest{flagInd} = innerRetinaSU.mosaic{1}.get('psth');

% clear params 
% param.FrameRate = 5; params.step = 2; params.show = true;
% 
% % % View movie of RGC linear response
% % vcNewGraphWin; ieMovie(innerRetinaSU.mosaic{1}.responseLinear(:,:,1:120));
% 
% % View movie of PSTH for mosaic
if showMovieFlag
    steadyStateFrame = 50; % Get rid of transient spiking
    vcNewGraphWin; ieMovie(psthTest{1,flagInd}(:,:,steadyStateFrame:end),params);
end

%%

end

load([isetbioRootPath '/local/v_rgcModels_output.mat'],'psth');
for flagInd = 1:length(osFlag)
    validationTemp = psth{1,flagInd}  - psthTest{1,flagInd};
    validationRun(flagInd) = sum(validationTemp(:))./sum(psth{1,flagInd}(:));
end

validationRun

% Unit test framework
% u = UnitTest;
% u.assertIsZero(0,'psth',1e-12);

% Save validation data
% psth = psthTest; save([isetbioRootPath '/local/v_rgcModels_output.mat'],'psthTest');

%% Estimate of bar position
% clear psthColAvg
% vcNewGraphWin
% for flagInd = 8%:8
%     psthSingle = psth{1,flagInd};
%     psthColAvg(:,:) = mean(psthSingle,1);
%     hold on;
%     plot(psthColAvg');
%     
% end