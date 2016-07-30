%% Use the coneMosaic object to simulate responses of foveal RGC mosaic
%
% This tutorial generates RGC responses to static image. This tutorial also
% includes
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

clx; ieInit;

% Initialize parameters of simulated retinal patch
ecc = [0,0]*1e-3;   % Cone mosaic eccentricity in meters from fovea
fov = 2;            % Scene Field of view in degrees

sceneType = 'rings rays';
cellType = 'on parasol';

%% Build a scene and oi for computing

s = sceneCreate(sceneType);
s = sceneSet(s,'fov',fov);
s = sceneAdjustLuminance(s,10);
vcAddObject(s);

oi = oiCreate;
oi = oiCompute(oi,s);
vcAddObject(oi); % oiWindow;

%% Build a default cone mosaic and compute the OI

cMosaic = coneMosaic('center',[0 0]*1e-3);  % Create the object
% cMosaic.rows = 100; cMosaic.cols = 120;
% cMosaic.rows = 144; cMosaic.cols = 176;
cMosaic.emGenSequence(500);

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

%% Set RGC mosaic parameters
% 
% experimentID = 'RPE_201602171';
% stimulusTest = 'bar';

% experimentI = 1;
% stimulusTestI = 1;
% cellTypeI = 1;
% 
% % Switch on the conditions indices
% % Experimental dataset
% switch experimentI
%     case 1; experimentID = 'RPE_201602171';
%     otherwise; error('Data not yet available');
% end
% % The other experimental data will be added to the RDT in the future.
% 
% % Stimulus: white noise or natural scene movie with eye movements
% switch stimulusTestI
%     case 1; stimulusTest = 'bar';
% end
% 
% % Cell type: ON or OFF Parasol
% switch cellTypeI
%     case 1; 
%         cellType = 'On Parasol RPE';      
%     case 2; 
%         cellType = 'Off Parasol RPE';        
%     case 3; 
%         cellType = 'On Midget RPE';        
%     case 4; 
%         cellType = 'Off Midget RPE';
%     case 5; 
%         cellType = 'SBC RPE';
%     case 6;
%         cellType = 'On Parasol Apricot';        
%     case 7;
%         cellType = 'Off Parasol Apricot';        
%     case 8;
%         cellType = 'On Midget Apricot';        
%     case 9;
%         cellType = 'Off Midget Apricot';        
%     case 10;
%         cellType = 'SBC Apricot';
%     case 11 
%         cellType = 'on parasol';
%     case 12 
%         cellType = 'off parasol';
%     otherwise;
%         cellType = 'On Parasol RPE';
% end

%% Set other RGC mosaic parameters

clear params innerRetinaSU
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = sqrt(sum(ecc.^2)); 
params.fov = fov;
params.eyeAngle = 0; ntrials = 0;

% Determined at beginning by user
% params.experimentID = experimentID; % Experimental dataset
% params.stimulusTest = stimulusTest; % WN or NSEM
params.cellType = cellType;         % ON or OFF Parasol;

% Set flag for average mosaic
params.averageMosaic = 1;

% Tell the RGC mosaic about how many bipolars per cone
params.inputScale = size(bp.sRFcenter,1);
params.inputSize = size(bp.responseCenter);

% Create RGC object
innerRetinaSU = ir(bp, params);
innerRetinaSU.mosaicCreate('type',cellType,'model','LNP');

%%
nTrials = 10; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);

%% Plot the cone, bipolar and RGC mosaics

% mosaicPlot(innerRetinaSU,bp,sensor,params,cellType,ecc);

%% Compute the inner retina response

innerRetinaSU = irCompute(innerRetinaSU, bp); 

% Get the PSTH from the object
innerRetinaSUPSTH = mosaicGet(innerRetinaSU.mosaic{1},'responsePsth');

% Plot all of the PSTHs together
figure; plot(vertcat(innerRetinaSUPSTH{:})')
title(sprintf('%s Simulated Mosaic at %1.1f\\circ Ecc\nMoving Bar Response',cellType(1:end-4),ecc));
xlabel('Time (msec)'); ylabel('PSTH (spikes/sec)');
set(gca,'fontsize',14);
lastSpikeTime=innerRetinaSU.mosaic{1}.mosaicGet('lastspiketime')
axis([0 lastSpikeTime 0 max(max(vertcat(innerRetinaSUPSTH{:})))]);
grid on;

%% Make a movie of the PSTH response

psthMovie = mosaicMovie(innerRetinaSUPSTH,innerRetinaSU, params);
% figure; ieMovie(psthMovie);