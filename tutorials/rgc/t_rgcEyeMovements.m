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

% clx; ieInit;

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
cMosaic.rows = 144; cMosaic.cols = 176;
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
% params.fov = fov;
params.eyeAngle = 0; ntrials = 0;

% Create RGC object
innerRetinaSU = ir(bp, params);
innerRetinaSU.mosaicCreate('type',cellType,'model','GLM');

%%
nTrials = 1; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);

%% Plot the cone, bipolar and RGC mosaics

% mosaicPlot(innerRetinaSU,bp,sensor,params,cellType,ecc);

%% Compute the inner retina response

innerRetinaSU = irCompute(innerRetinaSU, bp); 
lastTime = innerRetinaSU.mosaic{1}.get('last spike time');

%%
innerRetinaSU.mosaic{1}.set('dt',1);
psth = innerRetinaSU.mosaic{1}.get('psth');

clear params
params.vname = fullfile(isetbioRootPath,'local','vernier.avi'); 
param.FrameRate = 5; params.step = 2; params.show = false;
ieMovie(psth,params);

figure; imagesc(mean(psth,3))