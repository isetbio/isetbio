% t_pixiumWhiteNoise
% 
% Run binary white noise through the RGC array set up for the Pixium
% prosethesis.


%% Initialize
clear;
ieInit;

%% Parameters to alter

% Retinal patch eccentricity
patchEccentricity = 12; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 2.7;

% Stimulus length
nSteps = 1200;

% Activation curve

% Spatial activation of electrode

% Electrode PWM 


%% Load image
clear params
% One frame of a moving bar stimulus
% Set parameters for size
params.nSteps = nSteps;
params.row = 96;
params.col = 96;
params.fov = fov;
% % params.vfov = 0.7;

for blockNum = 201:300

clear psthNorm spikesout spikesoutM spikesoutsm whiteNoiseSmall whiteNoise iStim absorptions

blockNum
%%% Grating subunit stimulus

iStim = ieStimulusBinaryWhiteNoise(params);
absorptions = iStim.sensor;
whiteNoise = iStim;
%% Show raw stimulus for osIdentity
% % figure;
% % for frame1 = 1:size(whiteNoise.sceneRGB,3)
% %     imagesc(squeeze(whiteNoise.sceneRGB(:,:,frame1,:)));
% %     colormap gray; drawnow;
% % end
% % % close;

%% Outer segment calculation
% There is no simulated outer segment, this identity outer segment acts as
% a pass-through holder of the stimulus intensity information.

% Input = RGB
os = osCreate('displayrgb');

sceneSize = sceneGet(whiteNoise.scene,'size');
retinalPatchWidth = sensorGet(whiteNoise.sensor,'width','m');
% retinalPatchHeight = sensorGet(whiteNoise.absorptions,'height','m');
retinalPatchHeight = (sceneSize(1)/sceneSize(2))*retinalPatchWidth;

% % % coneSpacing = scene.wAngular*300
% coneSpacing = sensorGet(sensor,'dimension','um');
os = osSet(os, 'patchSize', retinalPatchWidth);

% timeStep = sensorGet(whiteNoise.sensor,'time interval','sec');
timeStep = 1/120;
os = osSet(os, 'timeStep', timeStep);

os = osSet(os, 'rgbData', whiteNoise.sceneRGB);
% os = osCompute(absorptions);

% % Plot the photocurrent for a pixel
% osPlot(os,absorptions);

retinalPatchSize = osGet(os,'size');

%% Build RGC array

clear paramsIR innerRetina
paramsIR.name    = 'Macaque inner retina 1'; % This instance
paramsIR.eyeSide   = 'left';   % Which eye
paramsIR.eyeRadius = 6;        % Radius in mm
paramsIR.eyeAngle  = 90;       % Polar angle in degrees

model   = 'LNP';    % Computational model
innerRetina = irCreate(os,paramsIR);
% innerRetina = rgcMosaicCreate(innerRetina,'type','onMidget','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','offMidget','model',model);
innerRetina = rgcMosaicCreate(innerRetina,'type','onParasol','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','offParasol','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','sbc','model',model);

irPlot(innerRetina,'mosaic');
% % figure;

innerRetina = irSet(innerRetina,'numberTrials',1);

innerRetina = irCompute(innerRetina,os);

% irPlot(innerRetina, 'linear');
% irPlot(innerRetina, 'psth');

%% Look at covariance matrix
% load('/Users/james/Documents/MATLAB/isetbio misc/optimal linear decoder/ws_pixiumWhiteNoise_May4.mat')
psthstruct = mosaicGet(innerRetina.mosaic{1},'responsePsth');

psth1 = psthstruct.psth;
spikesout = psthstruct.spikes;

szCells = size(psth1);
cellCtr = 0;
szEnd = length(psth1{1,1});
for i2 = 1:szCells(1)
    for j2 = 1:szCells(2)
        cellCtr = cellCtr+1;
        minend = min([szEnd length(psth1{i2,j2})]);
        psthM(cellCtr,1:minend) = psth1{i2,j2}(1:minend);
    end
end

for i2 = 1:szCells(1)*szCells(2)
psthM(i2,:) = psthM(i2,:) - mean(psthM(i2,:));
end

psthNorm=psthM'*diag(1./sqrt(sum(psthM'.*psthM')));
% psthNorm = psthNorm - mean(psthNorm(:));
psthCov = psthNorm'*psthNorm;
% figure; imagesc(psthCov)
% xlabel('Cell number'); ylabel('Cell number');
% title('Covariance matrix');


%%%%%
% load('/Users/james/Documents/MATLAB/isetbio misc/optimal linear decoder/WNstim_response_OnParasol_spikes.mat')
for i2 = 1:szCells(1)*szCells(2)
    spikesoutM(i2,:) = spikesout(i2,:) - mean(spikesout(i2,:));
end

psthNorm=spikesoutM'*diag(1./sqrt(sum(spikesoutM'.*spikesoutM')));
% psthNorm = psthNorm - mean(psthNorm(:));
psthCov = psthNorm'*psthNorm;
% figure; imagesc(psthCov)
% xlabel('Cell number'); ylabel('Cell number');
% title('Covariance matrix');

%%%%%

% PSTHs are here:
psthM;

% Stimulus here:
whiteNoise.sceneRGB;

whiteNoiseSmall = uint8(squeeze(whiteNoise.sceneRGB(:,:,:,1)));
% responseSpikes = mosaicGet(innerRetina.mosaic{1},'responseSpikes');   

spikesoutsm = uint8(spikesout);
filename1 = ['/Users/james/Documents/MATLAB/isetbio misc/optimal linear decoder/WNstim_response_OnParasol_block' num2str(blockNum) '.mat'];
save(filename1, 'whiteNoiseSmall','spikesoutsm');

end

% save('WNstim_response_OnParasol.mat','whiteNoiseSmall','responseSpikes');

% 
% rf = zeros(96,96);
% for fr = 1:1197
%     rf = rf+(sum(psth1{8,4}((fr-1)*100+1:fr*100)))*(iStim.sceneRGB(:,:,fr,1));
% %     rf = rf+(1*(iStim.sceneRGB(:,:,fr,1)));
% end
% figure; imagesc(rf)