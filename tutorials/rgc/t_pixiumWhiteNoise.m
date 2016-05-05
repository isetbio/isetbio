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


%%% Grating subunit stimulus

iStim = ieStimulusBinaryWhiteNoise(params);
absorptions = iStim.sensor;
whiteNoise = iStim;
%% Show raw stimulus for osIdentity
figure;
for frame1 = 1:size(whiteNoise.sceneRGB,3)
    imagesc(squeeze(whiteNoise.sceneRGB(:,:,frame1,:)));
    colormap gray; drawnow;
end
% close;

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

timeStep = sensorGet(whiteNoise.sensor,'time interval','sec');
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

innerRetina = irCompute(innerRetina,os);

irPlot(innerRetina, 'linear');
% irPlot(innerRetina, 'psth');

%%

psth1 = mosaicGet(innerRetina.mosaic{1},'responsePsth');

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

psthNorm=psthM'*diag(1./sqrt(sum(psthM'.*psthM')));

psthCov = psthNorm'*psthNorm;
figure; imagesc(psthCov)
xlabel('Cell number'); ylabel('Cell number');
title('Covariance matrix');


rf = zeros(96,96);
for fr = 1:1197
    rf = rf+(sum(psth1{8,4}((fr-1)*100+1:fr*100)))*(iStim.sceneRGB(:,:,fr,1));
%     rf = rf+(1*(iStim.sceneRGB(:,:,fr,1)));
end
figure; imagesc(rf)

% %% Compute RGC spiking
% numberTrials = 1;
% for tr = 1:numberTrials
%     innerRetina = irComputeSpikes(innerRetina);
% end
% %% Invert representation to form image/movie
% clear stimulusReconstruction
% [stimulusReconstruction, paramsRec] = irReconstruct(innerRetina);
% 
% 
% %% Build RGC array for healthy retina
% 
% clear paramsIR innerRetinaHealthy
% paramsIR.name    = 'Macaque inner retina 1'; % This instance
% paramsIR.eyeSide   = 'left';   % Which eye
% paramsIR.eyeRadius = 6;        % Radius in mm
% paramsIR.eyeAngle  = 90;       % Polar angle in degrees
% 
% model   = 'LNP';    % Computational model
% innerRetinaHealthy = irCreate(os,paramsIR);
% innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','onMidget','model',model);
% innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','offMidget','model',model);
% innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','onParasol','model',model);
% innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','offParasol','model',model);
% 
% %%
% innerRetinaHealthy = irComputeContinuous(innerRetinaHealthy,os);
% numberTrials = 1;
% for tr = 1:numberTrials
%     innerRetinaHealthy = irComputeSpikes(innerRetinaHealthy);
% end
% 
% irPlot(innerRetinaHealthy, 'linear');
% % irPlot(innerRetinaHealthy, 'mosaic');
% %% Invert representation to form image/movie
% clear stimulusReconstructionHealthy
% [stimulusReconstructionHealthy, paramsRecHealthy] = irReconstruct(innerRetinaHealthy);
% 
% %%
% name_str = 'gratingH_20Hz_Width12_ONM_ONP_OFFM_OFFP_NoSub.mp4';
% path_str = '/Users/james/Documents/MATLAB/isetbio misc/pixium_videos/meeting_april29/';
% vObj = VideoWriter([path_str name_str],'MPEG-4');
% vObj.FrameRate = 10;
% vObj.Quality = 100;
% open(vObj);
% 
% sizeScene = size(whiteNoise.sceneRGB(:,:,frame1,:));
% 
% % Play the movie with the stimulus
% for loopv = 1%:10
% h1=figure; set(gcf,'position',[160 60 1070 740]);
% hold on;
% for frame1 = 1:params.nSteps%size(whiteNoise.sceneRGB,3)
%     subplot(221);
%     imagesc(squeeze(whiteNoise.sceneRGB(:,:,frame1,:)));
%     colormap gray; 
%     
%     subplot(222);
%     for xPos = 1:numberElectrodesX
%         for yPos = 1:numberElectrodesY
%             hold on;
%             fill(xh+electrodeArray.center(xPos,numberElectrodesY+1-yPos,1),yh+electrodeArray.center(xPos,numberElectrodesY+1-yPos,2),electrodeArray.activation(xPos,yPos,frame1))
%         end
%     end
%     caxis([0 1]);
%     
%     
%     subplot(223);    
% %     imagesc((stimulusReconstructionHealthy(1:paramsRecHealthy.maxx,1:paramsRecHealthy.maxy,frame1)));
%     imagesc((stimulusReconstructionHealthy(1:sizeScene(1),1:sizeScene(2),frame1)));
%      colormap gray
% %      caxis([1*paramsRecHealthy.minR 1*paramsRecHealthy.maxR]);
%     caxis([.5*paramsRecHealthy.minR .5*paramsRecHealthy.maxR]);
% %     caxis([0 .5*paramsRecHealthy.maxR]);
%     title('Healthy');
%     
%     subplot(224);    
% %     imagesc((stimulusReconstruction(1:paramsRec.maxx,1:paramsRec.maxy,frame1)));
%     imagesc((stimulusReconstruction(1:sizeScene(1),1:sizeScene(2),frame1)));
% 
%      colormap gray
% % %     caxis([.5*paramsRec.minR .5*paramsRec.maxR]);
%     caxis([0 1*paramsRec.maxR]);
%     title('Prosthetic');
% %     pause(0.1);
% drawnow
% 
%     F = getframe(h1);
%     writeVideo(vObj,F);
% end
% end
% 
% 
% close(vObj)