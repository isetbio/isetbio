% t_bipolar
% 
% Implement full retinal pathway with sequential computations of the cone,
% bipolar and RGC responses.
% 
% The stimulus is an impulse. The goal is to examine the impulse repsonse
% of the bipolar object with the differentiator 
% 
% 5/2016 JRG (c) isetbio team

%% Initialize
clear
% ieInit

%% Load image sequence

% WN
load('/Users/james/Documents/MATLAB/akheitman/WN_mapPRJ/Stimuli/BW-8-1-0.48-11111_RNG_16807/testmovie_8pix_Identity_8pix.mat')
testmovieRand = double(testmovie.matrix(:,:,1:1200));
% 
% NSEM
% load('/Users/james/Documents/MATLAB/akheitman/NSEM_mapPRJ/Stimuli/NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat')
% testmovieRand = double(testmovie.matrix(:,:,1:3600));
% testmovieRand = zeros(size(testmovie.matrix,1),size(testmovie.matrix,2),36000);
% for frnum = 1:36000-1
%     testmovieRand(:,:,frnum) = double(testmovie.matrix(:,:,1+floor(frnum/10)));
% end
clear testmovie.matrix

testmovie.matrix = testmovieRand;

% testmovieAvg = zeros(size(testmovieRand));
% for frnum = 1:3600
%     testmovieRandFr = testmovieRand(:,:,frnum);
%     testmovieAvg(:,:,frnum) = mean(testmovieRandFr(:));
% end
% 
% testmovieRand = testmovieAvg;
% testmovie.matrix = testmovieRand;

% 
% paramsMovie.timeInterval = .001;
% paramsMovie.expTime = .001;
% paramsMovie.fov = 1;
% 
% nFrames = 40;
% nSteps = ((1/125)/paramsMovie.timeInterval)*nFrames
% 
% paramsMovie.nSteps = nSteps;
% 



% % iStim = ieStimulusMovie(testmovie.matrix(:,:,61:61+nFrames),paramsMovie);
% load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/iStim_NSEM/WN_320.mat');
% %  save(['/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/iStim_NSEM/samp_p001_frames_61_420.mat'],'iStim');

% % % iStim = data.iStim;
% % % load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/iStim_NSEM/samp_25_frames_2.mat')
% % load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/iStim_NSEM/frames_1.mat')
% absorptions = iStim.sensor; % cone isomerizations


% Set up parameters for stimulus.
timeStep = 1e-3;        % time step
nSamples = .4/timeStep;        % 2000 samples

flashIntens = 50000;    % flash intensity in R*/cone/sec (maintained for 1 bin only)

% Create human sensor.
sensor = sensorCreate('human');
% sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'size', [80 40]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', timeStep);


% Create stimulus
% sensor = sensorSet(sensor, 'photon rate', 1+4000*(testmovie.matrix));

% randMat = rand(size(testmovie.matrix));
% randMat = rand(80,40,3600);
% randMat = zeros(size(testmovie.matrix));
% randMat(22,22,1) = 1;
% % randMat(2,2,:) = rand(size(randMat(2,2,:)));
% % % randMat = repmat(rand(size(randMat(2,2,:))),size(randMat,1),size(randMat,2),1);
% testmovieRand = randMat>0.5; testmovie.matrix = 1*testmovieRand;
sensor = sensorSet(sensor, 'photon rate', 1+4000*(1*testmovieRand));
% sensor = sensorSet(sensor, 'photon rate', 1*(20*testmovieRand));

ct = sensorGet(sensor,'cone type');
sensor = sensorSet(sensor, 'cone type', 3*ones(size(ct)));

vcNewGraphWin; 
% plot(squeeze(iStim.sceneRGB(4,10,:,1)))
% plot(squeeze(sensor.data.volts(4,10,:)))
pRate = sensorGet(sensor,'photon rate');
plot(squeeze(pRate(41,21,:)))
xlabel('Time (msec)','fontsize',18); 
ylabel('Stimulus intensity','fontsize',18)

%% Outer segment calculation - linear model
% The iStim structure generates the movie, the scene, the oi and the
% cone absorptions. The next step is to get the outer segment current. The
% linear outer segment model is employed here.

% Initialize
osL = osCreate('linear');

% % % % % Set eccentricity of os patch here

% Set size of retinal patch based on absorptions sensor object
patchSize = sensorGet(sensor,'width','m');
osL = osSet(osL, 'patch size', patchSize);

% Set time step of simulation equal to absorptions sensor object
timeStep = sensorGet(sensor,'time interval','sec');
osL = osSet(osL, 'time step', timeStep);

% Set circular convolution, only steady state
paramsOSL.convolutionType = 1; 

% Compute the outer segment response to the absorptions with the linear
% model.
osL = osCompute(osL,sensor,paramsOSL);

% % Plot the photocurrent for a pixel.
% osPlot(osL,sensor);

%% osBioPhys

% % os = osBioPhys(sensor); 
% % paramsOS.bgVolts = 10*mean(vectorize(sensorGet(sensor,'volts')));
% % os = osSet(os, 'noiseFlag', noiseFlag);
% % os = osCompute(os, sensor, paramsOS);
% % 
% % % Set size of retinal patch
% % patchSize = sensorGet(sensor,'width','m');
% % os = osSet(os, 'patch size', patchSize);
% 
% % Initialize
osB = osCreate('biophys');

% % % % % Set eccentricity of os patch here

% Set size of retinal patch based on absorptions sensor object
patchSize = sensorGet(sensor,'width','m');
osB = osSet(osB, 'patch size', patchSize);

% Set time step of simulation equal to absorptions sensor object
timeStep = sensorGet(sensor,'time interval','sec');
osB = osSet(osB, 'time step', timeStep);

sensorVolts = sensorGet(sensor,'volts');
paramsOS.bgVolts = 10*mean(sensorVolts(:));

osBSub = osB;
% Compute the outer segment response to the absorptions with the linear
% model.
osB = osCompute(osB,sensor,paramsOS);

% % Plot the photocurrent for a pixel.
% osPlot(osB,sensor);

% osBSub.osSet('coneCurrentSignal',0);

% osBSub.osSet('coneCurrentSignal',osB.coneCurrentSignal(:,:,1:10:end));
%% Find bipolar responses
clear bp os
os = osL;
% os = osBSub;

% bp = bipolar(osL);

bp = bipolar(os);
% bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
% bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);

% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

bp = bipolarCompute(bp, os);

% bipolarPlot(bp,'response');

%% Find RGC responses
% Build and IR object that takes as input the bipolar mosaic.

% Initialize.
clear params 
clear innerRetinaBpSu
params.name      = 'Bipolar with nonlinear subunits'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees
% bp = bipolarSet(bp,'patchSize',2e-4);
innerRetinaBpSu = irCreate(bp, params);

% Create a subunit model for the on midget ganglion cell parameters
innerRetinaBpSu.mosaicCreate('model','Subunit','type','off parasol');
innerRetinaBpSu.mosaic{1}.mosaicSet('numberTrials',40);
% % Uncomment to get rid of spatial nonlinearity
newRectifyFunction = @(x) x;
innerRetinaBpSu.mosaic{1}.mosaicSet('rectifyFunction',newRectifyFunction);

tonicDriveOrig = 2.27;
% tD = innerRetinaBpSu.mosaic{1}.mosaicGet('tonicDrive');
% innerRetinaBpSu.mosaic{1}.mosaicSet('tonicDrive',0.001);

% irPlot(innerRetinaBpSu,'mosaic');

% Compute RGC mosaic responses
innerRetinaBpSu = irComputeContinuous(innerRetinaBpSu, bp);
irPlot(innerRetinaBpSu,'linear')
% irPlot(innerRetinaBpSu,'linear','cell',[4 4])
% axis([0 1200 2 2.5])
% axis([0 1200 -8 10])
% axis([0 1200 22.1 22.8])
% axis([0 1200 -18 20])
% axis([0 1200 -6.9 -6.5])

rLinSU = mosaicGet(innerRetinaBpSu.mosaic{1},'responseLinear');
rLinearSU{1,1,1} = 5*rLinSU{1,1}./max(rLinSU{1,1}) + tonicDriveOrig;
innerRetinaBpSu.mosaic{1}.mosaicSet('responseLinear', rLinearSU);
%%
numberTrials = 40;
for tr = 1:numberTrials
    innerRetinaBpSu = irComputeSpikes(innerRetinaBpSu);
end

% irPlot(innerRetinaBpSu,'psth')

%%
% tStart = 0;
% tEnd = 25;
% 
% vcNewGraphWin([],'upperleftbig'); 
% subplot(211); hold on;
% irPlot(innerRetinaBpSu,'raster','hold','on','color','b')
% title('Cascade Diff, WN, off parasol cell [1 1]');
% set(gca,'fontsize',14);
% axis([tStart tEnd 0 40]);
% subplot(212); 
% irPlot(innerRetinaBpSu,'psth','hold','on')
% ax2 = axis;
% axis([tStart tEnd ax2(3) ax2(4)]);
% set(gca,'fontsize',14);
% set(gcf,'position',[ 0.0063    0.4556    0.9813    0.4378]);
%% Measure the response for the original GLM model
% This response can be compared to the above response 

% Create outer segment
osD = osCreate('displayRGB');

% Set os parameters
coneSpacing = sensorGet(sensor,'width','m');
osD = osSet(osD, 'patchSize', coneSpacing);

coneSampling = sensorGet(sensor,'time interval');
osD = osSet(osD, 'timeStep', coneSampling);

% Set os stimulus
osD = osSet(osD, 'rgbData', testmovieRand);

% Build inner retina
clear params innerRetinaRGB
params.name      = 'Macaque inner retina 1 impulse'; % This instance

params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetinaRGB = irCreate(osD, params);

innerRetinaRGB.mosaicCreate('model','LNP','type','off parasol');

% Compute response
innerRetinaRGB = irComputeContinuous(innerRetinaRGB, osD);

% Plot response
irPlot(innerRetinaRGB,'linear');
%%
% irPlot(innerRetinaRGB,'linear','hold','on');%,'cell',[1 2])
% irPlot(innerRetinaRGB,'linear','cell',[1 2])

rLin = mosaicGet(innerRetinaRGB.mosaic{1},'responseLinear');
rLin{1,1,1} = 5*rLin{1,1}./max(rLin{1,1}) + tonicDriveOrig;
innerRetinaRGB.mosaic{1}.mosaicSet('responseLinear', rLin);

numberTrials = 40;
for tr = 1:numberTrials
    innerRetinaRGB = irComputeSpikes(innerRetinaRGB);
end
% irPlot(innerRetinaRGB,'psth','hold','on');
% % irPlot(innerRetinaRGB,'raster');

% legend('Cascade Model','Black Box Model')
% title('NSEM Off Parasol [1 1]');
% set(gca,'fontsize',14)
% axis([0 30 0 30])

%%
% vcNewGraphWin([],'upperleftbig'); 
% subplot(211); hold on;
% irPlot(innerRetinaRGB,'raster','hold','on','color','r')
% title('Black Box, WN, off parasol cell [1 1]');
% axis([tStart tEnd 0 40]);
% set(gca,'fontsize',14);
% subplot(212); 
% irPlot(innerRetinaRGB,'psth','hold','on')
% ax2 = axis;
% axis([tStart tEnd ax2(3) ax2(4)]);
% set(gca,'fontsize',14);
% set(gcf,'position',[ 0.0063    0.4556    0.9813    0.4378]);


%%
tStart = 0;
tEnd = 35;

vcNewGraphWin([],'upperleftbig'); 
subplot(311); hold on;
irPlot(innerRetinaRGB,'raster','hold','on','color','r')
title('Black Box, WN, off parasol cell [1 1]');
set(gca,'fontsize',14);
axis([tStart tEnd 0 40]);
subplot(312); hold on;
irPlot(innerRetinaBpSu,'raster','hold','on','color','b')
title('Cascade Conv, WN, off parasol cell [1 1]');
set(gca,'fontsize',14);
axis([tStart tEnd 0 40]);
subplot(313); 
irPlot(innerRetinaBpSu,'psth','hold','on','color','b')
irPlot(innerRetinaRGB,'psth','hold','on','color','r')
legend('Cascade Conv','Black Box');
grid on
ax2 = axis;
axis([tStart tEnd ax2(3) ax2(4)]);
set(gca,'fontsize',14);
% set(gcf,'position',[   0.0063    0.2356    0.6861    0.3578]);
set(gcf,'position',[ 0.0063    0.2354    0.7219    0.4549]);

%% Impulse response
% vcNewGraphWin([],'upperleftbig'); 
% % plot(rLin{1,1,1}(1:400)./max(abs(rLin{1,1,1}(1:400))));
% rLinV = rLin{1,1,1}(1:400);
% plot(-(rLinV - mean(rLinV))./max(abs((rLinV - mean(rLinV)))),'linewidth',4)
% hold on;
% rLinSUV = rLinSU{1,1,1}(1:400);
% plot(16+[1:length(rLinSUV)],-(rLinSUV - mean(rLinSUV))./max(abs((rLinSUV - mean(rLinSUV)))),':r','linewidth',4)
% grid on;
% xlabel('Time (msec)'); ylabel('Normalized Linear Response');
% title('Impulse Response for BB and Diff Cascade RGC Models');
% set(gca,'fontsize',14);
% legend('Black Box','Diff Cascade');

%% PSTHs

psth1 = mosaicGet(innerRetinaRGB.mosaic{1},'responsePsth');
psth2 = mosaicGet(innerRetinaBpSu.mosaic{1},'responsePsth');

binsize = 100;
for binval = 1:3500
    psthSU(binval) = mean(psth2.psth{1}((binval-1)*binsize + 1 :binval*binsize));
    psthRGB(binval) = mean(psth1.psth{1}((binval-1)*binsize + 1 :binval*binsize));
    
end
% figure; plot(psthRGB); hold on; plot(psthSU);

figure; scatter(psthRGB,psthSU);

%% Plot outer seg current

% figure; plot(mean(reshape(osL.coneCurrentSignal,40*80,3600))')
% osPlot(osB,sensor);
% osPlot(osB,sensor);
% sensor = sensorSet(sensor, 'photon rate', 1*(20*testmovieRand));
% osPlot(osB,sensor);
% osBSub.osSet('coneCurrentSignal',osB.coneCurrentSignal(:,:,1:10:end));
% figure; plot(mean(reshape(osBSub.coneCurrentSignal,40*80,3600))')
% meanosb =(mean(reshape(osBSub.coneCurrentSignal,40*80,3600))');
% figure; plot((meanosb - mean(meanosb))./max(abs((meanosb - mean(meanosb)))))
% meanlin = (mean(reshape(osL.coneCurrentSignal,40*80,3600))');
% figure; plot((meanlin - mean(meanlin))./max(abs((meanlin - mean(meanlin)))))
% hold on; plot((meanlin - mean(meanlin))./max(abs((meanlin - mean(meanlin)))))
% title('Comparison of OS Current averaged over all cones');

xval = ceil(79*rand); yval = ceil(39*rand);
figure; plot((squeeze( (osBSub.coneCurrentSignal(xval,yval,:) - mean((osBSub.coneCurrentSignal(xval,yval,:))))./max(abs((osBSub.coneCurrentSignal(xval,yval,:) - mean((osBSub.coneCurrentSignal(xval,yval,:))))))) ))
hold on; plot((squeeze( (osL.coneCurrentSignal(xval,yval,:) - mean((osL.coneCurrentSignal(xval,yval,:))))./max(abs((osL.coneCurrentSignal(xval,yval,:) - mean((osL.coneCurrentSignal(xval,yval,:))))))) ))
legend('osBioPhys','osLinear')
xlabel('Time (msec)'); ylabel('Normalized Outer Segment Current');
title('Comparison of OS Current for a single cone');
set(gca,'fontsize',14)
grid on
