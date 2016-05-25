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
% 
% % NSEM
% % load('/Users/james/Documents/MATLAB/akheitman/NSEM_mapPRJ/Stimuli/NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat')
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

randMat = rand(size(testmovie.matrix));
% randMat = zeros(size(testmovie.matrix));
% randMat(2,2,1) = 1;
% randMat(2,2,:) = rand(size(randMat(41,21,:)));
testmovieRand = randMat>0.5; testmovie.matrix = testmovieRand;
sensor = sensorSet(sensor, 'photon rate', 1+4000*(testmovieRand));

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
paramsOSL.convolutionType = 0; 

% Compute the outer segment response to the absorptions with the linear
% model.
osL = osCompute(osL,sensor,paramsOSL);

% % Plot the photocurrent for a pixel.
osPlot(osL,sensor);
%% Find bipolar responses

bp = bipolar(osL);
bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);

% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

bp = bipolarCompute(bp, osL);

bipolarPlot(bp,'response');

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
innerRetinaBpSu.mosaic{1}.mosaicSet('numberTrials',10);
% % Uncomment to get rid of spatial nonlinearity
newRectifyFunction = @(x) x;
innerRetinaBpSu.mosaic{1}.mosaicSet('rectifyFunction',newRectifyFunction);

% innerRetinaBpSu.mosaic{1}.mosaicSet('tonicDrive',0.01);

% irPlot(innerRetinaBpSu,'mosaic');

% Compute RGC mosaic responses
innerRetinaBpSu = irCompute(innerRetinaBpSu, bp);
irPlot(innerRetinaBpSu,'linear')
% irPlot(innerRetinaBpSu,'linear','cell',[4 4])
irPlot(innerRetinaBpSu,'psth')
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
%  irPlot(innerRetinaRGB,'linear');
% irPlot(innerRetinaRGB,'linear','hold','on');%,'cell',[1 2])
% irPlot(innerRetinaRGB,'linear','cell',[1 2])

rLin = mosaicGet(innerRetinaRGB.mosaic{1},'responseLinear');
% figure;
hold on;
plot(24+(1:length(rLin{1,1}(:))),-.15*squeeze(rLin{1,1}(:))+2.7,'r')
% plot(24+(1:length(rLin{1,1}(:))),squeeze(rLin{1,1}(:))-.05,'r')
% plot(squeeze(.1*rLin{1,1}(24:end-24))+2,'r')

