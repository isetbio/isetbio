% t_bipolar
% 
% Implement full retinal pathway with sequential computations of the cone,
% bipolar and RGC responses.
% 
% This is the first implementation of the bioplar object.
% 
% 5/2016 JRG (c) isetbio team

%% Initialize
clear
ieInit

%% Load image sequence

% This is how we compute ... there are other routines for creating the
% iStim data.  We should spend some time moving them up to RDT and then
% documenting.
% rdt = RdtClient('isetbio');
% rdt.crp('/resources/data/istim');
% data = rdt.readArtifact('iStim_NSEM_800fr','type','mat');
% iStim = data.iStim;


load('/Users/james/Documents/MATLAB/akheitman/NSEM_mapPRJ/Stimuli/NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat')

paramsMovie.nSteps = nSteps;
paramsMovie.timeInterval = .001;
paramsMovie.expTime = .001;
paramsMovie.fov = 1;

nFrames = 360;
nSteps = ((1/125)/paramsMovie.timeInterval)*nFrames

iStim = ieStimulusMovie(testmovie.matrix(:,:,61:61+nFrames),paramsMovie);

%  save(['/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/iStim_NSEM/samp_p001_frames_61_420.mat'],'iStim');

% % iStim = data.iStim;
% % load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/iStim_NSEM/samp_25_frames_2.mat')
% load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/iStim_NSEM/frames_1.mat')
absorptions = iStim.sensor; % cone isomerizations

vcNewGraphWin; 
plot(squeeze(iStim.sceneRGB(4,10,:,1)))
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
patchSize = sensorGet(absorptions,'width','m');
osL = osSet(osL, 'patch size', patchSize);

% Set time step of simulation equal to absorptions sensor object
timeStep = sensorGet(absorptions,'time interval','sec');
osL = osSet(osL, 'time step', timeStep);

% Set circular convolution, only steady state
paramsOSL.convolutionType = 0; 

% Compute the outer segment response to the absorptions with the linear
% model.
osL = osCompute(osL,absorptions,paramsOSL);

% % Plot the photocurrent for a pixel.
osPlot(osL,absorptions);

%% Outer segment calculation - biophysical model
% % As an alternative to the linear model, we can employ Fred's biophysical
% % model that captures properties of temporal adaptation.
% osBp = osCreate('BioPhys');
% 
% % Set size of retinal patch
% patchSize = sensorGet(absorptions,'width','m');
% osBp = osSet(osBp, 'patch size', patchSize);
% 
% % Set time step of simulation equal to absorptions
% timeStep = sensorGet(absorptions,'time interval','sec');
% osBp = osSet(osBp, 'time step', timeStep);
% 
% % Compute the outer segment response to the absorptions with the
% % biophysical model.
% 
% absorptionsVolts = sensorGet(absorptions,'volts');
% paramsOS.bgVolts = 10*mean(absorptionsVolts(:));
% osBp = osCompute(osBp,absorptions,paramsOS);
% 
% % Plot the photocurrent for a pixel.
% 
% % Some changes needed to the biophys plot function ...
% % Let's change the R* to R*/sec.
% % Use isomerizations and current where appropriate
% osPlot(osBp,absorptions);
%% Find bipolar responses
% The bipolar object takes as input the outer segment current. Bipolar
% processing consists of a spatial convolution and a temporal
% differentiator that matches the IR of the measured RGC output.

% Choose which outer segment to use for input.
os = osL;
% os = osBp;

% Create the bipolar object
bp = bipolar(os);

% Compute the bipolar response
bp = bipolarCompute(bp, os);

% Plot the response of every bipolar cell.
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

innerRetinaBpSu = irCreate(bp, params);

% Create a subunit model for the on midget ganglion cell parameters
innerRetinaBpSu.mosaicCreate('model','Subunit','type','off parasol');
innerRetinaBpSu.mosaic{1}.mosaicSet('numberTrials',10);
% % Uncomment to get rid of spatial nonlinearity
newRectifyFunction = @(x) x;
innerRetinaBpSu.mosaic{1}.mosaicSet('rectifyFunction',newRectifyFunction);

innerRetinaBpSu.mosaic{1}.mosaicSet('tonicDrive',0.01);

% irPlot(innerRetinaBpSu,'mosaic');

% Compute RGC mosaic responses
innerRetinaBpSu = irCompute(innerRetinaBpSu, bp);
% for tr = 1:58
%     tic
% innerRetinaBpSu = irComputeSpikes(innerRetinaBpSu, bp);
% toc
% end

% irPlot(innerRetinaBpSu, 'linear');
irPlot(innerRetinaBpSu, 'psth');
% irPlot(innerRetinaBpSu, 'raster','cell',[2 2]);

% irPlot(innerRetinaBpSu, 'raster','cell',[5 1]);
% irPlot(innerRetinaBpSu, 'psth','cell',[5 1]);

toc

%%
% figure; 
% cellNumber = [5, 3];
% minlen = min([length(rgcLpsth.psth{cellNumber(1),cellNumber(2)}) length(rgcBPpsth2.psth{cellNumber(1),cellNumber(2)})]); 
% scatter(rgcLpsth.psth{cellNumber(1),cellNumber(2)}(1:100:minlen),rgcBPpsth2.psth{cellNumber(1),cellNumber(2)}(1:100:minlen))
% 
% % rgcLpsth = mosaicGet(innerRetinaBpSu.mosaic{1},'psthresponse');
% % rgcBPpsth2 = mosaicGet(innerRetinaBpSu.mosaic{1},'psthresponse');
% for cellNumber1 = 1:7
%     for cellNumber2 = 1:7
%         cellNumber = [cellNumber1 cellNumber2];
%         minlen = min([length(rgcLpsth.psth{cellNumber(1),cellNumber(2)}) length(rgcBPpsth2.psth{cellNumber(1),cellNumber(2)})]); 
% 
%         ccMat = corrcoef([rgcLpsth.psth{cellNumber(1),cellNumber(2)}(1:minlen)' rgcBPpsth2.psth{cellNumber(1),cellNumber(2)}(1:minlen)'])
%         cc(cellNumber(1),cellNumber(2)) = ccMat(2,1);
%     end
% end
% 
% figure; stem(sort(cc(:)))

