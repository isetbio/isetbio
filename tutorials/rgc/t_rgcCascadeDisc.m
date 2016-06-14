
% t_rgcCascadeDisc
% 
% A comparison of RGC firing from "black box" GLM simulations and "cascade
% model" GLM simulations for the linear equivalent disc stimulus. The
% experimental recordings and the black box simulations are found in
% Brackbill, Heitman, ..., Chichilnisky, J. Neurophys 2016. 

% The theory behind the linear equivalent disc stimulus is laid out in
% Turner & Rieke, 2016. If an RGC has purely linear spatial processing, its
% response to the natural scene movie and disc movie will be the same. If
% there is a large difference, that indicates the spatial processing is
% nonlinear, likely due to bipolar subunits.
% 
% See also t_rgcCascade.m, EJLFigureReproduction/t_rgcNaturalScenesFig2.m
% 
% 6/2016 JRG (c) isetbio team

for cellInd = 10% :10

    
    
experimentI   = 1;       % Choose dataset to load parameters and spikes
cellTypeI     = 2;%:2    % Choose On Parasol (1) or Off Parasol (2)
stimulusTestI = 2;%:2     % Choose WN test stimulus (1) or NSEM test stimulus (2)
    
% Switch on the conditions indices
% Experimental dataset
switch experimentI
    case 1; experimentID = '2013-08-19-6';
    otherwise; error('Data not yet available');
%     case 2; experimentID = '2012-08-09-3'; % Data not yet available
%     case 3; experimentID = '2013-10-10-0';
%     case 4; experimentID = '2012-09-27-3';
end
% The other experimental data will be added to the RDT in the future.

% Stimulus: white noise or natural scene movie with eye movements
switch stimulusTestI
    case 1; stimulusTest = 'WN';
    case 2; stimulusTest = 'NSEM';
end

% Cell type: ON or OFF Parasol
switch cellTypeI
    case 1; cellType = 'On Parasol';
    case 2; cellType = 'Off Parasol';
end
    
%%
% load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/saved_iStim_linEqDisc/linEqDisc_10_notemp.mat')

rdt = RdtClient('isetbio');
rdt.crp('resources/data/rgc');
data = rdt.readArtifact('linEqDisc_10_notemp', 'type', 'mat');
linEqDisc = data.linEqDisc;


% Length of WN movie is 1200, take nFrames to limit natural movie to same length
nFrames = 3600; 
testmovieshort = linEqDisc{cellInd}; 
% ieMovie(testmovieshort(:,:,1:100))
frRS = 8;
testmovieRS = zeros(size(testmovieshort,1),size(testmovieshort,2),frRS*nFrames);
for frnum = 1:nFrames
    for frrep = 1:frRS
        testmovieRS(:,:,(frnum-1)*frRS+frrep) = testmovieshort(:,:,frnum);
    end
end


%% Show test movie
showFrames = 200;
% ieMovie(testmovieshort(:,:,1:showFrames));

%% Generate display, scene, oi, sensor
paramsStim.nsteps = 1;%nFrames;%size(testmovieshort,3);
%  Bipolar filter is setfor 0.001 sec, so it needs to be 0.001
paramsStim.timeInterval = 0.001; % sec
paramsStim.expTime = 0.001; % sec

% % For 2013-08-19-6
% r = 12 mm = 36 degs;
% theta = 330 degs;
% side = 'left';
paramsStim.fov = 4;
paramsStim.radius = 36;
paramsStim.theta = 330;
paramsStim.side = 'left';

% iStim = ieStimulusMovie(testmovieshort(:,:,1:nFrames),paramsStim);
iStimDisc = ieStimulusMovie(testmovieshort(:,:,1),paramsStim);
sensorDisc = iStimDisc.sensor;

% sensor.data.volts = 5e-4*double(testmovieshort)./255;
% sensorDisc.data.volts = 5e-4*testmovieRS./max(testmovieRS(:));
sensorDisc.data.volts = 5e-4*-1*(testmovieRS);
% clear testmovieRS testmovieshort
%% Outer segment calculation - linear model
% The iStim structure generates the movie, the scene, the oi and the
% cone absorptions. The next step is to get the outer segment current. The
% linear outer segment model is employed here.

% Initialize
osLDisc = osCreate('linear');
% osL = osSet(osL,'sConeFilter',1);
% osL = osSet(osL,'mConeFilter',1);
% osL = osSet(osL,'lConeFilter',1);
% % % % % Set eccentricity of os patch here

% Set size of retinal patch based on absorptions sensor object
patchSize = sensorGet(sensorDisc,'width','m');
osLDisc = osSet(osLDisc, 'patch size', patchSize);

% Set time step of simulation equal to absorptions sensor object
timeStep = sensorGet(sensorDisc,'time interval','sec');
osLDisc = osSet(osLDisc, 'time step', timeStep);

% Set circular convolution, only steady state
paramsOSL.convolutionType = 1; 

% osLSubDisc = osLDisc;
% 
% % Compute the outer segment response to the absorptions with the linear
% % model.
osLDisc = osCompute(osLDisc,sensorDisc,paramsOSL);
% osL = osSet(osL,'coneCurrentSignal',sensor.data.volts);
% % Plot the photocurrent for a pixel.
% osPlot(osL,sensor);

osLSubDisc = osSet(osLDisc, 'time step', 8*timeStep);
osLSubDisc.osSet('coneCurrentSignal',osLDisc.coneCurrentSignal(:,:,1:8:end));

%% osBioPhys

% % Initialize
osB = osCreate('biophys');

% % % % % Set eccentricity of os patch here

% Set size of retinal patch based on absorptions sensor object
patchSize = sensorGet(sensorDisc,'width','m');
osB = osSet(osB, 'patch size', patchSize);

% Set time step of simulation equal to absorptions sensor object
timeStep = sensorGet(sensorDisc,'time interval','sec');
osB = osSet(osB, 'time step', timeStep);

sensorVolts = sensorGet(sensorDisc,'volts');
paramsOS.bgVolts = 10*mean(sensorVolts(:));
clear sensorVolts

osBSub = osB;
% Compute the outer segment response to the absorptions with the linear
% model.
% osB = osCompute(osB,sensorDisc,paramsOS);
% 
% % % Plot the photocurrent for a pixel.
% % osPlot(osB,sensor);
% 
% % osBSub.osSet('coneCurrentSignal',0);
% 
% osBSub.osSet('coneCurrentSignal',osB.coneCurrentSignal(:,:,1:80:end));
% clear osB
% % osPlot(osBSub,sensorDisc);

%% Find bipolar responses
clear bp os
% % osDisc = osLDis;
osDisc = osLSubDisc;
% osDisc = osBSub;

% bpDisc = bipolar(osDisc);
bpDisc = bipolar(osDisc,cellType,2);
bpDisc = bipolarCompute(bpDisc, osDisc);

% bipolarPlot(bp,'response');

%%
% Set parameters
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = 36; 
params.eyeAngle = 0; ntrials = 0;

% Determined at beginning to allow looping
params.experimentID = experimentID; % Experimental dataset
params.stimulusTest = stimulusTest; % WN or NSEM
params.cellType = cellType;         % ON or OFF Parasol;

% Create object
innerRetinaDisc = irPhys(bpDisc, params);
nTrials = 57; innerRetinaDisc = irSet(innerRetinaDisc,'numberTrials',nTrials);

innerRetinaDisc = irCompute(innerRetinaDisc, bpDisc);
% innerRetinaDisc = irComputeContinous(innerRetinaDisc, bpDisc);
% innerRetinaDisc = irNormalize(innerRetinaDisc, innerRetinaDiscRGB);
% innerRetinaDisc = irComputeSpikes(innerRetinaDisc, bpDisc);
innerRetinaDisc = irSet(innerRetinaDisc,'timing',.0083);

irPlot(innerRetinaDisc,'linear','cell',[7 1]);

innerRetinaDiscPSTH = mosaicGet(innerRetinaDisc.mosaic{1},'responsePsth');


%% Generate outer segment object

% In this case, the RGC GLM calculation converts from the frame buffer
% values in the movie to the spiking responses.  For this calculation, we
% store the movie stimulus in the the outer segment object 'displayRGB'.

os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep', 1/120);

% Attach the movie to the object
os1 = osSet(os1, 'rgbData', -1*double(testmovieshort));

%% Generate RGC object for simulated GLM prediction of response
% Set the parameters for the inner retina RGC mosaic. For the inner retina
% type irPhys, the values for eyeSide, eyeRadius and eyeAngle have no
% effect, because those are dependent on the properties of the retinal
% piece used in the Chichilnisky Lab experiment.

% Set parameters
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = 12; 
params.eyeAngle = 0; ntrials = 0;

% Determined at beginning to allow looping
params.experimentID = experimentID; % Experimental dataset
params.stimulusTest = stimulusTest; % WN or NSEM
params.cellType = cellType;         % ON or OFF Parasol

% Create object
innerRetinaDiscRGB = irPhys(os1, params);
nTrials = 57; innerRetinaDiscRGB = irSet(innerRetinaDiscRGB,'numberTrials',nTrials);

%% Compute the inner retina response

innerRetinaDiscRGB = irCompute(innerRetinaDiscRGB, os1);

% Get the PSTH from the object
innerRetinaPSTH = mosaicGet(innerRetinaDiscRGB.mosaic{1},'responsePsth');

end

%% Generate outer segment object for GLM from RGB scene data

% In this case, the RGC GLM calculation converts from the frame buffer
% values in the movie to the spiking responses.  For this calculation, we
% store the movie stimulus in the the outer segment object 'displayRGB'.

os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep', 1/120);

% Attach the movie to the object
os1 = osSet(os1, 'rgbData', double(testmovieshort));
%% Create a new inner retina object and attach the recorded spikes
% We also want to compare the spikes recorded in the experiment to those
% from the simulation. We create a separate inner retina object with
% isetbio to store these spikes. This makes it easy to use isetbio to plot
% aspects of the response.

[testmovie, xval_mosaic] =  loadDataRGCFigure2(experimentI,stimulusTestI,cellTypeI);

% Create the object.
innerRetinaRecorded = irPhys(os1, params);  

% innerRetinaRecorded = irSet(innerRetinaRecorded,'timing',.008);
innerRetinaRecorded = irSet(innerRetinaRecorded,'numberTrials',nTrials);

% Set the recorded spikes that we got from the RDT into the object.
innerRetinaRecorded = irSet(innerRetinaRecorded,'recordedSpikes',xval_mosaic);
% Get the PSTH using an isetbio routine.
innerRetinaRecordedPSTH = mosaicGet(innerRetinaRecorded.mosaic{1},'responsePsth');

%%

tStart = 1.5;
tEnd = 21;%1*8.5;
cellNum = cellInd;

vcNewGraphWin([],'upperleftbig'); 
subplot(312); hold on;
irPlot(innerRetinaDiscRGB,'raster','cell',[cellNum 1],'hold','on','color','r')
title(sprintf('Black Box, NSEM Disc, off parasol cell [%d 1]',cellNum));
set(gca,'fontsize',14);
axis([tStart tEnd 0 57]);
axis off

subplot(313); hold on;
irPlot(innerRetinaDisc,'raster','cell',[cellNum 1],'hold','on','color','b')
title(sprintf('Cascade Conv, NSEM Disc, off parasol cell [%d  1]',cellNum));
set(gca,'fontsize',14);
% axis([tStart-.04 tEnd-.04 0 57]); % when using theoretical irGLM
axis([tStart tEnd 0 57]);
axis off


subplot(311); hold on;
irPlot(innerRetinaRecorded,'raster','cell',[cellNum 1],'hold','on','color','k')
title(sprintf('Recorded, NSEM Disc, off parasol cell [%d  1]',cellNum));
set(gca,'fontsize',14);

switch stimulusTestI
    case 1
        axis([tStart-0.5 tEnd-0.5 0 57]);
    case 2
        axis([tStart-1 tEnd-1 0 57]);
end
 set(gcf,'position',[ 0.0063   -0.0444    0.8819    0.9378]);
