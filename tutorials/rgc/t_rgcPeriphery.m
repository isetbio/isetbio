% t_rgcPeriphery
%
% A comparison of RGC firing from experimental recordings, "black box" GLM
% simulations and "cascade model" GLM simulations. The experimental
% recordings are from the dataset RPE 2016-02-17-1.
% 
% The tutorial begins with the user selecting the RGC type. The fitted
% parameters and cross validation data are then loaded from the remote data
% toolbox. The isetbio simulation starts with the movie stimulus, and
% generates the dynamic scene, optical image, cone isomerizations, bipolar
% resopnse and RGC spikes for the cascade model. Isetbio also handles the
% black box model.
% 
% This demonstrates that the linear cascade model matches the output of the
% black box model. The cascade model allows us to turn on nonlinearities in
% the outer segment current computation as well as spatial bipolar
% subunits. We can quantify whether these nonlinearities improve the
% prediction of the firing in response to natural scene movies and compare
% them with the experimental recordings.
% 
% There is a dependency on the repository isetbio/EJLFigureReproduction.
% This must be added to the Matlab path.
% 
% See also EJLFigureReproduction/t_rgcNaturalScenesFig2.m
% 
% 6/2016 JRG (c) isetbio team

%% Initialize 
clear

% Switch on input type
% White noise (WN) or natural scenes with eye movements (NSEM)

experimentI   = 1;       % Choose dataset to load parameters and spikes
cellTypeI     = 1;       % Choose On Parasol (1) or Off Parasol (2)
stimulusTestI = 1;       % Choose WN test stimulus (1) 
    
% Switch on the conditions indices
% Experimental dataset
switch experimentI
    case 1; experimentID = 'RPE_201602171';
    case 2; experimentID = '2013-10-19-6';
    otherwise; error('Data not yet available');
end
% Other experimental data will be added to the RDT in the future.

% Stimulus: white noise or natural scene movie with eye movements
switch stimulusTestI
    case 1; stimulusTest = 'WN';
    case 2; stimulusTest = 'NSEM';
    otherwise
        error('NSEM not yet implemented');
end

% Initialize the RDT object
rdt = RdtClient('isetbio');

% Cell type: ON or OFF Parasol
switch cellTypeI
    case 1; 
        cellType = 'On Parasol RPE';      
        % Load the cross validation data with the RDT
        rdt.crp('resources/data/rgc/rpe_dataset');
        data = rdt.readArtifact('xval_mosaic_WN_ONParasol_201602171', 'type', 'mat');
        xval_mosaic = data.xval_mosaic;
    case 2; 
        cellType = 'Off Parasol RPE';    
        % Load the cross validation data with the RDT
        rdt.crp('resources/data/rgc/rpe_dataset');
        data = rdt.readArtifact('xval_mosaic_WN_OFFParasol_201602171', 'type', 'mat');
        xval_mosaic = data.xval_mosaic;
    case 3; 
        cellType = 'On Midget RPE';
        % Load the cross validation data with the RDT
        rdt.crp('resources/data/rgc/rpe_dataset');
        data = rdt.readArtifact('xval_mosaic_WN_ONMidget_201602171', 'type', 'mat');
        xval_mosaic = data.xval_mosaic;
    case 4; 
        cellType = 'Off Midget RPE';     
        % Load the cross validation data with the RDT
        rdt.crp('resources/data/rgc/rpe_dataset');
        data = rdt.readArtifact('xval_mosaic_WN_OFFMidget_201602171', 'type', 'mat');
        xval_mosaic = data.xval_mosaic;
    case 5; 
        cellType = 'SBC RPE';     
        % Load the cross validation data with the RDT
        rdt.crp('resources/data/rgc/rpe_dataset');
        data = rdt.readArtifact('xval_mosaic_WN_ONSBC_201602171', 'type', 'mat');
        xval_mosaic = data.xval_mosaic;
end
frameRate = 121;
% Load OS from RDT (1) or run from scratch (0)
loadOS = 0;
%% Load stimulus movie and fit/spiking data using RemoteDataToolbox

% Loads the appropriate movie and spiking data for the experimental
% conditions.

rdt = RdtClient('isetbio');
rdt.crp('resources/data/rgc/rpe_dataset');
data = rdt.readArtifact('WN_testmovie_lnfit', 'type', 'mat');
WN_testmovie_lnfit = data.WN_testmovie_lnfit;

% Length of WN movie is 1100
nFrames = 1100; 
testmovieshort = WN_testmovie_lnfit(:,:,1:nFrames); 


if ~loadOS 
    
% Resample the movie to allow the biophysical outer segment current computation
% The biophysical model becomes unstable if the temporal step is too large    
frRS = 8;
testmovieRS = zeros(size(testmovieshort,1),size(testmovieshort,2),frRS*size(testmovieshort,3));
for frnum = 1:nFrames
    for frrep = 1:frRS
        testmovieRS(:,:,(frnum-1)*frRS+frrep) = testmovieshort(:,:,frnum);
    end
end

%% Show test movie
showFrames = 50;
ieMovie(testmovieshort(:,:,1:showFrames));

%% Generate display, scene, oi, sensor
paramsStim.nsteps = 1;

%  Bipolar filter is setfor 0.001 sec, so it needs to be 0.001
paramsStim.timeInterval = 1/frameRate;%0.001; % sec
paramsStim.expTime = 1/frameRate;%0.001; % sec

% % For 2013-08-19-6
% r = 12 mm = 36 degs;
% theta = 330 degs;
% side = 'left';
paramsStim.fov = 8;
paramsStim.radius = 36*0.3;
paramsStim.theta = 330;
paramsStim.side = 'left';

% Get the cone absorptions for the movie with the ieStimulusMovie function
iStim = ieStimulusMovie(permute(testmovieshort(:,:,1:nFrames), [2 1 3]),paramsStim);
sensor = iStim.sensor; sensor.data.volts = 5e-4*double(testmovieshort)./255;
sensor = sensorSet(sensor, 'cone type',3*ones([size(testmovieshort,1) size(testmovieshort,2)]));

clear testmovieRS % testmovieshort

%% Outer segment calculation - linear model
% The iStim structure generates the movie, the scene, the oi and the
% cone absorptions. The next step is to get the outer segment current. The
% linear outer segment model is employed here.

% Initialize
osL = osCreate('linear');
% osL = osSet(osL,'sConeFilter',1);

% Set size of retinal patch based on absorptions sensor object
patchSize = sensorGet(sensor,'width','m');
osL = osSet(osL, 'patch size', patchSize);

% Set time step of simulation equal to absorptions sensor object
timeStep = sensorGet(sensor,'time interval');
osL = osSet(osL, 'time step', timeStep);


% Set circular convolution, only steady state
paramsOSL.convolutionType = 1; 
paramsOSL.ecc = 15; % mm
paramsOSL.singleType = 1;

osLSub = osL;

% Compute the outer segment response to the absorptions with the linear
% model.
osL = osCompute(osL,sensor,paramsOSL);
osPlot(osL,sensor);
%% osBioPhys

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
clear sensorVolts

osBSub = osB;
% Compute the outer segment response to the absorptions with the linear
% model.
% osB = osCompute(osB,sensor,paramsOS);
% 
% % Plot the photocurrent for a pixel.
% osPlot(osB,sensor);

else
    rdt = RdtClient('isetbio');
    rdt.crp('resources/data/rgc');
    data = rdt.readArtifact('full_osLSub_sensor_NSEM', 'type', 'mat');
    osLSub = data.osLSub;
    
end


%% Find bipolar responses
clear bp os
% os = osL;
os = osLSub;
% os = osBSub;

% Select bipolar filter type
bpParams.cellType = 'onDiffuse';
% sets filter as theoretical, mean physiology, or individual phys:
bpParams.filterType = 1; 
% sets linear, on half-wave rectification, or on and off half-wave rect
bpParams.rectifyType = 1;
% bpParams.rectifyType = 3;

bp = bipolar(os, bpParams);

bp.bipolarSet('sRFcenter',[1]);
bp.bipolarSet('sRFsurround',[1]);

bp = bipolarCompute(bp, os);

bipolarPlot(bp,'response');

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

% Choose the cell to plot
params.cellIndices = 7;

% Create object
innerRetinaSU = irPhys(bp, params);
nTrials = 30; innerRetinaSU.mosaic{1} = innerRetinaSU.mosaic{1}.set('numberTrials',nTrials);
%% Compute the inner retina response

% Linear convolution
innerRetinaSU = irCompute(innerRetinaSU, bp); 

innerRetinaSU = irSet(innerRetinaSU,'timing',1/frameRate);

% Get the PSTH from the object
innerRetinaSUPSTH = mosaicGet(innerRetinaSU.mosaic{1},'responsePsth');

%% Generate outer segment object for GLM from RGB scene data

% In this case, the RGC GLM calculation converts from the frame buffer
% values in the movie to the spiking responses.  For this calculation, we
% store the movie stimulus in the the outer segment object 'displayRGB'.

os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep', 1/frameRate);

% Attach the movie to the object
os1 = osSet(os1, 'rgbData', double(testmovieshort));

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

params.cellIndices = 7;

% Create object
innerRetina = irPhys(os1, params);
nTrials = 30; innerRetina.mosaic{1} = innerRetina.mosaic{1}.set('numberTrials',nTrials);

%% Compute the inner retina response

innerRetina = irCompute(innerRetina, os1);

% Get the PSTH from the object
innerRetinaPSTH = mosaicGet(innerRetina.mosaic{1},'responsePsth');

% innerRetinaSU = irNormalize(innerRetinaSU,innerRetina);

%% Create a new inner retina object and attach the recorded spikes
% We also want to compare the spikes recorded in the experiment to those
% from the simulation. We create a separate inner retina object with
% isetbio to store these spikes. This makes it easy to use isetbio to plot
% aspects of the response.

% Create the object.
innerRetinaRecorded = irPhys(os1, params);  

% innerRetinaRecorded = irSet(innerRetinaRecorded,'timing',.008);
nTrials = 30; innerRetinaRecorded.mosaic{1} = innerRetinaRecorded.mosaic{1}.set('numberTrials',nTrials);

% Set the recorded spikes that we got from the RDT into the object.
innerRetinaRecorded = irSet(innerRetinaRecorded,'recordedSpikes',xval_mosaic);
% Get the PSTH using an isetbio routine.
innerRetinaRecordedPSTH = mosaicGet(innerRetinaRecorded.mosaic{1},'responsePsth');



%% Plot the rasters

% Set the time and cell number
tStart = 1.5;% 9%1.5;
tEnd = 9;%21;%18%21;%1*8.5;
cellNum = 1;

% Plot the original GLM prediction
% vcNewGraphWin([],'upperleftbig'); 
h1=figure; hold on;
subplot(312); hold on;
% subplot(211); hold on;
irPlot(innerRetina,'raster','cell',[1 1],'hf',h1,'dt',0.1,'color','r')
% irPlot(innerRetina,'raster','cell',[cellNum 1],'hold','on','color','r')
title(sprintf('Black Box, WN, off parasol cell [%d 1]',cellNum));
% set(gcf,'position',[0.0069    0.6667    0.9750    0.2267]);
set(gca,'fontsize',14);
axis([tStart tEnd 0 30]);
% axis off

% Plot the biophys/subunit prediction
subplot(313); hold on;
% subplot(212); hold on;
irPlot(innerRetinaSU,'raster','cell',[1 1],'hf',h1,'dt',0.1,'color','b');
% irPlot(innerRetinaSU,'raster','cell',[cellNum 1],'hold','on','color','b')
title(sprintf('Cascade Conv, WN, off parasol cell [%d  1]',cellNum));
% set(gcf,'position',[0.0069    0.6667    0.9750    0.2267]);
set(gca,'fontsize',14);

axis([tStart tEnd 0 30]);

% % Plot the recorded spikes
subplot(311); hold on;
irPlot(innerRetinaRecorded,'raster','cell',[params.cellIndices 1],'hf',h1,'dt',0.1,'color','k')
% irPlot(innerRetinaRecorded,'raster','cell',[params.cellIndices 1],'hold','on','color','k')
title(sprintf('Recorded, WN, off parasol cell [%d  1]',params.cellIndices));
% set(gcf,'position',[0.0069    0.6667    0.9750    0.2267]);
set(gca,'fontsize',14);

switch stimulusTestI
    case 1
        axis([tStart tEnd 0 30]);
    case 2
        axis([tStart tEnd 0 30]);
end
%  set(gcf,'position',[ 0.0063   -0.0444    0.8819    0.9378]);

 %% Plot the PSTHs
 
figure;

% subplot(414)
minlen = min([length(innerRetinaPSTH{cellNum}) length(innerRetinaRecordedPSTH{cellNum}) length(innerRetinaSUPSTH{cellNum}) ]);
hold off
switch stimulusTestI
    case 1
        plot((00+[1:minlen-1200])./1208, innerRetinaPSTH{cellNum}(0+(1:minlen-1200)),'r','linewidth',3);
        
        hold on;
        plot([1:minlen-1200]./1208,innerRetinaRecordedPSTH{cellNum}((0+(1:minlen-1200))),':k','linewidth',2);
        
        plot((00+[1:minlen-1200])./1208, innerRetinaSUPSTH{cellNum}(0-0+(1:minlen-1200)),':b','linewidth',3);
        
        ax3 = axis;
        axis([0 8.5 ax3(3) ax3(4)])

    case 2
        plot((00+[1:minlen-1200])./1208, innerRetinaPSTH{cellNum}(0+(1:minlen-1200)),'r','linewidth',3);
        hold on;
        plot([1:minlen-1200]./1208,innerRetinaRecordedPSTH{cellNum}((0+(1:minlen-1200))),':k','linewidth',2);
        hold on;
        plot((00+[1:minlen-1200])./1208, innerRetinaSUPSTH{cellNum}(0-0+(1:minlen-1200)),':b','linewidth',3);
        
        
        ax3 = axis;
        axis([tStart-1 tEnd-1 ax3(3) ax3(4)/2])

end

legend('Black Box','Recorded','Cascade Conv');
grid on
set(gca,'fontsize',14);
xlabel('Time (sec)'); ylabel('Response (spikes/sec)');
% % set(gcf,'position',[   0.0063    0.2356    0.6861    0.3308]);
% set(gcf,'position',[ 0.0063    0.2354    0.7219    0.4549]);
