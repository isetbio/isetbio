% t_rgcCascade
%
% A comparison of RGC firing from experimental recordings, "black box" GLM
% simulations and "cascade model" GLM simulations. The experimental
% recordings and the black box simulations are found in Brackbill, Heitman,
% ..., Chichilnisky, J. Neurophys 2016. 
% 
% This demonstrates that the linear cascade model matches the output of the
% black box model. The cascade model allows us to turn on nonlinearities in
% the outer segment current computation as well as spatial bipolar
% subunits. We can quantify whether these nonlinearities improve the
% prediction of the firing in response to natural scene movies and compare
% them with the experimental recordings.
% 
% See also EJLFigureReproduction/t_rgcNaturalScenesFig2.m
% 
% 6/2016 JRG (c) isetbio team

%% Initialize 
clear
% ieInit;

%% Switch on input type
% White noise (WN) or natural scenes with eye movements (NSEM)

% In order to generate the plot of the fractional variance for all of the
% selected cells, the script must be run with the three for loops (lines
% 46-48), and the two stimulusTestI values must be iterated through. To run the
% tutorial step by step, uncomment lines 40-42 and execute each cell,
% instead of running the script through.


plotFracFlag = 0;

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
%% Load stimulus movie and fit/spiking data using RemoteDataToolbox

% Loads the appropriate movie and spiking data for the experimental
% conditions.
[testmovie, xval_mosaic] =  loadDataRGCFigure2(experimentI,stimulusTestI,cellTypeI);

% Length of WN movie is 1200, take nFrames to limit natural movie to same length
nFrames = 3600; 
testmovieshort = double(testmovie.matrix(:,:,1:nFrames)); 

frRS = 8;
testmovieRS = zeros(size(testmovieshort,1),size(testmovieshort,2),frRS*size(testmovieshort,3));
for frnum = 1:nFrames
    for frrep = 1:frRS
        testmovieRS(:,:,(frnum-1)*frRS+frrep) = testmovieshort(:,:,frnum);
    end
end


% range = max(testmovieRS(:)) - min(testmovieRS(:));
% testmovieRS = 1*(testmovieRS./range - mean(testmovieRS(:))/range);
% 
% 
% range = max(testmovieshort(:)) - min(testmovieshort(:));
% testmovieRS = 1*(testmovieshort./range - mean(testmovieshort(:))/range);

% testmovieshort = zeros(80,40,10);
% testmovieshort(:,20:30,:) = 1;
%% Show test movie
showFrames = 50;
% ieMovie(testmovieshort(:,:,1:showFrames));

%% Generate display, scene, oi, sensor
paramsStim.nsteps = 1;%size(testmovieshort,3);
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

iStim = ieStimulusMovie(testmovieshort(:,:,1),paramsStim);
sensor = iStim.sensor;

% sensor.data.volts = 5e-4*double(testmovieshort)./255;
sensor.data.volts = 5e-4*testmovieRS./max(testmovieRS(:));
% clear testmovieRS testmovieshort
%% Outer segment calculation - linear model
% The iStim structure generates the movie, the scene, the oi and the
% cone absorptions. The next step is to get the outer segment current. The
% linear outer segment model is employed here.

% Initialize
osL = osCreate('linear');
% osL = osSet(osL,'sConeFilter',1);
% osL = osSet(osL,'mConeFilter',1);
% osL = osSet(osL,'lConeFilter',1);
% % % % % Set eccentricity of os patch here

% Set size of retinal patch based on absorptions sensor object
patchSize = sensorGet(sensor,'width','m');
osL = osSet(osL, 'patch size', patchSize);

% Set time step of simulation equal to absorptions sensor object
timeStep = sensorGet(sensor,'time interval','sec');
osL = osSet(osL, 'time step', timeStep);

% Set circular convolution, only steady state
paramsOSL.convolutionType = 1; 

osLSub = osL;

% Compute the outer segment response to the absorptions with the linear
% model.
osL = osCompute(osL,sensor,paramsOSL);
% osL = osSet(osL,'coneCurrentSignal',sensor.data.volts);
% % Plot the photocurrent for a pixel.
% osPlot(osL,sensor);

osLSub = osSet(osL, 'time step', 8*timeStep);
osLSub.osSet('coneCurrentSignal',osLSub.coneCurrentSignal(:,:,1:8:end));
% osPlot(osLSub,sensor)
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

osBSub = osB;
% Compute the outer segment response to the absorptions with the linear
% model.
% osB = osCompute(osB,sensor,paramsOS);
% 
% % % Plot the photocurrent for a pixel.
% % osPlot(osB,sensor);
% 
% % osBSub.osSet('coneCurrentSignal',0);
% 
% osBSub.osSet('coneCurrentSignal',osB.coneCurrentSignal(:,:,1:80:end));
% clear osB
%% Find bipolar responses
clear bp os
os = osL;
% os = osLSub;
% os = osBSub;

% bp = bipolar(osL);

bp = bipolar(os);
% bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
% bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);

% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

% Need to fix bp compute bc it copies sensor
bp = bipolarCompute(bp, os);

% bipolarPlot(bp,'response');

% bp.responseCenter = sensor.data.volts;
% bp = bpSet(bp,'responseCenter',sensor.data.volts);
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
innerRetinaSU = irPhys(bp, params);
nTrials = 57; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);
%% Compute the inner retina response
% 
% tCenterNew{1,1} = 1; tSurroundNew{1,1} = 0;
% innerRetina.mosaic{1}.mosaicSet('tCenter',tCenterNew);
% innerRetina.mosaic{1} = mosaicSet(innerRetina.mosaic{1},'tSurround',tSurroundNew);

% Linear convolution
% NEED TO DO MEAN SCALING AND CENTERING LIKE RGC
innerRetinaSU = irCompute(innerRetinaSU, bp);
innerRetinaSU = irSet(innerRetinaSU,'timing',.008);

irPlot(innerRetinaSU,'linear','cell',[7 1]);
irPlot(innerRetinaSU,'raster','cell',[7 1]);
irPlot(innerRetinaSU,'psth','cell',[7 1]);

% rLinSU = mosaicGet(innerRetinaSU.mosaic{1},'responseLinear');
% rLinearSU{1,1,1} = 5*rLinSU{1,1}./max(rLinSU{1,1}) + 2.1;
% innerRetinaSU.mosaic{1}.mosaicSet('responseLinear', rLinearSU);

% Spike computation
% for tr = 1:nTrials%ntrials
%     innerRetinaSU = irComputeSpikes(innerRetinaSU, bp);
% end

irPlot(innerRetinaSU,'linear','cell',[7 1]);

% Get the PSTH from the object
innerRetinaSUPSTH = mosaicGet(innerRetinaSU.mosaic{1},'responsePsth');

%%
% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 
%% Generate outer segment object

% In this case, the RGC GLM calculation converts from the frame buffer
% values in the movie to the spiking responses.  For this calculation, we
% store the movie stimulus in the the outer segment object 'displayRGB'.

os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep', 1/120);

% Attach the movie to the object
os1 = osSet(os1, 'rgbData', double(testmovieshort));

% 
% os1 = osSet(os1, 'timeStep', .001);
% 
% % Attach the movie to the object
% os1 = osSet(os1, 'rgbData', double(testmovieshort));

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
innerRetina = irPhys(os1, params);
nTrials = 57; innerRetina = irSet(innerRetina,'numberTrials',nTrials);

% Plot a few simple properties

% % Select the cells used in the actual paper
% switch cellTypeI
%     case 1; cellInd = 2;
%     case 2; cellInd = 31;
% end
% 
% % Plot the spatial RF, temporal IR and post-spike filter
% if plotFracFlag == 0
%     irPlotFig2Linear(innerRetina,cellInd);
% end
%% Compute the inner retina response

% Linear convolution
% innerRetina = irComputeContinuous(innerRetina, os1);

innerRetina = irCompute(innerRetina, os1);
irPlot(innerRetina,'linear','cell',[7 1]);

irPlot(innerRetina,'raster','cell',[7 1]);
irPlot(innerRetina,'psth','cell',[7 1]);

% Spike computation
for tr = 1:ntrials
    innerRetina = irComputeSpikes(innerRetina, os1);
end

% Get the PSTH from the object
innerRetinaPSTH = mosaicGet(innerRetina.mosaic{1},'responsePsth');

%% Create a new inner retina object and attach the recorded spikes
% We also want to compare the spikes recorded in the experiment to those
% from the simulation. We create a separate inner retina object with
% isetbio to store these spikes. This makes it easy to use isetbio to plot
% aspects of the response.

% Create the object.
innerRetinaRecorded = irPhys(os1, params);  

% innerRetinaRecorded = irSet(innerRetinaRecorded,'timing',.008);
innerRetinaRecorded = irSet(innerRetinaRecorded,'numberTrials',nTrials);

% Set the recorded spikes that we got from the RDT into the object.
innerRetinaRecorded = irSet(innerRetinaRecorded,'recordedSpikes',xval_mosaic);
% Get the PSTH using an isetbio routine.
innerRetinaRecordedPSTH = mosaicGet(innerRetinaRecorded.mosaic{1},'responsePsth');

irPlot(innerRetinaRecorded,'raster','cell',[7 1]);


% Compare Rasters and PSTHs for a particular cell
% % Chnage the cellInd value to plot a different cell from the dataset
% switch cellTypeI
%     case 1; cellInd = 2;
%     case 2; cellInd = 31;
% end
% % Plot the simulated and recorded rasters for a cell
% irPlotFig2Raster(innerRetina, innerRetinaRecorded,cellInd,stimulusTestI);
% % Plot the simulated and recorded PSTHs for a cell
% irPlotFig2PSTH(innerRetina, innerRetinaPSTH, innerRetinaRecordedPSTH,cellInd,stimulusTestI);


%%


%%
tStart = 0.5;
tEnd = 29;%1*8.5;
cellNum = 7;

vcNewGraphWin([],'upperleftbig'); 
subplot(411); hold on;
irPlot(innerRetina,'raster','cell',[cellNum 1],'hold','on','color','r')
title(sprintf('Black Box, WN, off parasol cell [%d 1]',cellNum));
set(gca,'fontsize',14);
axis([tStart tEnd 0 57]);
axis off

subplot(413); hold on;
irPlot(innerRetinaSU,'raster','cell',[cellNum 1],'hold','on','color','b')
title(sprintf('Cascade Conv, WN, off parasol cell [%d  1]',cellNum));
set(gca,'fontsize',14);
% axis([tStart-.04 tEnd-.04 0 57]); % when using theoretical irGLM
axis([tStart tEnd 0 57]);
axis off


subplot(412); hold on;
irPlot(innerRetinaRecorded,'raster','cell',[cellNum 1],'hold','on','color','k')
title(sprintf('Recorded, WN, off parasol cell [%d  1]',cellNum));
set(gca,'fontsize',14);

switch stimulusTestI
    case 1
        axis([tStart-0.5 tEnd-0.5 0 57]);
    case 2
        axis([tStart-1 tEnd-1 0 57]);
end
% axis off

% 
% subplot(414); 
% irPlot(innerRetina,'psth','cell',[cellNum 1],'hold','on','color','r')
% irPlot(innerRetinaSU,'psth','cell',[cellNum 1],'hold','on','color','b')
% irPlot(innerRetinaRecorded,'psth','cell',[cellNum 1],'hold','on','color','g')

% % irPlot(innerRetinaRecorded,'psth','cell',[cellNum 1],'hold','on','color','g')
% legend('Cascade Conv','Black Box','Recorded');
% grid on
% ax2 = axis;
% axis([tStart/10 tEnd/10 ax2(3) ax2(4)]);
% set(gca,'fontsize',14);
% % set(gcf,'position',[   0.0063    0.2356    0.6861    0.3578]);
% set(gcf,'position',[ 0.0063    0.2354    0.7219    0.4549]);

%
% figure;

subplot(414)
minlen = min([length(innerRetinaPSTH{cellNum}) length(innerRetinaRecordedPSTH{cellNum}) length(innerRetinaSUPSTH{cellNum}) ]);
hold off
switch stimulusTestI
    case 1
        plot((00+[1:minlen-1200])./1208, innerRetinaPSTH{cellNum}(600+(1:minlen-1200)),'r','linewidth',3);
        
        hold on;
        plot([1:minlen-1200]./1208,innerRetinaRecordedPSTH{cellNum}((0+(1:minlen-1200))),':k','linewidth',2);
        
        plot((00+[1:minlen-1200])./1208, innerRetinaSUPSTH{cellNum}(600-0+(1:minlen-1200)),':b','linewidth',3);
        
        ax3 = axis;
        axis([0 8.5 ax3(3) ax3(4)])

    case 2
%         plot((00+[1:minlen-1200])./1208, psthSim{cellNum}(1200+(1:minlen-1200)),'r','linewidth',3);
%         
%         hold on;
%         plot([1:minlen-1200]./1208,psthRecorded{cellNum}((1:minlen-1200)),':k','linewidth',2);

%         minlen = minlen - 1000;
%         plot((00+[1:minlen-1200])./1208, innerRetinaPSTH{cellNum}(1200+(1:minlen-1200)),'r','linewidth',3);
%         
%         hold on;
%         plot([1:minlen-1200]./1208,innerRetinaRecordedPSTH{cellNum}((000+(1:minlen-1200))),':k','linewidth',2);
%         
%         plot((00+[1:minlen-1200])./1208, innerRetinaSUPSTH{cellNum}(000-36+(1:minlen-1200)),':b','linewidth',3);

%         plot((00+[1:minlen-1200])./1208, innerRetinaPSTH{cellNum}(600+(1:minlen-1200)),'r','linewidth',3);
        
        plot((00+[1:minlen-1200])./1208, innerRetinaPSTH{cellNum}(1200+(1:minlen-1200)),'r','linewidth',3);
        hold on;
        plot([1:minlen-1200]./1208,innerRetinaRecordedPSTH{cellNum}((0+(1:minlen-1200))),':k','linewidth',2);
        hold on;
        plot((00+[1:minlen-1200])./1208, innerRetinaSUPSTH{cellNum}(1200-0+(1:minlen-1200)),':b','linewidth',3);
        
        
        ax3 = axis;
%         axis([tStart-1 tEnd-1 ax3(3) ax3(4)/4])
        axis([tStart-1 tEnd-1 ax3(3) 100])
        
        % axis([0-.5 8-.5 0 100])
end

legend('Black Box','Recorded','Cascade Conv');
grid on
set(gca,'fontsize',14);
xlabel('Time (sec)'); ylabel('Response (spikes/sec)');
% % set(gcf,'position',[   0.0063    0.2356    0.6861    0.3578]);
% set(gcf,'position',[ 0.0063    0.2354    0.7219    0.4549]);


%% Calculate fractional variance predicted
% The fraction of explained variance is used as a mesaure of the accuracy
% of the simulated PSTH.
fractionalVariance{experimentI,stimulusTestI,cellTypeI} = ...
    calculateFractionalVariance(innerRetinaPSTH, innerRetinaRecordedPSTH, stimulusTestI);

% Plot the explained variance for WN and NSEM against each other
if (stimulusTestI == 2) && (plotFracFlag == 1); 
    irPlotFig2FracVar(experimentI,cellTypeI,fractionalVariance); 
    
end;

plotFracFlag = 1;
%%%
