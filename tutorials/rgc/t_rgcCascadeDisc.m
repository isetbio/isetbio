
%  Need to have run t_rgcCascade and keep it in memory to do comparison
%  plots at end.

for cellInd = 5% :10

    
    
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
load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/saved_iStim_linEqDisc/linEqDisc_10_notemp.mat')

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
sensorDisc.data.volts = 5e-4*ieScale(testmovieRS);
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

bpDisc = bipolar(osDisc);
% bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
% bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);

% bipolarThreshold = -40;
% bp = bipolarSet(bp,'threshold',bipolarThreshold);

% Need to fix bp compute bc it copies sensor
bpDisc = bipolarCompute(bpDisc, osDisc);

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
innerRetinaDisc = irPhys(bpDisc, params);
nTrials = 57; innerRetinaDisc = irSet(innerRetinaDisc,'numberTrials',nTrials);

% Compute the inner retina response
% 
% tCenterNew{1,1} = 1; tSurroundNew{1,1} = 0;
% innerRetina.mosaic{1}.mosaicSet('tCenter',tCenterNew);
% innerRetina.mosaic{1} = mosaicSet(innerRetina.mosaic{1},'tSurround',tSurroundNew);

% Linear convolution
% innerRetinaDisc = irComputeContinuous(innerRetinaDisc, bpDisc);
innerRetinaDisc = irCompute(innerRetinaDisc, bpDisc);
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
% irPlot(innerRetina,'linear','cell',[7 1]);
% 
% irPlot(innerRetina,'raster','cell',[7 1]);
% irPlot(innerRetina,'psth','cell',[7 1]);

% Spike computation
% for tr = 1:ntrials
%     innerRetina = irComputeSpikes(innerRetina, os1);
% end

% Get the PSTH from the object
innerRetinaPSTH = mosaicGet(innerRetina.mosaic{1},'responsePsth');



% innerRetinaRecordedPSTH = mosaicGet(innerRetinaRecorded.mosaic{1},'responsePsth');

end

%%

tStart = 1.5;
tEnd = 21;%1*8.5;
cellNum = cellInd;

vcNewGraphWin([],'upperleftbig'); 
subplot(411); hold on;
irPlot(innerRetina,'raster','cell',[cellNum 1],'hold','on','color','r')
title(sprintf('Black Box, NSEM, off parasol cell [%d 1]',cellNum));
set(gca,'fontsize',14);
axis([tStart tEnd 0 57]);
axis off

subplot(412); hold on;
irPlot(innerRetinaDisc,'raster','cell',[cellNum 1],'hold','on','color','b')
title(sprintf('Cascade Conv, NSEM, off parasol cell [%d  1]',cellNum));
set(gca,'fontsize',14);
% axis([tStart-.04 tEnd-.04 0 57]); % when using theoretical irGLM
axis([tStart tEnd 0 57]);
axis off


subplot(413); hold on;
irPlot(innerRetinaRecorded,'raster','cell',[cellNum 1],'hold','on','color','k')
title(sprintf('Recorded, NSEM, off parasol cell [%d  1]',cellNum));
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
% irPlot(innerRetinaDisc,'psth','cell',[cellNum 1],'hold','on','color','b')
% irPlot(innerRetinaRecorded,'psth','cell',[cellNum 1],'hold','on','color','g')

% % irPlot(innerRetinaRecorded,'psth','cell',[cellNum 1],'hold','on','color','g')
% legend('Cascade Conv','Black Box','Recorded');
% grid on
% ax2 = axis;
% axis([tStart/10 tEnd/10 ax2(3) ax2(4)]);
% set(gca,'fontsize',14);
% % set(gcf,'position',[   0.0063    0.2356    0.6861    0.3578]);
% set(gcf,'position',[ 0.0063    0.2354    0.7219    0.4549]);
 set(gcf,'position',[ 0.0063   -0.0444    0.8819    0.9378]);

%
% figure;

subplot(414)
% minlen = min([  length(innerRetinaDiscPSTH{cellNum}) ]);
minlen = min([length(innerRetinaPSTH{cellNum}) length(innerRetinaRecordedPSTH{cellNum}) length(innerRetinaDiscPSTH{cellNum}) ]);
hold off
switch stimulusTestI
    case 1
        plot((00+[1:minlen-1200])./1208, innerRetinaPSTH{cellNum}(600+(1:minlen-1200)),'r','linewidth',3);
        
        hold on;
        plot([1:minlen-1200]./1208,innerRetinaRecordedPSTH{cellNum}((0+(1:minlen-1200))),':k','linewidth',2);
        
        plot((00+[1:minlen-1200])./1208, innerRetinaDiscPSTH{cellNum}(600-0+(1:minlen-1200)),':b','linewidth',3);
        
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
%         plot((00+[1:minlen-1200])./1208, innerRetinaDiscPSTH{cellNum}(000-36+(1:minlen-1200)),':b','linewidth',3);

%         plot((00+[1:minlen-1200])./1208, innerRetinaPSTH{cellNum}(600+(1:minlen-1200)),'r','linewidth',3);
        
        plot((00+[1:minlen-1200])./1208, innerRetinaPSTH{cellNum}(1200+(1:minlen-1200)),'r','linewidth',3);
        hold on;
        plot([1:minlen-1200]./1208,innerRetinaRecordedPSTH{cellNum}((0+(1:minlen-1200))),':k','linewidth',2);
        hold on;
        plot((00+[1:minlen-1200])./1208, innerRetinaDiscPSTH{cellNum}(1200-0+(1:minlen-1200)),':b','linewidth',3);
        
        
        ax3 = axis;
        axis([tStart-1 tEnd-1 ax3(3) ax3(4)/2])
%         axis([tStart-1 tEnd-1 ax3(3) 100])
        
        % axis([0-.5 8-.5 0 100])
end

legend('Black Box','Recorded','Cascade Conv');
grid on
set(gca,'fontsize',14);
xlabel('Time (sec)'); ylabel('Response (spikes/sec)');
% % set(gcf,'position',[   0.0063    0.2356    0.6861    0.3578]);
% set(gcf,'position',[ 0.0063    0.2354    0.7219    0.4549]);

%%

savename1 = 'innerRetinaDiscLinear5';
save(['/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/innerRetina/' savename1 '.mat'],'innerRetinaDisc');


