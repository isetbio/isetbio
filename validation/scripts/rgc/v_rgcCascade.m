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
% There is a dependency on the repository isetbio/EJLFigureReproduction.
% This must be added to the Matlab path.
% 
% See also EJLFigureReproduction/t_rgcNaturalScenesFig2.m
% 
% 6/2016 JRG (c) isetbio team

%% Initialize 
clear
% ieInit;

%% Switch on input type
% White noise (WN) or natural scenes with eye movements (NSEM)

experimentI   = 1;       % Choose dataset to load parameters and spikes
cellTypeI     = 2;%:2    % Choose On Parasol (1) or Off Parasol (2)
stimulusTestI = 1;%:2     % Choose WN test stimulus (1) or NSEM test stimulus (2)
    
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

% Load OS from RDT (1) or run from scratch (0)
loadOS = 0;
%% Load stimulus movie and fit/spiking data using RemoteDataToolbox

% Loads the appropriate movie and spiking data for the experimental
% conditions.
[testmovie, xval_mosaic] =  loadDataRGCFigure2(experimentI,stimulusTestI,cellTypeI);


% Length of WN movie is 1200, take nFrames to limit natural movie to same length
nFrames = 1200; 
testmovieshort = double(testmovie.matrix(:,:,1:nFrames)); 

if ~loadOS 
    
frRS = 8;
testmovieRS = zeros(size(testmovieshort,1),size(testmovieshort,2),frRS*size(testmovieshort,3));
for frnum = 1:nFrames
    for frrep = 1:frRS
        testmovieRS(:,:,(frnum-1)*frRS+frrep) = testmovieshort(:,:,frnum);
    end
end
%% Show test movie
showFrames = 50;
% ieMovie(testmovieshort(:,:,1:showFrames));

%% Generate display, scene, oi, sensor
paramsStim.nsteps = 1;%nFrames;%size(testmovieshort,3);
%  Bipolar filter is setfor 0.001 sec, so it needs to be 0.001
paramsStim.timeInterval = 1/125;%0.001; % sec
paramsStim.expTime = 1/125;%0.001; % sec

% % For 2013-08-19-6
% r = 12 mm = 36 degs;
% theta = 330 degs;
% side = 'left';
paramsStim.fov = 8;
paramsStim.radius = 36/3;
paramsStim.theta = 330;
paramsStim.side = 'left';

% iStim = ieStimulusMovie(testmovieshort(:,:,1:nFrames),paramsStim);
iStim = ieStimulusMovie(testmovieshort(:,:,1:nFrames),paramsStim);
sensor = iStim.sensor;

sensor.data.volts = 5e-4*double(testmovieshort)./255;
% sensor.data.volts = 5e-4*testmovieRS./max(testmovieRS(:));
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

% osLSub = osSet(osL, 'time step', 8*timeStep);
% osLSub.osSet('coneCurrentSignal',osL.coneCurrentSignal(:,:,1:8:end));
% osLSub.osSet('coneCurrentSignal',sensor.data.volts);


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
clear sensorVolts

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

else
    rdt = RdtClient('isetbio');
    rdt.crp('resources/data/rgc');
    data = rdt.readArtifact('full_osLSub_sensor_NSEM', 'type', 'mat');
    osLSub = data.osLSub;
    
end
tic

%% Generate outer segment object for GLM from RGB scene data

% In this case, the RGC GLM calculation converts from the frame buffer
% values in the movie to the spiking responses.  For this calculation, we
% store the movie stimulus in the the outer segment object 'displayRGB'.

os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep', 1/120);

% Attach the movie to the object
os1 = osSet(os1, 'rgbData', double(testmovieshort));

%% Set RGC parameters
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = 12; 
params.eyeAngle = 0; ntrials = 0;

% Determined at beginning to allow looping
params.experimentID = experimentID; % Experimental dataset
params.stimulusTest = stimulusTest; % WN or NSEM
params.cellType = cellType;         % ON or OFF Parasol

%% Create a new inner retina object and attach the recorded spikes
% We also want to compare the spikes recorded in the experiment to those
% from the simulation. We create a separate inner retina object with
% isetbio to store these spikes. This makes it easy to use isetbio to plot
% aspects of the response.

% Create the object.
innerRetinaRecorded = irPhys(os1, params);  
nTrials = 57; 
% innerRetinaRecorded = irSet(innerRetinaRecorded,'timing',.008);
innerRetinaRecorded = irSet(innerRetinaRecorded,'numberTrials',nTrials);

% Set the recorded spikes that we got from the RDT into the object.
innerRetinaRecorded = irSet(innerRetinaRecorded,'recordedSpikes',xval_mosaic);
% Get the PSTH using an isetbio routine.
innerRetinaRecordedPSTH = mosaicGet(innerRetinaRecorded.mosaic{1},'responsePsth');


for cellNumber = 1:118

%% Find bipolar responses
clear bp os
% os = osL;
os = osLSub;
% os = osBSub;

bpParams.cellType = 'offDiffuse';
% sets filter as theoretical, mean physiology, or individual phys:
bpParams.filterType = 3; 
% sets linear, on half-wave rectification, or on and off half-wave rect
bpParams.rectifyType = 1;
% bpParams.rectifyType = 3;

bpParams.cellLocation = cellNumber;

bp = bipolar(os, bpParams);

% bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
% bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);

% Need to fix bp compute bc it copies sensor
bp = bipolarCompute(bp, os);

% bipolarPlot(bp,'response');

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

params.cellIndices = cellNumber;

% Create object
innerRetinaSU = irPhys(bp, params);
nTrials = 57; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);
%% Compute the inner retina response

% Linear convolution
innerRetinaSU = irCompute(innerRetinaSU, bp); 

% innerRetinaSU = irComputeContinuous(innerRetinaSU, bp); 
% innerRetinaSU = irNormalize(innerRetinaSU, innerRetina);
% innerRetinaSU = irComputeSpikes(innerRetinaSU); 

innerRetinaSU = irSet(innerRetinaSU,'timing',1/120);

% Get the PSTH from the object
innerRetinaSUPSTH = mosaicGet(innerRetinaSU.mosaic{1},'responsePsth');



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

params.cellIndices = cellNumber;

% Create object
innerRetina = irPhys(os1, params);
nTrials = 57; innerRetina = irSet(innerRetina,'numberTrials',nTrials);

%% Compute the inner retina response

innerRetina = irCompute(innerRetina, os1);

% Get the PSTH from the object
innerRetinaPSTH = mosaicGet(innerRetina.mosaic{1},'responsePsth');


%% Plot the rasters

% % Set the time and cell number
% tStart = 1.5;% 9%1.5;
% tEnd = 21;%18%21;%1*8.5;
% cellNum = 1;
% 
% % Plot the original GLM prediction
% vcNewGraphWin([],'upperleftbig'); 
% subplot(312); hold on;
% % subplot(211); hold on;
% irPlot(innerRetina,'raster','cell',[cellNum 1],'hold','on','color','r')
% title(sprintf('Black Box, NSEM, off parasol cell [%d 1]',cellNum));
% set(gca,'fontsize',14);
% axis([tStart tEnd 0 57]);
% axis off
% 
% % Plot the biophys/subunit prediction
% subplot(313); hold on;
% % subplot(212); hold on;
% irPlot(innerRetinaSU,'raster','cell',[cellNum 1],'hold','on','color','b')
% title(sprintf('Cascade Conv, NSEM, off parasol cell [%d  1]',cellNum));
% set(gca,'fontsize',14);
% % axis([tStart-.04 tEnd-.04 0 57]); % when using theoretical irGLM
% % axis([tStart tEnd 0 57]);
% % switch stimulusTestI
% %     case 1
% %         axis([tStart+1 tEnd+1 0 57]);
% %     case 2
%         axis([tStart tEnd 0 57]);
% % end
% axis off
% 
% % % Plot the recorded spikes
% subplot(311); hold on;
% irPlot(innerRetinaRecorded,'raster','cell',[cellNum 1],'hold','on','color','k')
% title(sprintf('Recorded, NSEM, off parasol cell [%d  1]',cellNum));
% set(gca,'fontsize',14);
% 
% switch stimulusTestI
%     case 1
%         axis([tStart-0.5 tEnd-0.5 0 57]);
%     case 2
%         axis([tStart-1 tEnd-1 0 57]);
% end
%  set(gcf,'position',[ 0.0063   -0.0444    0.8819    0.9378]);

 %%
 
%  calculateFractionalVariance(innerRetinaPSTH, innerRetinaRecordedPSTH, stimulusTestI);
 
 %%
 
figure;
cellNum = cellNumber;
% subplot(414)
minlen = min([length(innerRetinaPSTH{1}) length(innerRetinaRecordedPSTH{cellNum}) length(innerRetinaSUPSTH{1 }) ]);
hold off
switch stimulusTestI
    case 1
        plot((00+[1:minlen-1200])./1208, innerRetinaPSTH{1}(600+(1:minlen-1200)),'r','linewidth',3);
        
        hold on;
        plot([1:minlen-1200]./1208,innerRetinaRecordedPSTH{cellNum}((0+(1:minlen-1200))),':k','linewidth',2);
        
        plot((00+[1:minlen-1200])./1208, innerRetinaSUPSTH{1}(600-0+(1:minlen-1200)),':b','linewidth',3);
        
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
        
        plot((00+[1:minlen-1200])./1208, innerRetinaPSTH{1}(1200+(1:minlen-1200)),'r','linewidth',3);
        hold on;
        plot([1:minlen-1200]./1208,innerRetinaRecordedPSTH{cellNum}((0+(1:minlen-1200))),':k','linewidth',2);
        hold on;
        plot((00+[1:minlen-1200])./1208, innerRetinaSUPSTH{1}(1200-0+(1:minlen-1200)),':b','linewidth',3);
        
        
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
% figure; scatter(innerRetinaPSTH{cellNum}(1200+(1:minlen-1200)),innerRetinaRecordedPSTH{cellNum}((0+(1:minlen-1200))))
% figure; scatter(innerRetinaPSTH{cellNum}(1200+(1:minlen-1200)),innerRetinaSUPSTH{cellNum}(1200-0+(1:minlen-1200)))
 
for cellNum2 = 1%length(innerRetinaPSTH)
    minlen = min([length(innerRetinaPSTH{cellNum2}) length(innerRetinaRecordedPSTH{cellNum2}) length(innerRetinaSUPSTH{cellNum2}) ]);
switch stimulusTestI
    

    case 1
        
        rsim = innerRetinaSUPSTH{cellNum2}(600+(1:minlen-1200));
        rrec = innerRetinaRecordedPSTH{cellNumber}((0+(1:minlen-1200)));
%         rrec = innerRetinaSUPSTH{cellNum2}(600-0+(1:minlen-1200));
    case 2
        rsim = innerRetinaSUPSTH{cellNum2}(1200+(1:minlen-1200));
        % rrec = innerRetinaSUPSTH{cellNum2}(1200-0+(1:minlen-1200));
        rrec = innerRetinaRecordedPSTH{cellNumber}((0+(1:minlen-1200)));
end

fracVar(cellNum2) = 1 - sum((rsim-rrec).^2)/sum((rrec-mean(rrec)).^2);
switch stimulusTestI
    case 1
        
        rsim = innerRetinaPSTH{cellNum2}(600+(1:minlen-1200));
        rrec = innerRetinaRecordedPSTH{cellNumber}((0+(1:minlen-1200)));
%         rrec = innerRetinaSUPSTH{cellNum2}(600-0+(1:minlen-1200));
    case 2
        rsim = innerRetinaPSTH{cellNum2}(1200+(1:minlen-1200));
        % rrec = innerRetinaSUPSTH{cellNum2}(1200-0+(1:minlen-1200));
        rrec = innerRetinaRecordedPSTH{cellNumber}((0+(1:minlen-1200)));
end

fracVar2(cellNum2) = 1 - sum((rsim-rrec).^2)/sum((rrec-mean(rrec)).^2);

switch stimulusTestI
    case 1
        
        rsim = innerRetinaPSTH{cellNum2}(600+(1:minlen-1200));
        rrec = innerRetinaSUPSTH{cellNum2}((600+(1:minlen-1200)));
%         rrec = innerRetinaSUPSTH{cellNum2}(600-0+(1:minlen-1200));
    case 2
        rsim = innerRetinaPSTH{cellNum2}(0+(1:minlen-1200));
        % rrec = innerRetinaSUPSTH{cellNum2}(1200-0+(1:minlen-1200));
        rrec = innerRetinaSUPSTH{cellNum2}((0+(1:minlen-1200)));
end

rsim = innerRetinaPSTH{cellNum2}(600+(1:minlen-1200));
rrec = innerRetinaSUPSTH{cellNum2}(600+(1:minlen-1200));
fracVar3(cellNum2) = 1 - sum((rsim-rrec).^2)/sum((rrec-mean(rrec)).^2);

% figure;scatter(innerRetinaPSTH{cellNum2}(600+(1:minlen-1200)),innerRetinaSUPSTH{cellNum2}(600+(1:minlen-1200)))
% 
% figure;scatter(innerRetinaRecordedPSTH{cellNum2}(0+(1:minlen-1200)),innerRetinaSUPSTH{cellNum2}(600+(1:minlen-1200)))


end
% fracVar
% fracVar2
% figure; scatter(fracVar, fracVar2); 
% hold on; plot(.01:.01:1,.01:.01:1)
% xlabel('SU'); ylabel('GLM');


fracVar1Cascade(cellNumber) = fracVar(cellNum2);

fracVar2Cascade(cellNumber) = fracVar2(cellNum2);
fracVarCascade(cellNumber) = fracVar3(cellNum2)

end
toc

figure; scatter(fracVar1Cascade,fracVar2Cascade);