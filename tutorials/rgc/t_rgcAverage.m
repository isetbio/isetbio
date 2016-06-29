% t_rgcAverage
% 
% An experimental dataset that represents an RGC mosaic of a peripheral
% patch is loaded from Chichilnisky Lab data, and a mosaic of the average
% RGC is created that completely tiles the space. The cone spacing and RGC
% size can be scaled to an eccentricity input by the user.
% 
% The experimental recordings used are from the dataset RPE 201602171.
% 
% There is a dependency on the repository isetbio/EJLFigureReproduction.
% This must be added to the Matlab path.
% 
% See also:
%   EJLFigureReproduction/t_rgcNaturalScenesFig2.m
%   t_rgcCascade.m
%   t_rgcPeriphery.m
% 
% 6/2016 JRG (c) isetbio team

%% Initialize 
clear
% ieInit;

%% Switch on input type
% White noise (WN) or natural scenes with eye movements (NSEM)

% Radius of simulated patch
ecc = 15.75; % mm = 2.5 degrees visual angle
fov = 8;%0.5;

% ecc = 0.75;
% fov = 0.5;

experimentI   = 1;       % Choose dataset to load parameters and spikes
cellTypeI     = 1;       % Choose 1. OnPar, 2. OffPar, 3. OnMidg, 4. OffMidg, 5. SBC
stimulusTestI = 1;       % Choose WN test stimulus (1) or NSEM test stimulus (2)
    
% Switch on the conditions indices
% Experimental dataset
switch experimentI
    case 1; experimentID = 'RPE_201602171';
    otherwise; error('Data not yet available');
end
% The other experimental data will be added to the RDT in the future.

% Stimulus: white noise or natural scene movie with eye movements
switch stimulusTestI
    case 1; stimulusTest = 'WN';
    case 2; stimulusTest = 'NSEM';
end

% Cell type: ON or OFF Parasol
switch cellTypeI
    case 1; 
        cellType = 'On Parasol RPE';         
%         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/mosaicGLM_RPE_onPar.mat') 
    case 2; 
        cellType = 'Off Parasol RPE';              
%         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/mosaicGLM_RPE_offPar.mat') 
    case 3; 
        cellType = 'On Midget RPE';
%         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/mosaicGLM_RPE_onMid.mat') 
    case 4; 
        cellType = 'Off Midget RPE';
%         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/mosaicGLM_RPE_offMid.mat') 
    case 5; 
        cellType = 'SBC RPE';
%         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/mosaicGLM_RPE_onSBC.mat') 
end

% mosaicGLMaverage = mosaicAverage(mosaicGLM);

%% Load stimulus movie and fit/spiking data using RemoteDataToolbox

% Loads the appropriate movie and spiking data for the experimental
% conditions.
% [testmovie, xval_mosaic] =  loadDataRGCFigure2(experimentI,stimulusTestI,cellTypeI);
% Binary white noise test movie

load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/WN_testmovie_lnfit.mat')

% Length of WN movie is 1200, take nFrames to limit natural movie to same length
nFrames = 1; 
testmovieshort = WN_testmovie_lnfit(:,:,1:nFrames); 
 
upSampleFactor = 10;
testmovieRS = zeros(size(testmovieshort,1),size(testmovieshort,2),upSampleFactor*size(testmovieshort,3));
for frnum = 1:nFrames
    for frrep = 1:upSampleFactor
        testmovieRS(:,:,(frnum-1)*upSampleFactor+frrep) = testmovieshort(:,:,frnum);
    end
end
%% Show test movie
showFrames = 50;
% ieMovie(testmovieshort(:,:,1:showFrames));

%% Generate display, scene, oi, sensor
paramsStim.nsteps = nFrames;%size(testmovieshort,3);
%  Bipolar filter is setfor 0.001 sec, so it needs to be 0.001
% Downsample cone current signal to get to RGC sampling rate

paramsStim.timeInterval = (1/125)/upSampleFactor;%0.001; % sec
paramsStim.expTime      = (1/125)/upSampleFactor;%0.001; % sec

paramsStim.fov = fov;
paramsStim.radius = ecc;
paramsStim.theta = 0;  % 3 oclock on the retina
paramsStim.side = 'left';

% iStim = ieStimulusMovie(testmovieshort(:,:,1:nFrames),paramsStim);
% iStim = ieStimulusMovie(testmovieshort(:,:,1:nFrames),paramsStim);
iStim = ieStimulusMovie(permute(testmovieRS(:,:,1:nFrames), [2 1 3]),paramsStim);
sensor = iStim.sensor;

% sensor.data.volts = 5e-4*double(testmovieshort)./255;
% sensor.data.volts = 5e-4*testmovieRS./max(testmovieRS(:));
clear testmovieRS % testmovieshort
%% Outer segment calculation - biophysical model
% The iStim structure generates the movie, the scene, the oi and the
% cone absorptions. The next step is to get the outer segment current. The
% linear outer segment model is employed here.

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
paramsOS.bgVolts = 1*mean(sensorVolts(:));
paramsOS.ecc = ecc; % mm
clear sensorVolts

osBSub = osB;
% Compute the outer segment response to the absorptions with the linear
% model.
osB = osCompute(osB,sensor,paramsOS);

% % Plot the photocurrent for a pixel.
osPlot(osB,sensor);

% % osBSub.osSet('coneCurrentSignal',0);
% 
% osBSub.osSet('coneCurrentSignal',osB.coneCurrentSignal(:,:,1:80:end));
% clear osB

%% Find bipolar responses
clear bp
% os = osBSub;

bpParams.cellType = 'onDiffuse';
% sets filter as theoretical, mean physiology, or individual phys:
bpParams.filterType = 1; 
% sets linear, on half-wave rectification, or on and off half-wave rect
bpParams.rectifyType = 1;
% bpParams.rectifyType = 3;

% bpParams.cellLocation = cellNumber;

bp = bipolar(osB, bpParams);

% bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
% bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);

% Need to fix bp compute bc it copies sensor
bp = bipolarCompute(bp, osB);

% bipolarPlot(bp,'response');

% bp = bpSet(bp,'responseCenter',sensor.data.volts);
%%
% Set parameters
clear params
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = ecc; 
params.fov = fov;
params.eyeAngle = 0; ntrials = 0;

% Determined at beginning to allow looping
params.experimentID = experimentID; % Experimental dataset
params.stimulusTest = stimulusTest; % WN or NSEM
params.cellType = cellType;         % ON or OFF Parasol;

params.cellIndices = 10;

params.averageMosaic = 1;

% Create object
innerRetinaSU = irPhys(bp, params);
nTrials = 30; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);
%% Compute the inner retina response

% Linear convolution
innerRetinaSU = irCompute(innerRetinaSU, bp); 

% innerRetinaSU = irComputeContinuous(innerRetinaSU, bp); 
% innerRetinaSU = irNormalize(innerRetinaSU, innerRetina);
% innerRetinaSU = irComputeSpikes(innerRetinaSU); 

innerRetinaSU = irSet(innerRetinaSU,'timing',1/120);

% Get the PSTH from the object
innerRetinaSUPSTH = mosaicGet(innerRetinaSU.mosaic{1},'responsePsth');

figure; plot(vertcat(innerRetinaSUPSTH{:})')
xlabel('Time (msec)'); ylabel('PSTH (spikes/sec)');
set(gca,'fontsize',14);
axis([0 1000 0 130]);
grid on;

%%
figure;

cone_mosaic = sensor.human.coneType;
[xg yg] = meshgrid([1:size(cone_mosaic,1),1:size(cone_mosaic,1)]);
xg2 = xg(1:size(cone_mosaic,1),1:size(cone_mosaic,2)); yg2 = yg(1:size(cone_mosaic,1),1:size(cone_mosaic,2));
scale1 = 1;%size(cone_mosaic,1)/80; 
scale2 = 1;%size(cone_mosaic,2)/40; 
% figure; scatter(xg2(:),yg2(:),40,4-cone_mosaic(:),'o','filled'); colormap jet; set(gca,'color',[0 0 0])
% figure; 
hold on;
scatter(yg2(:)./scale1,xg2(:)./scale2,40,4-cone_mosaic(:),'o','filled'); colormap jet; set(gca,'color',[0 0 0])

%%
for i = 1:length(innerRetinaSU.mosaic{1}.cellLocation); loc(i,:) = innerRetinaSU.mosaic{1}.cellLocation{i}; end;
% figure; scatter(loc(:,1),loc(:,2),100); axis equal


c1 = contour(innerRetinaSU.mosaic{1}.sRFcenter{1},[-.19 -.19]);
% figure; 
hold on;
for i = 1:length(innerRetinaSU.mosaic{1}.cellLocation)

    plot(loc(i,2)+c1(2,2:end)-6.5,loc(i,1)+c1(1,2:end)-6.5,'m','linewidth',3); axis equal
    
end