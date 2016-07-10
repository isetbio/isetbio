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


%% Initialize parameters for RGC mosaic

% Specify retinal patch eccentricity and field of view;
ecc = 0.03;    % in mm
fov = 3*0.25;  % in degrees of visual angle

% Choose the experimental dataset, cell type, and stimulus
experimentI   = 1;       % Choose 1. RPE 2016-02-17-1, loads parameters and spikes
cellTypeI     = 7;       % Choose 1. OnPar, 2. OffPar, 3. OnMidg, 4. OffMidg, 5. SBC
stimulusTestI = 1;       % Choose 1. Moving bar
%% Load image sequence

stimulusSelect = 1;
nSteps = 100;
barWidth = 40;
freq = 5;

switch stimulusSelect
    case 1
        % Grating subunit stimulus
        clear params
        stimP.fov      = fov; % degrees
        stimP.barWidth = barWidth;
        stimP.nSteps   = nSteps;
        stimP.expTime = 0.001; % sec
        stimP.timeInterval = 0.001; % sec
        stimP.freq = freq;
        stimP.radius = ecc;
        stimP.theta = 0;
        stimP.side = 'left';
        stimP.downSampleFactor = 8;
        iStim = ieStimulusGratingSubunit(stimP);

        % load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/iStim_subunitGrating.mat');
%         load('/Users/james/Documents/MATLAB/isetbio misc/iStim_subunitGratingLong.mat')
        absorptions = iStim.absorptions; % cone isomerizations
%         ieMovie(iStim.sceneRGB);
    case 2
        % Natural scene with eye movements stimulus
        rdt = RdtClient('isetbio');
        rdt.crp('resources/data/rgc');
        data = rdt.readArtifact('testmovie_schemeA_8pix_Identity_8pix', 'type', 'mat');
        testmovie = data.testmovie;
        absorptions = ieStimulusMovie(testmovie(:,:,1:nSteps)); 
end


sensor = absorptions;

% ieMovie(iStim.sceneRGB)

% View movie of cone absorptions
ieMovie(sensor.data.volts);


%% Outer segment calculation - biophysical model
% The next step is to get the outer segment current. The biophysical outer
% segment model is employed here.

% % Initialize
osB = osCreate('biophys');

% Set size of retinal patch based on absorptions sensor object
patchSize = sensorGet(sensor,'width','m');
osB = osSet(osB, 'patch size', patchSize);

% Set time step of simulation equal to absorptions sensor object
timeStep = sensorGet(sensor,'time interval','sec');
osB = osSet(osB, 'time step', timeStep);

% Set mean current of outer segment as boundary condition
sensorVolts = sensorGet(sensor,'volts');
paramsOS.bgVolts = 1*mean(sensorVolts(:));
paramsOS.ecc = ecc; % mm
clear sensorVolts

% Compute the outer segment response to the absorptions with the linear
% model.
osB = osCompute(osB,sensor,paramsOS);

% % Plot the photocurrent for a pixel.
% osPlot(osB,sensor);

% Subsample outer segment current
% osB.osSet('coneCurrentSignal',osB.coneCurrentSignal(:,:,1:8:end));

ieMovie(osB.coneCurrentSignal);

downSampleFactor = 8;
osCurrentSize = osGet(osB,'size');
osCurrentNew = zeros([osCurrentSize(1) osCurrentSize(2) osCurrentSize(3)/downSampleFactor]);
for fr = 1:osCurrentSize(3)/downSampleFactor
    osCurrentNew(:,:,fr) = osB.coneCurrentSignal(:,:,downSampleFactor*(fr-1)+downSampleFactor);
end

osBDS = osB;
osBDS = osSet(osBDS,'coneCurrentSignal', osCurrentNew);
osBDS = osSet(osBDS,'timeStep',0.008);
%% Find bipolar responses
% Build the bipolar mosaic and compute the response

clear bp

% Select the correct type of bipolar cell based on the user-specified RGC mosaic.
switch cellTypeI
    case {1,6}
        bpParams.cellType = 'onDiffuse';
    case {2,7}
        bpParams.cellType = 'offDiffuse';
    case 3
        bpParams.cellType = 'onMidget';
    case 4
        bpParams.cellType = 'offMidget';
    case 5
        bpParams.cellType = 'onDiffuse';
    otherwise
end

% Sets filter as theoretical, mean physiology, or individual phys:
bpParams.filterType = 1; 

% Sets linear, on half-wave rectification, or on and off half-wave rect
bpParams.rectifyType = 1;
% bpParams.rectifyType = 3;

% Make sure bipolar cells have the right eccentricity
bpParams.ecc = ecc;

% Initialize the bipolar object
bp = bipolar(osB, bpParams);

% Compute the bipolar mosaic response
bp = bipolarCompute(bp, osB);

% Visualize the bipolar mosaic response
% bipolarPlot(bp,'response');
ieMovie(bp.responseCenter);

%% Set RGC mosaic parameters

% Switch on the conditions indices
% Experimental dataset
switch experimentI
    case 1; experimentID = 'RPE_201602171';
    otherwise; error('Data not yet available');
end
% The other experimental data will be added to the RDT in the future.

% Stimulus: white noise or natural scene movie with eye movements
switch stimulusTestI
    case 1; stimulusTest = 'bar';
end

% Cell type: ON or OFF Parasol
switch cellTypeI
    case 1; 
        cellType = 'On Parasol RPE';      
    case 2; 
        cellType = 'Off Parasol RPE';        
    case 3; 
        cellType = 'On Midget RPE';        
    case 4; 
        cellType = 'Off Midget RPE';
    case 5; 
        cellType = 'SBC RPE';
    case 6;
        cellType = 'On Parasol Apricot';
        
    case 7;
        cellType = 'On Parasol Apricot';
        
    case 8;
        cellType = 'Off Parasol Apricot';
        
    case 9;
        cellType = 'On Midget Apricot';
        
    case 10;
        cellType = 'Off Midget Apricot';
        
    otherwise;
        cellType = 'SBC Apricot';
end

%% Set other RGC mosaic parameters

clear params innerRetinaSU
params.name = 'macaque phys';
params.eyeSide = 'left'; 
params.eyeRadius = ecc; 
params.fov = fov;
params.eyeAngle = 0; ntrials = 0;

% Determined at beginning by user
params.experimentID = experimentID; % Experimental dataset
params.stimulusTest = stimulusTest; % WN or NSEM
params.cellType = cellType;         % ON or OFF Parasol;

% Set flag for average mosaic
params.averageMosaic = 1;

% Tell the RGC mosaic about how many bipolars per cone
params.inputScale = size(bp.sRFcenter,1);
params.inputSize = size(bp.responseCenter);

% Create RGC object
innerRetinaSU = irPhys(bp, params);
nTrials = 30; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);

%% Plot the cone, bipolar and RGC mosaics

% mosaicPlot(innerRetinaSU,bp,sensor,params,cellType,ecc);

%% Compute the inner retina response

innerRetinaSU = irCompute(innerRetinaSU, bp); 

% Get the PSTH from the object
innerRetinaSUPSTH = mosaicGet(innerRetinaSU.mosaic{1},'responsePsth');

% Plot all of the PSTHs together
% figure; plot(vertcat(innerRetinaSUPSTH{:})')
figure; plot(vertcat(innerRetinaSUPSTH{21})')
title(sprintf('%s Simulated Mosaic at %1.1f\\circ Ecc\nMoving Bar Response',cellType(1:end-4),ecc));
xlabel('Time (msec)'); ylabel('PSTH (spikes/sec)');
set(gca,'fontsize',14);
% axis([0 length(innerRetinaSUPSTH{1}) 0 130]);
grid on;

%% Make a movie of the PSTH response

psthMovie = mosaicMovie(innerRetinaSUPSTH,innerRetinaSU, params);
