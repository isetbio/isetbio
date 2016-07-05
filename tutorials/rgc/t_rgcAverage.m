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
% ecc = 15.75; % mm = 2.5 degrees visual angle
% fov = 8;%0.5;

ecc = 0.4;
fov = 0.3;

experimentI   = 1;       % Choose dataset to load parameters and spikes
cellTypeI     = 3;       % Choose 1. OnPar, 2. OffPar, 3. OnMidg, 4. OffMidg, 5. SBC
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

%% Moving bar stimulus


% Set up Gabor stimulus using sceneCreate('harmonic',params)
% fov = 0.6;%1.2;

paramsStim.barwidth = 2;
paramsStim.meanLuminance = 200;
paramsStim.row = 64; params.col = 64;

paramsStim.expTime = 0.001;
paramsStim.timeInterval = 0.001;
paramsStim.nSteps = 100;     % Number of stimulus frames

upSampleFactor = 10;
paramsStim.timeInterval = .001;%(1/125)/upSampleFactor;%0.001; % sec
paramsStim.expTime      = .001;%(1/125)/upSampleFactor;%0.001; % sec

paramsStim.fov = fov;
paramsStim.radius = ecc;
paramsStim.theta = 0;  % 3 oclock on the retina
paramsStim.side = 'left';

iStim = ieStimulusBar(paramsStim);
sensor = iStim.absorptions;

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
switch cellTypeI
    case 1
        bpParams.cellType = 'onDiffuse';
    case 2
        bpParams.cellType = 'offDiffuse';
    case 3
        bpParams.cellType = 'onMidget';
    case 4
        bpParams.cellType = 'offMidget';
    case 5
        bpParams.cellType = 'onDiffuse';
end
% sets filter as theoretical, mean physiology, or individual phys:
bpParams.filterType = 1; 
% sets linear, on half-wave rectification, or on and off half-wave rect
bpParams.rectifyType = 1;
% bpParams.rectifyType = 3;

% bpParams.cellLocation = cellNumber;

bpParams.ecc = ecc;

bp = bipolar(osB, bpParams);

% bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
% bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);

% Need to fix bp compute bc it copies sensor
bp = bipolarCompute(bp, osB);

% bipolarPlot(bp,'response');

% bp = bpSet(bp,'responseCenter',sensor.data.volts);
%%
% Set parameters
clear params innerRetinaSU
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

params.inputScale = size(bp.sRFcenter,1);
params.inputSize = size(bp.responseCenter);

% Create object
innerRetinaSU = irPhys(bp, params);
nTrials = 30; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);
%% Compute the inner retina response

% Linear convolution
innerRetinaSU = irCompute(innerRetinaSU, bp); 

% innerRetinaSU = irComputeContinuous(innerRetinaSU, bp); 
% innerRetinaSU = irNormalize(innerRetinaSU, innerRetina);
% innerRetinaSU = irComputeSpikes(innerRetinaSU); 

% Get the PSTH from the object
innerRetinaSUPSTH = mosaicGet(innerRetinaSU.mosaic{1},'responsePsth');

figure; plot(vertcat(innerRetinaSUPSTH{:})')
title(sprintf('%s Simulated Mosaic at %1.1f\\circ Ecc\nMoving Bar Response',cellType(1:end-4),fov));
xlabel('Time (msec)'); ylabel('PSTH (spikes/sec)');
set(gca,'fontsize',14);
axis([0 length(innerRetinaSUPSTH{1}) 0 130]);
grid on;

%%
figure;
% 
cone_mosaic = sensor.human.coneType;
% [xg yg] = meshgrid([1:size(cone_mosaic,1),1:size(cone_mosaic,1)]);
% xg2 = xg(1:size(cone_mosaic,1),1:size(cone_mosaic,2)); yg2 = yg(1:size(cone_mosaic,1),1:size(cone_mosaic,2));
scale1 = 1;%size(cone_mosaic,1)/80; 
scale2 = 1;%size(cone_mosaic,2)/40; 
% % figure; scatter(xg2(:),yg2(:),40,4-cone_mosaic(:),'o','filled'); colormap jet; set(gca,'color',[0 0 0])
% % figure; 
% hold on;
% scatter(yg2(:)./scale1,xg2(:)./scale2,40,4-cone_mosaic(:),'o','filled'); colormap jet; set(gca,'color',[0 0 0])

% figure; 
imagesc(4-cone_mosaic); colormap jet;

    % circle samples
bpRad = 1;%0.6;
circle_samples = 0:0.05:2*pi;
x_circle = bpRad*cos(circle_samples);
y_circle = bpRad*sin(circle_samples);
% 
szbp = size(bp.sRFcenter,1)-.5;
% figure;
hold on;
szmosaic = size(cone_mosaic);
ycind = 0; xcind = 0;
for yc = 1:szbp*2:szmosaic(1)
    ycind = ycind+1;
    for xc = 1:szbp*2:szmosaic(2)
        xcind = xcind+1;
%         plot((szbp*y_circle+yc)./scale2,(szbp*x_circle+xc+(szbp)*mod(ycind+1,2))./scale1,'linewidth',1,'color','c');
        plot((szbp*x_circle+xc+0*(szbp)*mod(ycind+1,2))./scale1,(szbp*y_circle+yc)./scale2,'linewidth',1,'color','c');
    end
end


%
for i = 1:length(innerRetinaSU.mosaic{1}.cellLocation); 
    loc(i,:) = params.inputScale.*innerRetinaSU.mosaic{1}.cellLocation{i}; 
end;
% figure; scatter(loc(:,1),loc(:,2),100); axis equal

rfRad = 1*innerRetinaSU.mosaic{1}.rfDiaMagnitude/2;
circle_samples = 0:0.05:2*pi;
c1(1,:) = rfRad*cos(circle_samples);
c1(2,:) = rfRad*sin(circle_samples);

% figure; 
hold on;
for i = 1:length(innerRetinaSU.mosaic{1}.cellLocation)
%     plot(loc(i,2)+c1(1,2:end)+0*rfRad-0*innerRetinaSU.mosaic{1}.cellLocation{1,1}(1),loc(i,1)+c1(2,2:end)+0*rfRad-0*innerRetinaSU.mosaic{1}.cellLocation{1,1}(2),'m','linewidth',3); axis equal
        plot(loc(i,2)+c1(1,2:end)+1*rfRad-1*innerRetinaSU.mosaic{1}.cellLocation{1,1}(1),loc(i,1)+c1(2,2:end)+rfRad-2*innerRetinaSU.mosaic{1}.cellLocation{1,1}(2),'m','linewidth',3); axis equal
%     plot(loc(i,1)+c1(1,2:end)+2*rfRad-1*innerRetinaSU.mosaic{1}.cellLocation{1,1}(1),loc(i,2)+c1(2,2:end)+rfRad-1*innerRetinaSU.mosaic{1}.cellLocation{1,1}(2),'m','linewidth',3); axis equal
%     plot(loc(i,2)+c1(2,2:end)-16.5,loc(i,1)+c1(1,2:end)-16.5,'m','linewidth',3); axis equal
    
end
axis equal
% axis([0 75 0 75])
title(sprintf('%s simulated mosaic at %1.1f\\circ Ecc',cellType(1:end-3),ecc));
set(gca,'fontsize',16)

%%

mloc = max(loc);
psthMovie = zeros(2*ceil(mloc(1)),2*ceil(mloc(2)),length(vertcat(innerRetinaSUPSTH{:})'));
for i = 1:length(innerRetinaSU.mosaic{1}.cellLocation); 
    loc(i,:) = params.inputScale.*innerRetinaSU.mosaic{1}.cellLocation{i}; 
    psthMovie(round(2*loc(i,1)),round(2*loc(i,2)),:) = innerRetinaSUPSTH{i};
end;
ieMovie(psthMovie(:,:,200:end));
%%
% figure; 
% subplot(131);
% imagesc(psthMovie(:,:,200));
% title(sprintf('Mosaic PSTH at RGC Location\nt = 100 msec'));
% set(gca,'fontsize',16)
% subplot(132);
% imagesc(psthMovie(:,:,500));
% title(sprintf('Mosaic PSTH at RGC Location\nt = 400 msec'));
% set(gca,'fontsize',16)
% subplot(133);
% imagesc(psthMovie(:,:,800));
% title(sprintf('Mosaic PSTH at RGC Location\nt = 700 msec'));
% set(gca,'fontsize',16)