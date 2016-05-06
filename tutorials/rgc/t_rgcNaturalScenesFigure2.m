% t_rgcNaturalScenesFigure2
% 
% Reproduce Fig. 2 of Heitman, Brackbill, Greschner, Litke, Sher &
% Chichilnisky, 2016 with isetbio.
% 
% http://biorxiv.org/content/early/2016/03/24/045336
% 
% Load a movie stimulus, generate an os object for the movie, generate an
% inner retina object, load parameters from a physiology experiment in the
% Chichilnisky Lab into the IR object and compute the response of the RGC
% mosaic.
% 
% Selected rasters and PSTHs for particular cells are shown which match the
% cells in the figure. The fractional variance explained is computed for
% each cell for a white noise test movie and a natural scenes test movie.
% These values are plotted against each other to show the failure of the
% GLM model to predict the cells' responses to natural scene stimuli.
% 
% Currently, only the dataset from the first experiment is available on the
% Remote Data Toolbox. This dataset includes parameters for GLM fits to
% cells in the On Parasol and Off Parasol mosaics as well as the recorded
% spikes from the experiment. 
% 
% 4/2016
% (c) isetbio team

%% Initialize 
clear
ieInit;

%% Switch on input type
% White noise (WN) or natural scenes with eye movements (NSEM)

% experimentI   = 1;
% cellTypeI     = 1;
% stimulusTestI = 1;

plotFracFlag = 0;

for experimentI   = 1       % Choose dataset to load parameters and spikes
for cellTypeI     = 1%:2    % Choose On Parasol (1) or Off Parasol (2)
for stimulusTestI = 1:2     % Choose WN test stimulus (1) or NSEM test stimulus (2)
    
switch experimentI
    case 1; experimentID = '2013-08-19-6';
    case 2; experimentID = '2012-08-09-3';
    case 3; experimentID = '2013-10-10-0';
    case 4; experimentID = '2012-09-27-3';
end

switch stimulusTestI
    case 1; stimulusTest = 'WN';
    case 2; stimulusTest = 'NSEM';
end

switch cellTypeI
    case 1; cellType = 'On Parasol';
    case 2; cellType = 'Off Parasol';
end
%% Load stimulus movie using RemoteDataToolbox

[testmovie, xval_mosaic] =  loadDataRGCFigure2(experimentI,stimulusTestI,cellTypeI);

nFrames = 1200; % Length of WN movie is 1200, take nFrames to limit natural movie to same length
testmovieshort = testmovie.matrix(:,:,1:nFrames); 

%% Show test movie
showFrames = 50;
ieMovie(testmovieshort(:,:,1:showFrames));

%% Generate outer segment object

% In this case, the RGC GLM calculation converts from the frame buffer
% values in the movie to the outer segment responses.  That form of the
% outer segment object is called 'displayRGB'.

os1 = osCreate('displayRGB'); 
os1 = osSet(os1, 'timeStep', 1/120);

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
params.experimentID = experimentID; 
params.cellType = cellType;
params.stimulusTest = stimulusTest; % WN or NSEM, from above

% Create object
innerRetina = irPhys(os1, params);
nTrials = 57; innerRetina = irSet(innerRetina,'numberTrials',nTrials);

%% Plot a few simple properties
switch cellTypeI
    case 1; cellInd = 2;
    case 2; cellInd = 31;
end
irPlotFig2Linear(innerRetina,cellInd);
%% Compute the inner retina response

% Lienar convolution
innerRetina = irCompute(innerRetina, os1);

% % % Spike computation
for tr = 1:ntrials
    innerRetina = irComputeSpikes(innerRetina, os1);
end

innerRetinaPSTH = mosaicGet(innerRetina.mosaic{1},'responsePsth');

%% Create a new inner retina object and attach the recorded spikes
innerRetinaRecorded = irPhys(os1, params);
innerRetinaRecorded = irSet(innerRetinaRecorded,'numberTrials',nTrials);

innerRetinaRecorded = irSet(innerRetinaRecorded,'recordedSpikes',xval_mosaic);
innerRetinaRecordedPSTH = mosaicGet(innerRetinaRecorded.mosaic{1},'responsePsth');

%% Compare Rasters and PSTHs for a particular cell
switch cellTypeI
    case 1; cellInd = 2;
    case 2; cellInd = 31;
end
irPlotFig2Raster(innerRetina, innerRetinaRecorded,cellInd,stimulusTestI);
irPlotFig2PSTH(innerRetina, innerRetinaPSTH, innerRetinaRecordedPSTH,cellInd,stimulusTestI);

%% Calculate fractional variance predicted

fractionalVariance{experimentI,stimulusTestI,cellTypeI} = ...
    calculateFractionalVariance(innerRetinaPSTH, innerRetinaRecordedPSTH, stimulusTestI);

if (stimulusTestI == 2) && (plotFracFlag == 1); 
    irPlotFig2FracVar(experimentI,cellTypeI,fractionalVariance); 
end;

plotFracFlag = 1;
%%%
end%stimulusTestI
end%cellTypeI
end%experimentI
