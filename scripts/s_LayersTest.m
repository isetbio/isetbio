%% s_Layers.m
%
% Testing the new layer architecture.
%
% BW ISETBIO Team, 2017

%%
ieInit;

%% Get the data from the RDT site

rdt = RdtClient('isetbio');
rdt.crp('/resources/data/cmosaics');
% rdt.listArtifacts('type','mat','print',true);

% coneMosaicDataFixedEye;
% vaConeMosaic
data = rdt.readArtifact('coneMosaicDataFixedEye', 'type', 'mat');
cMosaic   = data.cMosaic;
alignedC = data.alignedC;

% cMosaic.window;

%% Create a set of bipolar cell types in the bipolar mosaic

bpL = bipolarLayer(cMosaic);

% Make each type of bipolar mosaic
cellType = {'ondiffuse','offdiffuse','onmidget','offmidget','onSBC'};

clear bpMosaicParams
bpMosaicParams.rectifyType = 1;  % Experiment with this

bpMosaic  = cell(1,length(cellType));
bpNTrials = cell(1,length(cellType));
for ii = 1:length(cellType)
    
    bpMosaicParams.cellType = cellType{ii};
    
    bpMosaic{ii} = bipolarMosaic(cMosaic, bpMosaicParams);
    bpMosaic{ii}.set('sRFcenter',1);
    bpMosaic{ii}.set('sRFsurround',0);
    
    [~, bpNTrialsCenterTemp, bpNTrialsSurroundTemp] = ...
        bpMosaic{ii}.compute(cMosaic,'coneTrials',alignedC);
    bpNTrials{ii} = bpNTrialsCenterTemp - bpNTrialsSurroundTemp;
    
end
bpL.mosaic = bpMosaic;

% Try varying some experimental parameters
%
% bpL.mosaic{1}.set('sRFcenter',5); bpL.mosaic{1}.set('sRFsurround',1);
% Now compute and redisplay.
%
% disp('Computing bipolar responses'); [~, bpNTrialsCenter,
% bpNTrialsSurround] = bp.compute(cMosaic,'coneTrials',alignedC);
%
%
% TODO
%
% * After showing movie, the numbers on the mosaic axis are missing 
% * The units on the center size may not be correct. 
% * We need to allow changing the size of the center and surround on the bipolar.

%% Mosaics shown directly
% bpL.mosaic{1}.window;
% bpL.mosaic{2}.window;

%% The size of the RFs are surprising here

% The bipolar sizes we model are all the same for all the types.  The
% increaed RF size of the corresponding RGCs arise from spatial summation
% at the next synapse in the model.
bpL.window;

%% Retinal ganlion cell model

clear rgcLayer

% Create retina ganglion cell layer object
rgcL = rgcLayer(bpL);

% There are various parameters you could set.  We will write a script
% illustrating these later.  We need a description.
mosaicParams.centerNoise = 0;
mosaicParams.ellipseParams = [1 1 0];  % Principle, minor and theta
% mosaicParams.axisVariance = .1;
mosaicParams.type  = cellType;
mosaicParams.model = 'GLM';

diameters = [5 5 3 3 3];  % In microns.

cellType = {'on parasol','off parasol','on midget','off midget','smallbistratified'};
for ii = 1:length(cellType)
    mosaicParams.rfDiameter = diameters(ii);
    mosaicParams.type = cellType{ii};
    mosaicParams.inMosaic = 1;   % Could switch up and match inputs to outputs
    rgcL.mosaicCreate(mosaicParams);
end

nTrials = 1; rgcL.set('numberTrials',nTrials);

%% Compute the inner retina response and visualize

% Number of trials refers to number of repeats of the same stimulus
disp('Computing rgc responses');
[rgcL, nTrialsSpikes] = rgcL.compute(bpMosaic,'bipolarTrials',bpNTrials);


%% Show the layer

rgcL.window;

%% Retinal ganglion cell layer window

% TODO:
%   Put up 'ieInWindowMessage() when the movie is playing Label the
%   distances on the x and y axes
%

% I wish we could have two windows up at the same time.  Read the code to
% see why it is always the same window.
rgcL.mosaic{1}.window;
%%
rgcL.mosaic{3}.window;
%%
rgcL.mosaic{5}.window;  % This seems off to BW, slanted line???

%%