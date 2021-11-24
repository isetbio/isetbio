% s_rgcGridNtrial
%
% Description:
%    Calculating multiple trials using the s_rgcGrid calculation.
%
%    Illustrate the use of the multiple trials. There are two ways to run
%    multiple trials.
%
%    If the stimulus is an oiSequence, then you can set the eye movement
%    path variable to include multiple trials and these will be returned.
%
%    If the stimulus is a single oi, you need to execute a loop.
%
%    Here we show how to make an oiSequence for an oi. And we show you to
%    run in a loop.
%
% Notes:
%    * Programming TODO: We could push the loop on the cone mosaic into the
%      cone mosaic calculation. Right now we don't accept multiple trials
%      for the empath when there is a single oi.
%    * We could make the bipolar compute always return a cell array.
%
% See Also:
%   s_rgcGrid
%

% History:
%    XX/XX/17  BW   ISETBIO Team, 2017
%    06/11/19  JNM  Documentation pass

%%
ieInit

%% Scene, oi, grid, cone mosaic
imSize = 128;
lineSpacing = 48;
fov = 1.5; % deg
scene = sceneCreate('grid lines', imSize, lineSpacing);
scene = sceneSet(scene, 'fov', fov);
oi = oiCreate;  % Standard human optics
oi = oiCompute(oi, scene);

%% Run 2 trials with the same stimulus
cMosaic = coneMosaic;
cMosaic.setSizeToFOV(fov);

%% New method for multiple cone trials when not an oiSequence
nMovements = 25;
nTrials = 2;
emPaths = cMosaic.emGenSequence(nMovements, 'nTrials', nTrials);
[coneAbsorptions, coneCurrent] = cMosaic.compute(oi, 'emPath', emPaths, ...
    'currentFlag', true);

% vcNewGraphWin;
% ieMovie(squeeze(coneAbsorptions(1, :, :, :)));
% ieMovie(squeeze(coneAbsorptions(2, :, :, :)));
% ieMovie(squeeze(coneCurrent(2, :, :, :)));
% ieMovie(squeeze(coneCurrent(1, :, :, :)));

%% Loop (old) method for multiple cone trials when not an oiSequence
% Returns the same coneAbsorptions and coneCurrent
%
% nMovements = 25;
% nTrials = 2;
% emPaths = cMosaic.emGenSequence(nMovements, 'nTrials', nTrials);
% coneAbsorptions = zeros(nTrials, cMosaic.rows, cMosaic.cols, nMovements);
% coneCurrent = zeros(size(coneAbsorptions));
% for ii = 1:nTrials
%     cMosaic.compute(oi, 'emPath', emPaths, 'currentFlag', true);
%     coneAbsorptions(ii, :, :, :) = cMosaic.absorptions;
%     coneCurrent(ii, :, :, :) = cMosaic.current;
%     emPaths = cMosaic.emGenSequence(nMovements, 'nTrials', nTrials);
% end

% vcNewGraphWin;
% ieMovie(squeeze(coneAbsorptions(1, :, :, :)));
% ieMovie(squeeze(coneAbsorptions(2, :, :, :)));
% ieMovie(squeeze(coneCurrent(2, :, :, :)));
% ieMovie(squeeze(coneCurrent(1, :, :, :)));

%
%% Make the bipolar layer with just one mosaic
clear bpL bpMosaicParams bpTrials
bpL = bipolarLayer(cMosaic);

% Compute multiple trials based on the cone mosaic current
ii = 1;
bpMosaicParams.spread = 2;  % RF diameter w.r.t. input samples
bpMosaicParams.stride = 2;  % RF diameter w.r.t. input samples
bpL.mosaic{ii} = bipolarMosaic(cMosaic, 'on midget', bpMosaicParams);
bpTrials{1} = bpL.mosaic{ii}.compute('coneTrials', coneCurrent);
% bpL.window;

%% Make the RGC layer and show it
clear rgcL rgcParams
rgcL = rgcLayer(bpL);

% Spread and stride are not working
rgcParams.rfDiameter = 2;

% rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1}, 'on midget');
rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1}, 'on midget', rgcParams);
nTrialsSpikes = rgcL.compute('bipolarScale', 50, ...
    'bipolarContrast', 0.4, 'bipolarTrials', bpTrials, 'coupling', false);

rgcL.window;

%% Show the first trial
vcNewGraphWin([], 'tall');
subplot(2, 1, 1)
imagesc(sum(squeeze(nTrialsSpikes{1}(1, :, :, :)), 3));
subplot(2, 1, 2)
imagesc(sum(squeeze(nTrialsSpikes{1}(2, :, :, :)), 3));
