%%t_coneMosaicAbsorptionsEccBased  Ecc-based cone efficiency correction
%
% Description:
%    Demonstrates the effect of applying eccentricity-based cone efficiency
%    correction to a static stimulus. This demo includes eye movements.
%
% 
% NPC, ISETBIO Team, 2018
%
% 06/16/18  npc  Wrote it.

%% Initialize
ieInit; clear all;

%% Build a simple scene and oi (retinal image) for computing
%
% Make a gabor scene
fieldOfView = 0.5;
parms.freq = 4;
parms.contrast = 1;
parms.ph = 0;
parms.ang = 0;
parms.row = 128;
parms.col = 128;
parms.GaborFlag = 0;
[scene, parms] = sceneCreate('harmonic', parms);
scene = sceneSet(scene, 'fov', fieldOfView);

% Then the oi
oi = oiCreate;
oi = oiCompute(oi,scene);

%% Build a default cone mosaic and compute isomerizatoins
saveMosaic = ~true;
mosaicFileName = sprintf('mosaic.mat');
if (saveMosaic)
    cMosaic  = coneMosaicHex(7, ...
            'name', 'test', ...
            'fovDegs', fieldOfView, ...
            'eccBasedConeDensity', true, ...
            'eccBasedConeQuantalEfficiency', ~true, ...
            'sConeMinDistanceFactor', 2.0, ... 
            'sConeFreeRadiusMicrons', 0.15*300, ...                   
            'spatialDensity', [0 6/10 3/10 1/10], ...
            'latticeAdjustmentPositionalToleranceF', 0.5, ...         
            'latticeAdjustmentDelaunayToleranceF', 0.05, ... 
            'maxGridAdjustmentIterations', 20, ...
            'marginF', [] ... 
        );
    save(mosaicFileName, 'cMosaic', '-v7.3');
else
    load(mosaicFileName, 'cMosaic');
end

% Add photon noise
cMosaic.noiseFlag = 'none'; % 'random';

% Set the desired integration time, here 5 milliseconds
cMosaic.integrationTime = 5/1000;

cMosaic.eccBasedConeQuantalEfficiency = ~true;
cMosaic.visualizeGrid(...
        'generateNewFigure', true, ...
        'visualizedConeAperture', 'lightCollectingArea', ...
        'apertureShape', 'disks', ...
        'overlayHexMesh', true, ...
        'generateNewFigure', true);
    

%% Generate a sequence of 10 eye posistions.
eyeMovementsNum = 10;
cMosaic.emGenSequence(eyeMovementsNum);
cMosaic.emPositions = cMosaic.emPositions;
theEmPositions = cMosaic.emPositions;


%% Compute isomerizations
fprintf('Computing absorptions (WITHOUT ecc-based correction in absorption efficiency)\n');
absorptions1 = cMosaic.compute(oi);
visualizedTrialIndex = 1;
[~,idx] = max(absorptions1(:));
[~,maxConeRow,maxConeCol,timeBinOfMaxActivation] = ind2sub(size(absorptions1), idx);
visualizedActivationPattern1 = squeeze(absorptions1(visualizedTrialIndex ,:,:,timeBinOfMaxActivation));
signalRange1 = [min(visualizedActivationPattern1(:)) max(visualizedActivationPattern1(:))];


%% Compute isomerizations again, now with ecc-based corrections in cone efficiency
cMosaic2 = cMosaic.copy();
fprintf('Computing absorptions (WITH ecc-based correction in absorption efficiency)\n');
cMosaic2.eccBasedConeQuantalEfficiency = true;
cMosaic2.emPositions = theEmPositions;

cMosaic2.visualizeGrid(...
        'generateNewFigure', true, ...
        'visualizedConeAperture', 'lightCollectingArea', ...
        'apertureShape', 'disks', ...
        'overlayHexMesh', true, ...
        'generateNewFigure', true);
    
absorptions2 = cMosaic2.compute(oi);
visualizedActivationPattern2 = squeeze(absorptions2(visualizedTrialIndex ,:,:,timeBinOfMaxActivation));
signalRange2 = [min(visualizedActivationPattern2(:)) max(visualizedActivationPattern2(:))];

% Compute signal range
signalRange = [min([signalRange1(1) signalRange2(1)]) max([signalRange1(2) signalRange2(2)])];

% Generate figure
hFig = figure(100); clf;
set(hFig, 'Position', [10 10 1162 982], 'Color', [1 1 1]); 
subplot(2,2,3);
imagesc(sceneGet(scene, 'luminance'));
colormap(gca,gray(1024));
set(gca, 'XTick', [], 'YTick', [],  'FontSize', 12);
hC = colorbar();
hC.TickLabels = {}; 
axis('image');
title('scene luminance');


lut = brewermap(1024, '*greys');

generateVideo = false;
if (generateVideo)
   writerObj = VideoWriter('video.m4v', 'MPEG-4'); % H264 format
   writerObj.FrameRate = 30; 
   writerObj.Quality = 100;
   % Open video stream
   open(writerObj);
end
    

for k = 1:eyeMovementsNum

    subplot(2,2,1);
    visualizedActivationPattern1 = squeeze(absorptions1(visualizedTrialIndex ,:,:, k));
    cMosaic.renderActivationMap(gca, visualizedActivationPattern1, ...
             'signalRange', signalRange, ...
             'mapType', 'modulated disks', ...
             'colorMap',  lut);
    colorbar();
    set(gca, 'Color', [0 0 0],  'FontSize', 12)
    title('absorptions without ecc-based corrections');

    subplot(2,2,2);
    visualizedActivationPattern2 = squeeze(absorptions2(visualizedTrialIndex ,:,:, k));
    cMosaic2.renderActivationMap(gca, visualizedActivationPattern2, ...
             'signalRange', signalRange, ...
             'mapType', 'modulated disks', ...
             'colorMap',  lut);
    colorbar();
    set(gca, 'Color', [0 0 0],  'FontSize', 12);
    title('absorptions with ecc-based corrections');

    subplot(2,2,4);
    residual = visualizedActivationPattern2-visualizedActivationPattern1;
    cMosaic.renderActivationMap(gca, visualizedActivationPattern2-visualizedActivationPattern1, ...
             'signalRange', max(abs(residual(:)))*[-1 1], ...
             'mapType', 'modulated disks', ...
             'colorMap', brewermap(1024, '*RdBu'));
    hC = colorbar();
    hC.TickLabels = {}; 
    set(gca, 'Color', [1 1 1], 'FontSize', 12);
    title('differential absorptions');
    drawnow;
    if (generateVideo)
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end
end

if (generateVideo)
    close(writerObj);
end
    