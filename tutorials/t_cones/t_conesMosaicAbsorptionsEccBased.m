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
stimFOV = 0.4;
mosaicFOV = 0.6;
parms.freq = 5;
parms.contrast = 1.0;
parms.ph = 0;
parms.ang = 0;
parms.row = 128;
parms.GaborFlag = .12;
parms.col = 128;
[scene, parms] = sceneCreate('harmonic', parms);
scene = sceneSet(scene, 'fov', stimFOV);

% Then the oi
oi = oiCreate;
oi = oiCompute(oi,scene);

%% Build a default cone mosaic and compute isomerizatoins
loadPreviouslySavedMosaic = true;
mosaicFileName = sprintf('mosaic.mat');
if (~loadPreviouslySavedMosaic)
    cMosaic  = coneMosaicHex(9, ...
            'name', 'test', ...
            'fovDegs', mosaicFOV, ...
            'eccBasedConeDensity', true, ...
            'eccBasedConeQuantalEfficiency', ~true, ...
            'sConeMinDistanceFactor', 2.0, ... 
            'sConeFreeRadiusMicrons', 0.15*300, ...                   
            'spatialDensity', [0 6/10 3/10 1/10], ...
            'latticeAdjustmentPositionalToleranceF', 0.5, ...         
            'latticeAdjustmentDelaunayToleranceF', 0.05, ... 
            'maxGridAdjustmentIterations', 200, ...
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

% No correction for ecc-dependent change in cone efficiency
cMosaic.eccBasedConeQuantalEfficiency = ~true;

cMosaic2 = cMosaic.copy();
cMosaic2.eccBasedConeQuantalEfficiency = true;

%% Generate two em sequences each with  40 eye posistions.
eyeMovementsNum = 40;
cMosaic.emGenSequence(eyeMovementsNum, 'nTrials', 2);
theEmPositions = cMosaic.emPositions;

subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 3, ...
       'heightMargin',  0.03, ...
       'widthMargin',    0.03, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.03, ...
       'topMargin',      0.02);

% Generate figure
hFig = figure(1); clf;
set(hFig, 'Position', [10 10 1295 792], 'Color', [1 1 1]); 

generateVideo = true;
if (generateVideo)
   writerObj = VideoWriter('video.m4v', 'MPEG-4'); % H264 format
   writerObj.FrameRate = 30; 
   writerObj.Quality = 100;
   % Open video stream
   open(writerObj);
end

% Plot the scene
subplot('Position', subplotPosVectors(2,3).v);
xspace = (1:sceneGet(scene, 'cols'));
xspace = xspace - mean(xspace);
xspace = xspace/max(xspace) * stimFOV/2;
yspace = (1:sceneGet(scene, 'rows'));
yspace = yspace - mean(yspace);
yspace = yspace/max(yspace) * stimFOV/2;
luminanceMap = sceneGet(scene, 'luminance');

meanLum = sceneGet(scene, 'mean luminance');
imagesc(xspace, yspace, 0.5 + 0.5*(luminanceMap - meanLum)/meanLum);
set(gca, 'CLim', [0 1]);
axis('image');
hold on;
plot([xspace(1) xspace(end)], [0 0], 'g-', 'LineWidth', 1.5);
plot([0 0],[yspace(1) yspace(end)], 'g-', 'LineWidth', 1.5);
hold off;        

xlabel('space (degs)');
lut = gray(1024);
colormap(gca,lut);
set(gca, 'XTick', stimFOV/2*[-1 0 1], 'YTick', stimFOV/2*[-1 0 1],  ...
    'XLim', stimFOV/2*[-1 1],  'YLim', stimFOV/2*[-1 1], 'YTickLabels', {}, ...
    'FontSize', 12, 'Color', squeeze(lut(512,:)));
hC = colorbar();
hC.TickLabels = {};
title('scene luminance');
    
% Plot the two mosaics
ax = subplot('Position', subplotPosVectors(1,1).v);
cMosaic.visualizeGrid(...
    'axesHandle', ax, ...
    'visualizedConeAperture', 'lightCollectingArea', ...
    'apertureShape', 'disks', ...
    'overlayHexMesh', true);
set(ax, 'XLim', stimFOV/2*[-1 1]*cMosaic2.micronsPerDegree*1e-6, ...
     'YLim', stimFOV/2*[-1 1]*cMosaic2.micronsPerDegree*1e-6, ...
     'XTick', [], 'YTick', [], 'FontSize', 12);
title(ax,'mosaic with constant cone efficiency');
colormap(ax, gray(1024));
hC = colorbar(ax);
hC.Ticks = [];

ax2 = subplot('Position', subplotPosVectors(2,1).v);
cMosaic2.visualizeGrid(...
    'axesHandle', ax2, ...
    'visualizedConeAperture', 'lightCollectingArea', ...
    'apertureShape', 'disks', ...
    'overlayHexMesh', true);
set(ax2, 'XLim', stimFOV/2*[-1 1]*cMosaic2.micronsPerDegree*1e-6, ...
     'YLim', stimFOV/2*[-1 1]*cMosaic2.micronsPerDegree*1e-6, ...
     'XTick', [], 'YTick', [], 'FontSize', 12);
colormap(ax2, gray(1024));
hC = colorbar(ax2);
hC.Ticks = [];
title(ax2,'mosaic with ecc-based cone efficiency');
    
cMosaic.emPositions = theEmPositions;
cMosaic2.emPositions = theEmPositions;
    
%% Compute isomerizations
absorptions1 = cMosaic.compute(oi);
absorptions2 = cMosaic2.compute(oi);
    
[~,idx] = max(absorptions1(:));
[~,maxConeRow,maxConeCol,timeBinOfMaxActivation] = ind2sub(size(absorptions2), idx);
    
signalRange = [min(absorptions2(:)) max(absorptions2(:))];

% Plot the time course of the responses
for visualizedTrialIndex = 1:2

    % Plot mean activation during the course of the emPath
    for k = 1:eyeMovementsNum
        
        % Extract the activation patterns for current trial and time bin
        visualizedActivationPattern1 = squeeze(absorptions1(visualizedTrialIndex,:,:, k));
        visualizedActivationPattern2 = squeeze(absorptions2(visualizedTrialIndex,:,:, k));
        residual = visualizedActivationPattern2-visualizedActivationPattern1;
        
        % Plot activation of first mosaic
        subplot('Position', subplotPosVectors(1,2).v);
        cMosaic.renderActivationMap(gca, visualizedActivationPattern1, ...
             'signalRange', signalRange, ...
             'visualizedConeAperture', 'geometricArea', ...
             'mapType', 'modulated disks', ...
             'showColorBar', true, ...
             'labelColorBarTicks', true, ...
             'titleForColorBar', 'R*/cone/5msec', ...
             'colorMap', brewermap(1024, '*greys'), ...
             'crossHairPosition', squeeze(theEmPositions(visualizedTrialIndex,k,:))*cMosaic2.patternSampleSize(1), ...
             'visualizedFOV', stimFOV, ...
             'backgroundColor', [0 0 0]);
        xlabel(''); ylabel('');  set(gca, 'XTickLabels', {}, 'YTickLabels', {});
        title('absorptions in mosaic with constant cone efficiency');

        % Plot activation of second mosaic
        subplot('Position', subplotPosVectors(2,2).v);
        cMosaic2.renderActivationMap(gca, visualizedActivationPattern2, ...
             'signalRange', signalRange, ...
             'visualizedConeAperture', 'geometricArea', ...
             'mapType', 'modulated disks', ...
             'showColorBar', true, ...
             'labelColorBarTicks', true, ...
             'titleForColorBar', 'R*/cone/5msec', ...
             'colorMap', brewermap(1024, '*greys'), ...
             'crossHairPosition', squeeze(theEmPositions(visualizedTrialIndex,k,:))*cMosaic2.patternSampleSize(1), ...
             'visualizedFOV', stimFOV, ...
             'backgroundColor', [0 0 0]);
        ylabel(''); set(gca, 'YTickLabels', {});
        title('absorptions in mosaic with ecc-based cone efficiency');

        % Plot the difference in absorptions b/n the 2 mosaics
        subplot('Position', subplotPosVectors(1,3).v);

        cMosaic2.renderActivationMap(gca, residual, ...
             'signalRange', 50*[-1 1], ...
             'visualizedConeAperture', 'geometricArea', ...
             'mapType', 'modulated disks', ...
             'showColorBar', true, ...
             'labelColorBarTicks', true, ...
             'titleForColorBar', 'R*/cone/5msec', ...
             'crossHairPosition', squeeze(theEmPositions(visualizedTrialIndex,k,:))*cMosaic2.patternSampleSize(1), ...
             'visualizedFOV', stimFOV, ...
             'colorMap', brewermap(1024, '*RdBu'), ...
             'backgroundColor', [1 1 1]);
        ylabel(''); xlabel('');  set(gca, 'YTickLabels', {});
        title('absorptions difference');
    
        drawnow;
        if (generateVideo)
            frame = getframe(gcf);
            writeVideo(writerObj, frame);
        end
    end % for k
end % visualizedTrialIndex

if (generateVideo)
    close(writerObj);
end
    