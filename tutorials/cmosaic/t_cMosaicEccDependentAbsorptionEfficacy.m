% Demostrate different ecc-dependent capabilties of the new @cMosaic object
%
% Description:
%    Shows how to set different ecc-dependent properties of the @cMosaic,
%    and demonstrates their effect on the computed mean response. 
%    The following ecc-dependent properties are examined:
%    - eccVaryingConeAperture
%    - eccVaryingOuterSegmentLength
%    - eccVaryingConeBlur
%    - eccVaryingMacularPigmentDensity
%    - eccVaryingMacularPigmentDensityDynamic
%
% See Also:
%   t_cMosaicBasic
%   t_cMosaicBenchMark
%   t_cMosaicOffAxis
%   t_cMosaicFromConeMosaicHex

% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;
clear;
close all;

%% Generate a grating scene (3 deg FOV, 20 c/deg, with a res of 2048x2048 pixels)
pixelsNum = 2048;
fovDegs = 3.0;
stimFreqCyclesPerDeg = 20;

parms.freq = stimFreqCyclesPerDeg*fovDegs;
parms.contrast = 1;
parms.ph = 0;
parms.ang = 0;
parms.row = pixelsNum;
parms.col = pixelsNum;
parms.GaborFlag = 0;
scene = sceneCreate('harmonic', parms);
scene = sceneSet(scene, 'fov', fovDegs);

%% Compute the optical image
oi = oiCreate('human');
oi = oiCompute(oi,scene,'pad value','mean');

%% Set mosaic size and eccentricity
mosaicSize = [1.0 1.0];
mosaicEcc = [1 0];
figNo = 0;

%% Demonstrate the effect of ecc-dependent IS aperture and OS length
 cond1Struct = struct(...
    'name', 'ecc-dependent IS area / OS length', ...
    'eccVaryingConeAperture', true, ...
    'eccVaryingOuterSegmentLength', true);
cond2Struct = struct(...
    'name', 'median value of IS area / OS length', ...
    'eccVaryingConeAperture', ~true, ...
    'eccVaryingOuterSegmentLength', ~true);
figNo = figNo + 1;
runIt(oi, mosaicSize, mosaicEcc, cond1Struct, cond2Struct, figNo, 'IS aperture/OS length');


%% Demonstrate the effect of ecc-dependent IS aperture blur
% Note. The effect of ecc-dependent aperture blur will be captured more 
% precisely as the #of pixels in the scene (resolution) are increased. 
% If the pixel size is greater than the cone aperture, then this effect
% cannot be captured. You can ask the mosaic for a suggested pixel size (in
% degrees as follows):
%   eccVaryingConeBlur = ~true;
%   pixelSizeDegs = cMosaic.suggestedScenePixelSizeDegs(eccVaryingConeBlur);
cond1Struct = struct(...
    'name', 'ecc-dependent aperture blur', ...
    'eccVaryingConeBlur', true, ...
    'eccVaryingConeAperture', true, ...
    'eccVaryingOuterSegmentLength', true);
cond2Struct = struct(...
    'name', 'median value of aperture blur', ...
    'eccVaryingConeBlur', ~true, ...
    'eccVaryingConeAperture', true, ...
    'eccVaryingOuterSegmentLength', true);
figNo = figNo + 1;
runIt(oi, mosaicSize, mosaicEcc, cond1Struct, cond2Struct, figNo, 'aperture blur');


%% Demonstrate the effect of aperture shape: pillbox vs Gaussian
c1 = struct('smoothLocalVariations', true, 'shape','Pillbox');
c2 = struct('smoothLocalVariations', true, 'shape','Gaussian', 'sigma', 0.204);
cond1Struct = struct(...
    'name', 'pillbox aperture', ...
    'eccVaryingConeBlur', true, ...
    'eccVaryingConeAperture', true, ...
    'eccVaryingOuterSegmentLength', true, ...
    'coneApertureModifiers', c1);
cond2Struct = struct(...
    'name', 'Gaussian aperture', ...
    'eccVaryingConeBlur', true, ...
    'eccVaryingConeAperture', true, ...
    'eccVaryingOuterSegmentLength', true, ...
    'coneApertureModifiers', c2);
figNo = figNo + 1;
runIt(oi, mosaicSize, mosaicEcc, cond1Struct, cond2Struct, figNo, 'aperture shape: pillbox (default) vs. Gaussian');


%% Demonstrate the effect of ecc-dependent macular pigment density (no eye movements)
cond1Struct = struct(...
    'name', 'ecc-dependent MP density', ...
    'eccVaryingMacularPigmentDensity', true, ...
    'eccVaryingConeBlur', ~true, ...
    'eccVaryingConeAperture', true, ...
    'eccVaryingOuterSegmentLength', true);
cond2Struct = struct(...
    'name', 'median value of MP density', ...
    'eccVaryingMacularPigmentDensity', ~true, ...
    'eccVaryingConeBlur', ~true, ...
    'eccVaryingConeAperture', true, ...
    'eccVaryingOuterSegmentLength', true);
figNo = figNo + 1;
runIt(oi, mosaicSize, mosaicEcc, cond1Struct, cond2Struct, figNo, 'macular pigment');


%% Demonstrate the effect of correcting the ecc-dependent macular pigment (with eye movements)
cond1Struct = struct(...
    'name', 'ecc-dependent MP (correction for eye movements)', ...
    'eccVaryingMacularPigmentDensityDynamic', true, ...
    'eccVaryingMacularPigmentDensity', true, ...
    'eccVaryingConeBlur', ~true, ...
    'eccVaryingConeAperture', true, ...
    'eccVaryingOuterSegmentLength', true);
cond2Struct = struct(...
    'name', 'ecc-dependent MP (no correction for eye movements)', ...
    'eccVaryingMacularPigmentDensityDynamic', ~true, ...
    'eccVaryingMacularPigmentDensity', true, ...
    'eccVaryingConeBlur', ~true, ...
    'eccVaryingConeAperture', true, ...
    'eccVaryingOuterSegmentLength', true);

% Generate eye movement data. These are store
figNo = 1000;
figNo = figNo + 1;
mosaicEcc = [0 0];
runItWithEyeMovements(oi, mosaicSize, mosaicEcc, cond1Struct, cond2Struct, figNo, 'macular pigment dynamic');

%% END SCRIPT


function runItWithEyeMovements(oi, mosaicSize, mosaicEcc, cond1Struct, cond2Struct, figNo, plotTitle)
    % Generate an ON-axis mosaic (ecc: 0,0)
    cm = cMosaic(...
        'sizeDegs', mosaicSize, ...      
        'eccentricityDegs', mosaicEcc ... 
        );
    
    % 100 msec eye movement data
    emDurationSeconds = 100/1000;
    % Alter micro-saccade interval to shorten the simulation
    microSaccadeIntervalSeconds = 50/1000;
    % Generate eye movements
    cm.emGenSequence(emDurationSeconds, ...
        'microsaccadeType', 'heatmap/fixation based', ...
        'microsaccadeMeanIntervalSeconds', microSaccadeIntervalSeconds, ...
        'nTrials', 1, ...
        'randomSeed', 10);
    
    % Compute response for condition1 
    fNames = fieldnames(cond1Struct);
    for k = 1:numel(fNames)
        cm.(fNames{k}) = cond1Struct.(fNames{k});
    end
    r1 = cm.compute(oi, 'withFixationalEyeMovements', true);
    
    % Compute response for condition2 
    fNames = fieldnames(cond2Struct);
    for k = 1:numel(fNames)
        cm.(fNames{k}) = cond2Struct.(fNames{k});
    end
    r2 = cm.compute(oi, 'withFixationalEyeMovements', true);
 
    % Visualize response time series
    [instancesNum, timeSeriesLength, conesNum] = size(r1);
    
    for timePoint = 1:timeSeriesLength  
        figFrame(timePoint) = visualizeResponses(cm, squeeze(r1(1,timePoint,:)), squeeze(r2(1,timePoint,:)), ...
            struct('trial', 1, 'timePoints', 1:timePoint) , ...  % Show eye movement data up to this time point
            cond1Struct.name, cond2Struct.name, ...
            figNo+timePoint, sprintf('%s (time sample: %d)', plotTitle, timePoint));
    end
    
    % Animate the sequence of frames
    for k = 1:10
        for iFig = 1:numel(figFrame)
            figure(figFrame(iFig));
            pause(0.1);
        end
    end
end

%%
function runIt(oi, mosaicSize, mosaicEcc, cond1Struct, cond2Struct, figNo, plotTitle)
    cm = cMosaic(...
        'sizeDegs', mosaicSize, ...      
        'eccentricityDegs', mosaicEcc ... 
        );

    % Compute response for condition1 
    fNames = fieldnames(cond1Struct);
    for k = 1:numel(fNames)
        cm.(fNames{k}) = cond1Struct.(fNames{k});
    end
    r1 = cm.compute(oi);
    
    % Compute response for condition2 
    fNames = fieldnames(cond2Struct);
    for k = 1:numel(fNames)
        cm.(fNames{k}) = cond2Struct.(fNames{k});
    end
    r2 = cm.compute(oi);

    % Visualize responses
    displayedEyeMovementStruct = [];
    visualizeResponses(cm, r1, r2, displayedEyeMovementStruct, ...
        cond1Struct.name, cond2Struct.name, ...
        figNo, plotTitle);

end

%%
function hFig = visualizeResponses(cm, r1, r2, emStruct, t1, t2, figNo, figName)
   sv = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 3, ...
       'heightMargin',  0.01, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.06, ...
       'topMargin',      0.03);
   
   activationRange(1) = min([min(r1(:)) min(r2(:))]);
   activationRange(2) = max([max(r1(:)) max(r2(:))]);
   cMap = brewermap(1024, 'reds');
   
   hFig = figure(figNo); clf;
   set(hFig, 'Position', [10 10 1400 550], 'Name', figName);
   ax = subplot('Position', sv(1,1).v);
   cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
                 'activation', r1(:), ...
                 'activationRange', activationRange, ...
                 'activationColorMap', cMap, ...
                 'displayedEyeMovementData', emStruct, ...
                 'horizontalActivationColorBar', true, ...
                 'colorBarTickLabelPostFix', 'R*', ...
                 'crossHairsOnOpticalImageCenter', true, ...
                 'plotTitle', sprintf('r1 (%s)', t1));
             
   ax = subplot('Position', sv(1,2).v);
   cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
                 'activation', r2(:), ...
                 'activationRange', activationRange, ...
                 'activationColorMap', cMap, ...
                 'displayedEyeMovementData', emStruct, ...
                 'horizontalActivationColorBar', true, ...
                 'colorBarTickLabelPostFix', 'R*', ...
                 'crossHairsOnOpticalImageCenter', true, ...
                 'plotTitle', sprintf('r2 (%s)', t2));
             
   diff = 100*(r1-r2)./r1;
   cMap = brewermap(1024, '*RdBu');
   activationRange = max(abs(diff(:)))*[-1 1];
   ax = subplot('Position', sv(1,3).v);
   cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
                 'activation', diff(:), ...
                 'activationRange', activationRange, ...
                 'activationColorMap', cMap, ...
                 'displayedEyeMovementData', emStruct, ...
                 'horizontalActivationColorBar', true, ...
                 'colorBarTickLabelPostFix', '%', ...
                 'crossHairsOnOpticalImageCenter', true, ...
                 'backgroundColor', [0 0 0], ...
                 'plotTitle', '100*(r1-r2)/r1');
             
end

%% END FUNCTIONS