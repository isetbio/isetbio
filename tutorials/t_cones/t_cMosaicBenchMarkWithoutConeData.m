% Benchmark the new @cMosaic object vs the old @coneMosaicHex
%
% Description:
%    Compares the total excitation of the old vs. the new cone mosaic
%    objects. Unlike the situation in t_cMosaicBenchMark, where the 
%    @cMosaic object is constructed from the coneData exported by the 
%    old coneMosaic object (so identical to it), here the @cMosaic object
%    is just matched to the old mosaic only in size and eccentricity. 
%    The 2 mosaics have different numbers of cones.
%
% See Also:
%   t_cMosaicBenchMark

% History:
%    09/29/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.
%% Initialize
ieInit;
clear;
close all;


%% Generate test scene
fprintf('Computing scenes ...\n');
sceneFOVDegs = 2.0;
% A sinusoid
params.freq = 10;
params.row = 1024;
params.col = 1024;
spectrum = 400:2:700;
scene1 = sceneCreate('Harmonic', params, spectrum);
%scene1 = sceneCreate('Harmonic', params);
scene1 = sceneSet(scene1, 'fov', sceneFOVDegs);

%% Compute optical images
fprintf('Computing optical images ...\n');
oi1 = oiCreate;
oi1 = oiCompute(scene1, oi1);


%% Generate the old-style cone mosaic
fprintf('Computing mosaics ...\n');
coneDensities = [0.6 0.3 0.1];
regenerateMosaic = ~true;
resourcesDir = fullfile(strrep(isetRootPath, 'isettools', 'tutorials'), 't_cones', 'resources');
if (regenerateMosaic)
    mosaicFOVDegs = sceneFOVDegs*0.6;
    theOldConeMosaic = coneMosaicHex(7, ...
        'spatialDensity', [0 coneDensities(1) coneDensities(2) coneDensities(3)], ...
        'fovDegs', mosaicFOVDegs, ...
        'eccBasedConeDensity', true, ...
        'eccBasedConeQuantalEfficiency', true, ...
        'eccBasedMacularPigment', true, ...
        'maxGridAdjustmentIterations', 10, ...
        'integrationTime', 100/1000, ...
        'noiseFlag', 'none' ...
            );
    save(fullfile(resourcesDir, 'coneMosaicHexBenchMark1.mat'), 'theOldConeMosaic');
else
    load(fullfile(resourcesDir, 'coneMosaicHexBenchMark1.mat'), 'theOldConeMosaic');
end
fprintf(' Done. \n'); 


theNewConeMosaic = cMosaic(...
    'micronsPerDegree', theOldConeMosaic.micronsPerDegree, ...
    'sizeDegs', theOldConeMosaic.fov, ...
    'coneDensities', coneDensities, ...
    'integrationTime', theOldConeMosaic.integrationTime, ...
    'noiseFlag', 'none');
        
[lConeIndices, mConeIndices,sConeIndices] = theOldConeMosaic.indicesForCones();
oldConeMosaicConesNum = numel(lConeIndices) + numel(mConeIndices) + numel(sConeIndices);
newConeMosaicConesNum = numel(theNewConeMosaic.lConeIndices) + numel(theNewConeMosaic.mConeIndices) + numel(theNewConeMosaic.sConeIndices);

%% Visualize the 2 mosaics
fprintf('Visualizing mosaics ...\n');
hFig = figure(1); clf;
set(hFig, 'Position', [10 10 1100 900]);
ax = subplot(2,2,1); cla(ax);
theOldConeMosaic.visualizeGrid(...
    'axesHandle', ax, ...
    'ticksInMicrons', true, ...
    'visualizedConeAperture',  'geometricArea', ...
    'plotTitle', sprintf('@coneMosaicHex (%d cones)', oldConeMosaicConesNum));

ax = subplot(2,2,2); cla(ax);
theNewConeMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'domain', 'microns', ...
    'visualizedConeAperture', 'geometricArea', ...
    'domainVisualizationTicks', struct('x', -200:50:200, 'y', -200:50:200), ...
    'plotTitle', sprintf('@cMosaic (%d cones)', newConeMosaicConesNum));


oldConeMosaicExcitations = theOldConeMosaic.compute(oi1);
newConeMosaicExcitations = theNewConeMosaic.compute(oi1);

activationRange(1) = min([min(oldConeMosaicExcitations(:)) min(newConeMosaicExcitations(:))]);
activationRange(2) = max([max(oldConeMosaicExcitations(:)) max(newConeMosaicExcitations(:))]);

% Compare the overall excitations ratio of the 2 mosaics
excitationsRatio = sum(sum(oldConeMosaicExcitations(find(theOldConeMosaic.pattern>1)))) / sum(newConeMosaicExcitations(:))

ax = subplot(2,2,3);
theOldConeMosaic.renderActivationMap(ax, oldConeMosaicExcitations, ...
            'signalRange', activationRange, ...
            'mapType', 'modulated disks', ...
            'titleForMap',  sprintf('mean excitation level: %f', mean(oldConeMosaicExcitations(find(theOldConeMosaic.pattern>1)))));

ax = subplot(2,2,4);
theNewConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'activation', newConeMosaicExcitations, ...
             'activationRange', activationRange, ...
             'visualizedConeAperture', 'geometricArea', ...
             'backgroundColor', [0 0 0], ...
             'plotTitle',  sprintf('mean excitation level: %f', mean(newConeMosaicExcitations(:))));

         