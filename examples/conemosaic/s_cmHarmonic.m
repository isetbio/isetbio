%% cMosaic responses to a spatial harmonic
%
% Compute the responses of a small cMosaic to a 10 cycles/degree harmonic
% during a short fixational eye-movement sequence. The final plots show the
% response time series of one S cone, two L cones, and two M cones.
%
% See also
%   cMosaic, harmonicP, sceneCreate, cMosaic.compute, cMosaic.emGenSequence

%%
ieInit;

%% Create the harmonic scene and optical image

sceneFOV = 0.6;
harmonicSpatialFrequencyCPD = 10;

% harmonicP specifies frequency in cycles/image. The 0.5 degree mosaic
% therefore spans about five cycles of this harmonic.
hparams = harmonicP( ...
    'freq', harmonicSpatialFrequencyCPD*sceneFOV, ...
    'contrast', 1, ...
    'ph', pi/2, ...
    'ang', 0, ...
    'row', 256, ...
    'col', 256);

scene = sceneCreate('harmonic', hparams);
scene = sceneSet(scene, 'fov', sceneFOV);

oi = oiCreate('human');
oi = oiCompute(oi, scene, 'pad value', 'mean');

%% Load a small cMosaic and create a short eye-movement sequence

cm = mosaicLoad([0.5 0.5], [0 0]);
cm.integrationTime = 10/1000;
cm.noiseFlag = 'none';

eyeMovementDurationSeconds = 200/1000;
cm.emGenSequence(eyeMovementDurationSeconds, ...
    'microsaccadeType', 'none', ...
    'nTrials', 1, ...
    'randomSeed', 1);

%% Compute the noise-free cone excitations

[excitations, ~, ~, ~, timeAxis] = cm.compute(oi, ...
    'withFixationalEyeMovements', true);

%% Select cones whose responses vary most over the eye-movement sequence

temporalRange = squeeze(max(excitations, [], 2) - min(excitations, [], 2));

[~, order] = sort(temporalRange(cm.sConeIndices), 'descend');
sConeIndex = cm.sConeIndices(order(1));

[~, order] = sort(temporalRange(cm.lConeIndices), 'descend');
lConeIndices = cm.lConeIndices(order(1:2));

[~, order] = sort(temporalRange(cm.mConeIndices), 'descend');
mConeIndices = cm.mConeIndices(order(1:2));

selectedConeIndices = [sConeIndex; lConeIndices(:); mConeIndices(:)];
selectedConeLabels = {'S cone', 'L cone', 'L cone', 'M cone', 'M cone'};
selectedConeColors = [0 0 1; 1 0 0; 1 0 0; 0 0.55 0; 0 0.55 0];
selectedResponses = squeeze(excitations(1, :, selectedConeIndices));

%% Show the mosaic, selected cones, and eye-movement path

[mosaicPlotData, ~] = cm.plot('eye movement path', [], ...
    'time point', 'last', ...
    'plot title', 'Mosaic, selected cones, and eye-movement path');

ax = mosaicPlotData.axesHandle;
selectedConePositions = cm.coneRFpositionsDegs(selectedConeIndices, :);
hold(ax, 'on');
for ii = 1:numel(selectedConeIndices)
    plot(ax, selectedConePositions(ii, 1), selectedConePositions(ii, 2), ...
        'o', ...
        'Color', selectedConeColors(ii, :), ...
        'MarkerSize', 10, ...
        'LineWidth', 2);
    text(ax, selectedConePositions(ii, 1), selectedConePositions(ii, 2), ...
        sprintf('  %s %d', selectedConeLabels{ii}, selectedConeIndices(ii)), ...
        'Color', selectedConeColors(ii, :), ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w');
end
hold(ax, 'off');

%% Plot the selected cone response time series

hFig = ieFigure([], 'wide');
layout = tiledlayout(hFig, 3, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for ii = 1:numel(selectedConeIndices)
    ax = nexttile(layout);
    plot(ax, timeAxis*1000, selectedResponses(:, ii), ...
        'Color', selectedConeColors(ii, :), ...
        'LineWidth', 1.5);
    grid(ax, 'on');
    xlabel(ax, 'Time (ms)');
    ylabel(ax, 'Cone excitations');
    title(ax, sprintf('%s %d', selectedConeLabels{ii}, ...
        selectedConeIndices(ii)));
end

title(layout, sprintf(['Responses to a %g cycles/degree harmonic ' ...
    '(approximately %g cycles across the mosaic)'], ...
    harmonicSpatialFrequencyCPD, ...
    harmonicSpatialFrequencyCPD*cm.sizeDegs(1)));

%%
