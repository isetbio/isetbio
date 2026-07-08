%% Responses of filters in quadrature
%
% A quadrature pair contains two filters with the same spatial frequency
% and a 90 degree phase difference. The individual linear responses depend
% on stimulus phase. Their vector magnitude (amplitude) and squared
% magnitude (energy) do not.
%
% See also
%   harmonicP, imageHarmonic

%%
ieInit;

%% Create a one-dimensional stimulus and quadrature filter pair

nSamples = 128;
x = (0:nSamples-1)/nSamples;
frequencyCyclesPerImage = 2;

cosine = cos(2*pi*frequencyCyclesPerImage*x);
sine = sin(2*pi*frequencyCyclesPerImage*x);
stimulus = 0.4*square(2*pi*frequencyCyclesPerImage*x) + 0.5;

% Unit-length filters make their linear responses directly comparable.
cosineFilter = cosine/norm(cosine);
sineFilter = sine/norm(sine);

%% Compute responses as the stimulus phase changes

shifts = 0:nSamples-1;
cosineResponse = zeros(size(shifts));
sineResponse = zeros(size(shifts));

for ii = 1:numel(shifts)
    shiftedStimulus = circshift(stimulus, shifts(ii));
    cosineResponse(ii) = dot(shiftedStimulus, cosineFilter);
    sineResponse(ii) = dot(shiftedStimulus, sineFilter);
end

quadratureAmplitude = hypot(cosineResponse, sineResponse);
quadratureEnergy = quadratureAmplitude.^2;

fprintf('1D quadrature energy relative range: %.3g\n', ...
    range(quadratureEnergy)/mean(quadratureEnergy));

%% Plot the stimulus, linear responses, and combined quadrature response

hFig = ieFigure([], 'wide');
layout = tiledlayout(hFig, 1, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile(layout);
plot(ax, x, stimulus, 'k-', x, cosine, 'r-', x, sine, 'b-');
grid(ax, 'on');
xlabel(ax, 'Normalized position (image width)');
ylabel(ax, 'Value');
legend(ax, {'Stimulus', 'Cosine', 'Sine'});
title(ax, 'Stimulus and quadrature pair');

ax = nexttile(layout);
plot(ax, shifts, cosineResponse, 'r-', ...
    shifts, sineResponse, 'b-', ...
    shifts, quadratureAmplitude, 'k-', ...
    'LineWidth', 1.5);
grid(ax, 'on');
xlabel(ax, 'Circular shift (samples)');
ylabel(ax, 'Response');
legend(ax, {'Cosine response', 'Sine response', 'Combined amplitude'});
title(ax, 'Linear responses rotate with phase');

ax = nexttile(layout);
plot(ax, shifts, quadratureEnergy/quadratureEnergy(1), 'k-', ...
    'LineWidth', 1.5);
grid(ax, 'on');
xlabel(ax, 'Circular shift (samples)');
ylabel(ax, 'Relative energy');
ylim(ax, [0.95 1.05]);
title(ax, 'Combined energy is phase invariant');

%% Repeat the quadrature calculation for a two-dimensional image

imageRows = 128;
stimulusImage = repmat(stimulus, [imageRows 1]);
cosineImage = repmat(cosine, [imageRows 1]);
sineImage = repmat(sine, [imageRows 1]);

cosineImageFilter = cosineImage/norm(cosineImage(:));
sineImageFilter = sineImage/norm(sineImage(:));

cosineImageResponse = zeros(size(shifts));
sineImageResponse = zeros(size(shifts));
for ii = 1:numel(shifts)
    shiftedImage = circshift(stimulusImage, shifts(ii), 2);
    cosineImageResponse(ii) = dot(shiftedImage(:), cosineImageFilter(:));
    sineImageResponse(ii) = dot(shiftedImage(:), sineImageFilter(:));
end

imageQuadratureEnergy = cosineImageResponse.^2 + sineImageResponse.^2;
fprintf('2D quadrature energy relative range: %.3g\n', ...
    range(imageQuadratureEnergy)/mean(imageQuadratureEnergy));

%% Localize the quadrature filters with a Gaussian envelope

gaussianSpread = 32;
gaussianEnvelope = fspecial('gaussian', ...
    size(stimulusImage), gaussianSpread);

windowedCosineFilter = gaussianEnvelope .* cosineImage;
windowedSineFilter = gaussianEnvelope .* sineImage;
windowedCosineFilter = windowedCosineFilter/norm(windowedCosineFilter(:));
windowedSineFilter = windowedSineFilter/norm(windowedSineFilter(:));

windowedCosineResponse = zeros(size(shifts));
windowedSineResponse = zeros(size(shifts));
for ii = 1:numel(shifts)
    shiftedImage = circshift(stimulusImage, shifts(ii), 2);
    windowedCosineResponse(ii) = dot(shiftedImage(:), ...
        windowedCosineFilter(:));
    windowedSineResponse(ii) = dot(shiftedImage(:), ...
        windowedSineFilter(:));
end

windowedQuadratureEnergy = ...
    windowedCosineResponse.^2 + windowedSineResponse.^2;
fprintf('Gaussian-windowed quadrature energy relative range: %.3g\n', ...
    range(windowedQuadratureEnergy)/mean(windowedQuadratureEnergy));

%% Plot the two-dimensional and localized quadrature responses

hFig = ieFigure([], 'wide');
layout = tiledlayout(hFig, 1, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile(layout);
imagesc(ax, stimulusImage);
axis(ax, 'image');
colormap(ax, gray);
title(ax, 'Two-dimensional stimulus');

ax = nexttile(layout);
plot(ax, shifts, imageQuadratureEnergy/imageQuadratureEnergy(1), ...
    'k-', 'LineWidth', 1.5);
grid(ax, 'on');
xlabel(ax, 'Horizontal shift (samples)');
ylabel(ax, 'Relative energy');
ylim(ax, [0.95 1.05]);
title(ax, 'Unwindowed quadrature energy');

ax = nexttile(layout);
plot(ax, shifts, windowedQuadratureEnergy/windowedQuadratureEnergy(1), ...
    'k-', 'LineWidth', 1.5);
grid(ax, 'on');
xlabel(ax, 'Horizontal shift (samples)');
ylabel(ax, 'Relative energy');
title(ax, 'Gaussian-windowed quadrature energy');

%% Compute quadrature energy from cone excitations
%
% Repeat the calculation after the harmonic images pass through human
% optics and are sampled by a small cMosaic. The cone-excitation patterns
% produced by harmonics at 0 and 90 degrees serve as the quadrature filters.
%
% We use only L cones for this test. In a natural mixed L/M/S mosaic,
% neighboring cones have different spectral sensitivities and therefore
% different responses to the same image intensity. Those cone-class
% differences can confound this simple spatial quadrature calculation.
% Restricting the calculation to one cone class removes that source of
% error, although it does not guarantee lower error for every particular
% lattice and spatial frequency. Irregular cone positions and optics still
% prevent perfect phase invariance.

coneMosaicSizeDegs = [0.05 0.05];
coneSceneFOVDegs = 0.08;
coneHarmonicSpatialFrequencyCPD = 50;
coneStimulusPhases = linspace(0, 2*pi, 17);
coneStimulusPhases(end) = [];

cm = cMosaic( ...
    'sizeDegs', coneMosaicSizeDegs, ...
    'eccentricityDegs', [0 0], ...
    'coneDensities', [1 0 0 0], ...
    'integrationTime', 10/1000, ...
    'noiseFlag', 'none', ...
    'randomSeed', 1);

coneExcitations = zeros(cm.conesNum, numel(coneStimulusPhases));
for ii = 1:numel(coneStimulusPhases)
    hparams = harmonicP( ...
        'freq', coneHarmonicSpatialFrequencyCPD*coneSceneFOVDegs, ...
        'contrast', 1, ...
        'ph', coneStimulusPhases(ii), ...
        'ang', 0, ...
        'row', 256, ...
        'col', 256);

    scene = sceneCreate('harmonic', hparams);
    scene = sceneSet(scene, 'fov', coneSceneFOVDegs);

    oi = oiCreate('human');
    oi = oiCompute(oi, scene, 'pad value', 'mean');
    coneExcitations(:, ii) = squeeze(cm.compute(oi));
end

% Remove each cone's mean response across phase before using the excitation
% patterns at 0 and 90 degrees as the quadrature filters.
centeredConeExcitations = coneExcitations - mean(coneExcitations, 2);
quarterCycleIndex = find(coneStimulusPhases == pi/2, 1);

coneCosineFilter = centeredConeExcitations(:, 1);
coneSineFilter = centeredConeExcitations(:, quarterCycleIndex);
coneCosineFilter = coneCosineFilter/norm(coneCosineFilter);
coneSineFilter = coneSineFilter/norm(coneSineFilter);

coneCosineResponse = coneCosineFilter' * centeredConeExcitations;
coneSineResponse = coneSineFilter' * centeredConeExcitations;
coneQuadratureEnergy = coneCosineResponse.^2 + coneSineResponse.^2;

coneFilterCorrelation = dot(coneCosineFilter, coneSineFilter);
coneFilterNormRatio = norm(centeredConeExcitations(:, 1)) / ...
    norm(centeredConeExcitations(:, quarterCycleIndex));
coneEnergyRelativeRange = range(coneQuadratureEnergy) / ...
    mean(coneQuadratureEnergy);
fprintf('Cone quadrature-filter correlation: %.3g\n', ...
    coneFilterCorrelation);
fprintf('Cone quadrature-filter norm ratio: %.3g\n', ...
    coneFilterNormRatio);
fprintf('Cone quadrature energy relative range: %.3g\n', ...
    coneEnergyRelativeRange);

%% Plot the cone mosaic quadrature calculation

cm.plot('mosaic', [], ...
    'plot title', sprintf('Quadrature sampling by a %.2g degree all-L cMosaic', ...
    coneMosaicSizeDegs(1)));

hFig = ieFigure([], 'wide');
layout = tiledlayout(hFig, 1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile(layout);
plot(ax, rad2deg(coneStimulusPhases), coneCosineResponse, 'r-o', ...
    rad2deg(coneStimulusPhases), coneSineResponse, 'b-o', ...
    'LineWidth', 1.5);
grid(ax, 'on');
xlabel(ax, 'Stimulus phase (degrees)');
ylabel(ax, 'Response');
legend(ax, {'0 degree cone filter', '90 degree cone filter'});
title(ax, sprintf('Sampled filter correlation: %.3f', ...
    coneFilterCorrelation));

ax = nexttile(layout);
plot(ax, rad2deg(coneStimulusPhases), ...
    coneQuadratureEnergy/mean(coneQuadratureEnergy), ...
    'k-o', 'LineWidth', 1.5);
grid(ax, 'on');
xlabel(ax, 'Stimulus phase (degrees)');
ylabel(ax, 'Energy relative to mean');
title(ax, sprintf('Relative energy range: %.3f', ...
    coneEnergyRelativeRange));
set(ax, 'YLim', [0 2]);

%%
