function t_wavefrontSampling
% Show how PSF and OTF are effected by different wavelength sampling params
%
% Syntax:
%   t_wavefrontSampling
%
% Description
%    Examine the effects of different sampling parameters in a wavefront on
%    the resulting PSF and the OTF.
%
%    (1) Shows that by increasing the 'spatial samples' param, the
%    frequency resolution with which the OTF is computed is increased. For
%    stimuli in which there is significant power in sp. freq. < 10 c/deg, 
%    the default setting of 201 samples is inadequate as the resulting OTF
%    is sampled very sparsely in this frequency range. In such cases, it is
%    recommended that the user sets the value of 'spatial samples' to some
%    value >> 201, say 501. See figure (1).
%     
%    (2) Shows that by increasing the 'reference pupil size', the spatial
%    resolution with which the PSF is computed is increased, and the
%    spatial support of the OTF is expanded. So increasing the reference
%    pupil size may be advantageous in cases where stimuli have a lot of
%    energy in the very high frequency range (perhaps in Adaptive Optics
%    experiments with small spots/high frequency gratings). See figure (2).
%
%    (3) Shows that when the 'sample interval domain' is set to 'pupil' the
%    PSF is sampled in manner that depends on the wavelength (sampling is
%    finer in shorter wavelengths than in longer wavelengths. This also
%    results in different spatial supports at different wavelengths. See
%    figure (3). Generally, we use 'sample interval domain' set to 'psf', 
%    because that keeps the spatial sampling consistent in the places we
%    generally care about.
%

% History
%    01/29/18  npc  Wrote it.
%    01/30/18  npc  Added 2Dmap plotView option
%    10/01/18  jnm  Formatting

%% Wavelengths to plot
wavelengths = [450 500 550 600];

%% Type of plot to generate. Choose between {'1Dslice', '2Dmap'}
plotView = '1Dslice';

%% ------- Effect of changing the spatial samples -------
% 261*2 gives a sampling interval of 1.003 c/deg
spatialSamplesList = [261 * 2 + 1, 201];
% Constants
referencePupilSize = 16.5;
sampleDomain = 'psf';

% Figure stuff
figNo = 1;
legends = {};
figName = sprintf(strcat("'ref pupil plane size' = %2.2fmm, ", ...
    "'spatial samples' = VARY, 'sample interval domain' = '%s'"), ...
    referencePupilSize, sampleDomain);

for k = 1:numel(spatialSamplesList)
    spatialSamples = spatialSamplesList(k);
    legends{numel(legends) + 1} = ...
        sprintf('sp. samples: %2.0f', spatialSamples);
    theWVF = createWVF(sampleDomain, spatialSamples, ...
        referencePupilSize, wavelengths);
    plotWVFsampling(figNo, theWVF, figName, k == 1, plotView, legends, k);
end

%% ------- Effect of changing the referencePupilSize -------
referencePupilSizes = [16.5 30];

% Constants
spatialSamples = 401;
sampleDomain = 'psf';

% Figure stuff
figNo = 2;
legends = {};
figName = sprintf(strcat("'ref pupil plane size' = VARY, 'spatial ", ...
    "samples' = %2.0f, 'sample interval domain' = '%s'"), ...
    spatialSamples, sampleDomain);

for k = 1:numel(referencePupilSizes)
    referencePupilSize = referencePupilSizes(k);
    legends{numel(legends) + 1} = ...
        sprintf('ref. pupil size: %2.1f', referencePupilSize);
    theWVF = createWVF(sampleDomain, spatialSamples, ...
        referencePupilSize, wavelengths);
    plotWVFsampling(figNo, theWVF, figName, k == 1, plotView, legends, k);
end

%% ------- Effect of changing the sample interval domain  -------
sampleDomains = {'psf', 'pupil'};

% Constants
spatialSamples = 201;
referencePupilSize = 16.5;

% Figure stuff
figNo = 3;
legends = {};
figName = sprintf(strcat("'ref pupil plane size' = %2.1f, ", ...
    "'spatial samples' = %2.0f, 'sample interval domain' = VARY"), ...
    referencePupilSize, spatialSamples);

for k = 1:numel(referencePupilSizes)
    sampleDomain = sampleDomains{k};
    legends{numel(legends) + 1} = ...
        sprintf('sample domain: ''%s''', sampleDomain);
    theWVF = createWVF(sampleDomain, spatialSamples, ...
        referencePupilSize, wavelengths);
    plotWVFsampling(figNo, theWVF, figName, k == 1, plotView, legends, k);
end

end

%% Method to generate a custom WVF object
function theWVF = createWVF(sampleDomain, spatialSamples, ...
    referencePupilSize, wavelengths)
% Helper function to generate a custom wavefront object
%
% Syntax:
%   theWVF = createWVF(sampleDomain, spatialSamples, ...
%       referencePupilSize, wavelengths)
%
% Description:
%    This is an embeddded helper function within t_wavefrontSampling to
%    allow for the creation of a custom wavefront object.
%
% Inputs:
%    sampleDomain       - Vector. String vector containing the possible
%                         domain categories.
%    spatialSample      - Numeric. The number of spatial samples.
%    referencePupilSize - Numeric. The pupil size in mm.
%    wavelengths        - Vector. A vector containing the desired
%                         wavelengths in nanometers.
%
% Outputs:
%    theWVF             - Struct. The created wavefront structure.
%
% Optional key/value pairs:
%    None.
%
if (~isempty(referencePupilSize))
    theWVF = wvfCreate('calc wavelengths', wavelengths, ...
        'ref pupil plane size', referencePupilSize, ...
        'spatial samples', spatialSamples, ...
        'sample interval domain', sampleDomain);
else
    theWVF = wvfCreate('calc wavelengths', wavelengths, ...
        'spatial samples', spatialSamples, ...
        'sample interval domain', sampleDomain);
end
% Compute the PSF
theWVF = wvfComputePSF(theWVF);
end

%% Method to extract and plot the PSF/OTF out of the WVF
function plotWVFsampling(figNo, theWVF, figName, resetFigure, plotView, ...
    legends, plotID)
% Helper function to plot the WVF's PSF and OTF
%
% Syntax:
%   plotWVFsampling(figNo, theWVF, figName, resetFigure, plotView, ...
%       legends, plotID)
%
% Description:
%    This is an embeddded helper function within t_wavefrontSampling to
%    demonstrate the PSF and OTF plots for the provided wavefront.
%
% Inputs:
%    figNo       - Numeric. An integer denoting the figure number.
%    theWVF      - Struct. The wavefront structure.
%    figName     - String. A string containing the figure name.
%    resetFigure - Boolean. A boolean in integer format on whether or not
%                  to display the figure.
%    plotView    - String. A string describing the plot view. Options are
%                  '1Dslide' and '2Dmap'.
%    legends     - Vector. A vector containing the strings describing the
%                  pictured plots.
%    plotID      - Numeric. An integer denoting the plot ID.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
assert(ismember(plotView, {'1Dslice', '2Dmap'}), ...
    'Incorrect plotView: choose between ''1Dslice'' and ''2Dmap''.');
markerSizes = [10 6];
colors = brewermap(4, 'Set2');

if (strcmp(plotView, '2Dmap'))
    hFig = figure(figNo + plotID * 100);
    legendsToUse = {};
    resetFigure = true;
else
    hFig = figure(figNo);
    legendsToUse = legends;
end
if (resetFigure)
    clf;
    set(hFig, ...
        'Position', [10 + 100 * figNo, 10 + 100 * figNo, 1500, 800], ...
        'Color', [1 1 1], 'Name', figName);
end

wavelengths = wvfGet(theWVF, 'wave');
subplotPosVectors = NicePlot.getSubPlotPosVectors('rowsNum', 2, ...
    'colsNum', numel(wavelengths), 'heightMargin', 0.1, ...
    'widthMargin', 0.04, 'leftMargin', 0.04, 'rightMargin', 0.00, ...
    'bottomMargin', 0.07, 'topMargin', 0.04);

psfRange = 15;   % arc min
otfRange = 400;  % cycles/deg

for iWave = 1:numel(wavelengths)
    targetWavelength = wavelengths(iWave);
    psf = wvfGet(theWVF, 'psf', targetWavelength);
    otf = wvfGet(theWVF, 'otf', targetWavelength);
    psfSupport = wvfGet(theWVF, 'psf angular samples', ...
        'min', targetWavelength);
    otfSupport = wvfGet(theWVF, 'otf support', 'um', targetWavelength) ...
        * wvfGet(theWVF, 'um per degree');

    subplot('Position', subplotPosVectors(1, iWave).v);
    hold on
    plotPSF(psf, psfSupport, psfRange, plotView, ...
        targetWavelength, iWave, legendsToUse, ...
        squeeze(colors(plotID, :)), markerSizes(plotID));

    subplot('Position', subplotPosVectors(2, iWave).v);
    hold on
    plotOTF(otf, otfSupport, otfRange, plotView, ...
        targetWavelength, iWave, legendsToUse, ...
        squeeze(colors(plotID, :)), markerSizes(plotID));
end
end

%% Method to plot the PSF
function plotPSF(psf, support, range, plotView, targetWavelength, col, ...
    legends, markerColor, markerSize)
% A helper function to plot the PSF called by plotWVFsampling
%
% Syntax:
%   plotPSF(psf, support, range, plotView, targetWavelength, col, ...
%       legends, markerColor, markerSize)
%
% Description:
%    This is an embeddded helper function within t_wavefrontSampling to
%    plot only the PSF, which is called by plotWVFsampling.
%
% Inputs:
%    psf              - Matrix. A matrix of the point spread function.
%    support          - Vector. The vector containing the PSF support.
%    range            - Numeric. The range for the PSF.
%    plotView         - String. A string describing the plot view. Options
%                       are '1Dslice' and '2Dmap'.
%    targetWavelength - Numeric. The wavelength to use, in nanometers.
%    col              - Numeric. The figure column to display in.
%    legends          - Vector. The legend(s) to use.
%    markerColor      - Vector. RGB values indicating marker color.
%    markerSize       - Numeric. The desired marker size.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
psf = psf / max(psf(:));
centerPos = floor(size(psf, 1) / 2) + 1;
if (isempty(range)), range = max(support); end
switch plotView
    case '1Dslice'
        plot(support, squeeze(psf(centerPos, :)), 'ko-', ...
            'LineWidth', 1.5, 'MarkerSize', markerSize, ...
            'MarkerFaceColor', markerColor);
        set(gca, 'XLim', [0 range], 'YLim', [0 1]);
    case '2Dmap'
        imagesc(support, support, psf);
        axis 'xy';
        axis 'image'
        set(gca, 'XLim', [-range range], 'YLim', [-range range], ...
            'CLim', [0 1]);
end
if ~isempty(legends), legend(legends, 'Location', 'NorthEast'); end
grid on;
box on;
set(gca, 'FontSize', 14);
if (col > 1), set(gca, 'YTickLabel', {}); end
xlabel('space (arc min)');
title(sprintf('PSF (%2.0f nm)', targetWavelength));
end

%% Method to plot the OTF
function plotOTF(otf, support, range, plotView, targetWavelength, col, ...
    legends, markerColor, markerSize)
% A helper function to plot the OTF called by plotWVFsampling
%
% Syntax:
%   plotOTF(otf, support, range, plotView, targetWavelength, col, ...
%       legends, markerColor, markerSize)
%
% Description:
%    This is an embeddded helper function within t_wavefrontSampling to
%    plot only the OTF, which is called by plotWVFsampling.
%
% Inputs:
%    otf              - Matrix. A matrix of the optical transfer function.
%    support          - Vector. The vector containing the OTF support.
%    range            - Numeric. The range for the OTF.
%    plotView         - String. A string describing the plot view. Options
%                       are '1Dslice' and '2Dmap'.
%    targetWavelength - Numeric. The wavelength to use, in nanometers.
%    col              - Numeric. The figure column to display in.
%    legends          - Vector. The legend(s) to use.
%    markerColor      - Vector. RGB values indicating marker color.
%    markerSize       - Numeric. The desired marker size.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
otfMag = fftshift(abs(otf));
centerPos = floor(size(otfMag, 1) / 2) + 1;
if (isempty(range)), range = max(support); end
switch plotView
    case '1Dslice'
        plot(support, squeeze(otfMag(centerPos, :)), 'ko-', ...
            'LineWidth', 1.5, 'MarkerSize', markerSize, ...
            'MarkerFaceColor', markerColor);
        set(gca, 'XLim', [0.3 range], 'YLim', [0 1], 'XScale', 'log');
    case '2Dmap'
        imagesc(support, support, otfMag);
        axis 'xy';
        axis 'image'
        set(gca, 'XLim', [-range range], 'YLim', [-range range], ...
            'CLim', [0 1]);
end
if ~isempty(legends), legend(legends, 'Location', 'NorthEast'); end
grid on;
box on
set(gca, 'FontSize', 14);
if (col > 1), set(gca, 'YTickLabel', {}); end
xlabel('sp. freq. (c/deg)');
title(sprintf('OTF (%2.0f nm)', targetWavelength));
end
