function visualizePSFAtWavelengthOnTopOfConeMosaic(hFig, ax, theConeMosaic, thePSF, theVisualizedWavelengthIndex, maxPSF, ...
    theTitle, varargin)
% VISUALIZEPSFATWAVELENGTHONTOPOFCONEMOSAIC Visualizes the Point Spread Function (PSF)
% at a specified wavelength on top of a cone mosaic.
%
% (Written by CoPilot)
%
% Syntax:
%   visualizePSFAtWavelengthOnTopOfConeMosaic(hFig, ax, theConeMosaic, thePSF,
%   theVisualizedWavelengthIndex, maxPSF, theTitle, varargin)
%
% Description:
%   This function overlays the Point Spread Function (PSF) at a given wavelength
%   onto a cone mosaic visualization. It allows for customization of the visualization
%   limits, ticks, and additional information strings.
%
% Inputs:
%   hFig                     - Handle to the figure where the visualization will be drawn.
%   ax                       - Handle to the axes where the cone mosaic will be displayed.
%   theConeMosaic            - An instance of the cone mosaic object to visualize.
%   thePSF                   - The Point Spread Function data to visualize.
%   theVisualizedWavelengthIndex - Index of the wavelength to visualize from the PSF data.
%   maxPSF                   - Maximum value for scaling the PSF visualization.
%   theTitle                 - Title for the visualization plot.
%
% Optional Name-Value Pair Arguments:
%   'XLimsArcMin'            - Limits for the x-axis in arc minutes (default: [-10 10]).
%   'YLimsArcMin'            - Limits for the y-axis in arc minutes (default: [-10 10]).
%   'XTicksArcMin'           - Custom x-axis ticks in arc minutes (default: []).
%   'YTicksArcMin'           - Custom y-axis ticks in arc minutes (default: []).
%   'withInfoString'         - String to display additional information on the plot (default: '').
%
% Example:
%   visualizePSFAtWavelengthOnTopOfConeMosaic(hFig, ax, coneMosaic, psfData, 5, 1, 'PSF Visualization', ...
%   'XLimsArcMin', [-15 15], 'YLimsArcMin', [-15 15], 'withInfoString', 'Sample Info');
%
% See also: CONEMOSAIC, PSF


% Parse input
p = inputParser;
% Optional params
p.addParameter('XLimsArcMin', [-10 10], @isnumeric);
p.addParameter('YLimsArcMin', [-10 10], @isnumeric);
p.addParameter('XTicksArcMin', [], @isnumeric);
p.addParameter('YTicksArcMin', [], @isnumeric);
p.addParameter('withInfoString', '', @ischar);

p.parse(varargin{:});
XLimsArcMin = p.Results.XLimsArcMin;
YLimsArcMin = p.Results.YLimsArcMin;
XTicksArcMin = p.Results.XTicksArcMin;
YTicksArcMin = p.Results.YTicksArcMin;
infoString = p.Results.withInfoString;

thePSFData.data = squeeze(thePSF.data(:,:,theVisualizedWavelengthIndex));
thePSFData.supportXdegs = thePSF.supportX/60;
thePSFData.supportYdegs = thePSF.supportY/60;

cMap = 0.9*brewermap(1024, 'blues') + 0.1*brewermap(1024, 'greens');
cMap(2,:) = cMap(1,:);
cMap(1,:) = [0.9 .9 .9];

theConeMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', ax, ...
    'conesAlpha', 0.2, ...
    'visualizedConeAperture', 'lightCollectingArea4sigma', ...
    'activation', zeros(1,1,theConeMosaic.conesNum), ...
    'activationColorMap', cMap, ...
    'activationRange', [0 1], ...
    'visualizedConeApertureThetaSamples', 20, ...
    'withSuperimposedPSF', thePSFData, ...
    'withSuperimposedPSFcontourLineColor', [0 0.2 0.9], ...
    'domainVisualizationLimits', [...
    theConeMosaic.eccentricityDegs(1) + XLimsArcMin(1)/60 ...
    theConeMosaic.eccentricityDegs(1) + XLimsArcMin(2)/60 ...
    theConeMosaic.eccentricityDegs(2) + YLimsArcMin(1)/60 ...
    theConeMosaic.eccentricityDegs(2) + YLimsArcMin(2)/60], ...
    'domainVisualizationTicks', struct('x', XTicksArcMin/60, 'y', YTicksArcMin/60), ...
    'backgroundColor', [1 1 1], ...
    'fontSize', 12, ...
    'noXLabel', true, ...
    'noYLabel', true);

if (~isempty(infoString))
    xo = theConeMosaic.eccentricityDegs(1) + XLimsArcMin(1)/60 + (XLimsArcMin(2)-XLimsArcMin(1))/60*0.03;
    yo = theConeMosaic.eccentricityDegs(2) + YLimsArcMin(1)/60 + (YLimsArcMin(2)-YLimsArcMin(1))/60*0.07;
    t = text(ax, xo,yo, infoString, 'Color', [0.3 0.3 0.5]);
    t.FontSize = 20;
    t.BackgroundColor = 'none';
    t.EdgeColor = 'none';
end

title(ax, theTitle, 'FontWeight', 'normal');
drawnow;
end
