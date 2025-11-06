function visualizePSF(theOI, targetWavelength, psfRangeArcMin, varargin)
% visualizePSF - Visualizes the Point Spread Function (PSF) of an optical image.
%
% (Written by CoPilot)
%
% Syntax:
%   visualizePSF(theOI, targetWavelength, psfRangeArcMin, 'ParameterName', ParameterValue, ...)
%
% Inputs:
%   theOI              - Optical image object (OI) containing optical properties.
%   targetWavelength   - Wavelength (in nm) at which to visualize the PSF.
%   psfRangeArcMin     - Range of the PSF visualization in arc minutes.
%
% Optional Name-Value Pair Arguments:
%   'axesHandle'                        - Handle to the axes for plotting (default: new axes).
%   'withSuperimposedMosaic'            - Cone mosaic object to superimpose (default: []).
%   'figureTitle'                        - Title for the figure (default: '').
%   'fontSize'                           - Font size for labels (default: [] for automatic).
%   'contourLevels'                      - Levels for contour plotting (default: 0.1:0.1:1.0).
%   'includePupilAndInFocusWavelengthInTitle' - Logical to include pupil and wavelength in title (default: true).
%   'noXLabel'                           - Logical to suppress X-axis label (default: false).
%   'noYLabel'                           - Logical to suppress Y-axis label (default: false).
%   'psfColorMap'                        - Custom colormap for PSF visualization (default: [] for greyscale).
%
% Example:
%   visualizePSF(theOI, 550, 10, 'figureTitle', 'PSF Visualization');
%
% See also: opticsGet, coneMosaicHex, brewermap

p = inputParser;
p.addParameter('axesHandle', [], @ishandle);
p.addParameter('withSuperimposedMosaic', [], @(x)(isa(x, 'coneMosaicHex')));
p.addParameter('figureTitle', '', @ischar);
p.addParameter('fontSize', []);
p.addParameter('contourLevels', 0.1:0.1:1.0);
p.addParameter('includePupilAndInFocusWavelengthInTitle', true, @islogical);
p.addParameter('noXLabel', false, @islogical);
p.addParameter('noYLabel', false, @islogical);
p.addParameter('psfColorMap', [], @isnumeric);

% Parse input
p.parse(varargin{:});
contourLevels = p.Results.contourLevels;
axesHandle = p.Results.axesHandle;
theMosaic = p.Results.withSuperimposedMosaic;
figureTitle = p.Results.figureTitle;
noXLabel = p.Results.noXLabel;
noYLabel = p.Results.noYLabel;
psfColorMap = p.Results.psfColorMap;

psfRangeArcMin = 0.5*psfRangeArcMin;
psfTicksMin = (-30:5:30);
if (psfRangeArcMin <= 2)
    psfTicks = (-3:0.5:3);
elseif (psfRangeArcMin <= 5)
    psfTicks = 0.2*psfTicksMin;
elseif (psfRangeArcMin <= 10)
    psfTicks = psfTicksMin;
elseif (psfRangeArcMin <= 20)
    psfTicks = 2*psfTicksMin;
elseif (psfRangeArcMin <= 40)
    psfTicks = 4*psfTicksMin;
elseif (psfRangeArcMin <= 50)
    psfTicks = 5*psfTicksMin; 
elseif (psfRangeArcMin <= 60)
    psfTicks = 6*psfTicksMin; 
elseif (psfRangeArcMin <= 100)
    psfTicks = 10*psfTicksMin; 
elseif (psfRangeArcMin <= 200)
    psfTicks = 20*psfTicksMin;
elseif (psfRangeArcMin <= 400)
    psfTicks = 40*psfTicksMin; 
end

if (psfRangeArcMin <= 2)
    psfTickLabels = sprintf('%2.1f\n', psfTicks);
else
    psfTickLabels = sprintf('%2.0f\n', psfTicks);
end

optics = oiGet(theOI, 'optics');
fLengthMeters = opticsGet(optics, 'focalLength');
fN = opticsGet(optics, 'fnumber');
pupilDiameterMM = fLengthMeters / fN * 1000;

wavelengthSupport = opticsGet(optics, 'wave');
[~,idx] = min(abs(wavelengthSupport-targetWavelength));
targetWavelength = wavelengthSupport(idx);

% Get PSF slice at target wavelength
wavePSF = opticsGet(optics,'psf data',targetWavelength);

% Extract PSF support
%
% The commented code here is what we were doing.  But the wavePSF
% provides the support in microns, so we just use that directly.
% The .xy field of wavePSF is (n,m,2), while the comment out opticsGet
% returned a cell array.  We fixed this below.  We are not entirely sure
% we have the x and y straight if support ever becomes rectangular.
%
% nSamp = opticsGet(optics,'otf size');
% if (nSamp(1) ~= nSamp(2))
%     error('Do not yet handle non-square psf support');
% end
% psfSupportMicrons = opticsGet(optics,'psf support','um',oiGet(theOI,'fsupport'),nSamp(1));
psfSupportMicrons = wavePSF.xy;

% Not sure when optics structures have a micronsPerDegree field
% Not desired to set or access this directly.
%
% @Nicolas, once we figure out wvf2oiSpecial we should do this
% with sets and gets and make sure we enforce consistency.
if (isfield(optics, 'micronsPerDegree'))
    micronsPerDegree = optics.micronsPerDegree;
else
    focalLengthMeters = opticsGet(optics, 'focalLength');
    focalLengthMicrons = focalLengthMeters * 1e6;
    micronsPerDegree = focalLengthMicrons * tand(1);
end

% If the x and y support are different, we should worry as
% we are not entirely sure of conventions.  Plus, we think
% that would be really bad as we'd then need the scene to have
% different number of pixels/deg in x and y.
%
% The transposing of psfSupport2 is a little squirrly.
psfSupport1 = squeeze(psfSupportMicrons(:,:,1));
psfSupport2 = squeeze(psfSupportMicrons(:,:,2))';
if (any(psfSupport1(:) ~= psfSupport2(:)))
    error('Time to think hard about structure of returned psf support');
end
xGridMinutes = 60*psfSupport1/micronsPerDegree;
yGridMinutes = 60*psfSupport2'/micronsPerDegree;
xSupportMinutes = xGridMinutes(1,:);
ySupportMinutes = yGridMinutes(:,1);

% Extract slice through horizontal meridian
[~,idx] = min(abs(ySupportMinutes));
psfSlice = wavePSF.psf(idx,:)/max(wavePSF.psf(:));

if (isempty(axesHandle))
    figure(); clf;
    axesHandle = subplot('Position', [0.15 0.2 0.9 0.7]);
    fontSize = 20;
else
    fontSize = 12;
end
if (~isempty(p.Results.fontSize))
    fontSize = p.Results.fontSize;
end

axes(axesHandle);

if (~isempty(theMosaic))
   theMosaic.visualizeGrid('axesHandle', axesHandle, ...
       'backgroundColor', 0.6*[1 1 1], ...
       'visualizedConeAperture', 'lightCollectingArea', ...
       'labelConeTypes', false);
   hold(axesHandle, 'on'); 
   % transform minutes to meters
   xSupportMinutes = xSupportMinutes / 60 * theMosaic.micronsPerDegree * 1e-6;
   ySupportMinutes = ySupportMinutes / 60 * theMosaic.micronsPerDegree * 1e-6;
   psfRangeArcMin = psfRangeArcMin / 60 * theMosaic.micronsPerDegree * 1e-6;
   psfTicks = psfTicks / 60 * theMosaic.micronsPerDegree * 1e-6;
end

if (isempty(psfColorMap))
    cmap = brewermap(1024, 'greys');
else
    cmap = psfColorMap;
end

colormap(cmap);

if (~isempty(theMosaic))
    transparentContourPlot(axesHandle, xSupportMinutes, ySupportMinutes, wavePSF/max(wavePSF(:)), ...
        [0.1 0.3 0.5 0.7 0.9], cmap);
    plot(axesHandle, xSupportMinutes, psfRangeArcMin*(psfSlice-1), '-', 'Color', [0.1 0.3 0.3], 'LineWidth', 4.0);
    plot(axesHandle, xSupportMinutes, psfRangeArcMin*(psfSlice-1), '-', 'Color', [0.3 0.99 0.99], 'LineWidth', 2);
else
    contourf(axesHandle,xSupportMinutes, ySupportMinutes, wavePSF.psf/max(wavePSF.psf(:)), contourLevels, ...
        'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);

    %imagesc(xSupportMinutes, ySupportMinutes, wavePSF.psf/max(wavePSF.psf(:)));
    hold(axesHandle, 'on');
    plot(axesHandle, xSupportMinutes, psfRangeArcMin*(psfSlice-1), 'c-', 'LineWidth', 3.0);
    plot(axesHandle, xSupportMinutes, psfRangeArcMin*(psfSlice-1), 'b-', 'LineWidth', 1.0);

end


xtickangle(axesHandle, 0);
axis(axesHandle, 'image'); 
axis(axesHandle, 'ij');
grid(axesHandle, 'on'); box(axesHandle,  'on');

set(axesHandle, 'XLim', psfRangeArcMin*1.05*[-1 1], 'YLim', psfRangeArcMin*1.05*[-1 1], 'CLim', [0 1], ...
                'XTick', psfTicks, 'YTick', psfTicks, 'XTickLabel', psfTickLabels, 'YTickLabel', psfTickLabels);
set(axesHandle, 'XColor', [0 0 0], 'YColor', [0 0 0]);
if (~noXLabel)
    xlabel(axesHandle,'space (arc min)');
else
    set(axesHandle, 'XTickLabel', {});
end

if (~noYLabel)
    ylabel(axesHandle,'space (arc min)');
else
    set(axesHandle, 'YTickLabel', {});
end
set(axesHandle, 'FontSize', fontSize);


if (isempty(figureTitle))
    if (~p.Results.includePupilAndInFocusWavelengthInTitle)
         title(axesHandle, sprintf('%s', oiGet(theOI,'name')), ...
            'FontWeight', 'Normal', 'FontSize', fontSize);
    else
        title(axesHandle, sprintf('(%2.0fnm,%dmm pupil)\n%s', targetWavelength, pupilDiameterMM, oiGet(theOI,'name')), ...
            'FontWeight', 'Normal', 'FontSize', fontSize);
    end
else
    title(axesHandle, figureTitle, 'FontWeight', 'Normal', 'FontSize', fontSize);
end
end

function transparentContourPlot(axesHandle, xSupportMinutes, ySupportMinutes, zData, zLevels, cmap)
    C = contourc(xSupportMinutes, ySupportMinutes, zData, zLevels);
    dataPoints = size(C,2);
    startPoint = 1;
    hold(axesHandle, 'on');
    while (startPoint < dataPoints)
        theLevel = C(1,startPoint);
        theLevelVerticesNum = C(2,startPoint);
        x = C(1,startPoint+(1:theLevelVerticesNum));
        y = C(2,startPoint+(1:theLevelVerticesNum));
        v = [x(:) y(:)];
        f = 1:numel(x);
        patch(axesHandle, 'Faces', f, 'Vertices', v, 'EdgeColor', 0.5*(1-theLevel)*[1 1 1], ...
            'FaceColor', cmap(round(theLevel*size(cmap,1)),:), ...
            'FaceAlpha', 0.4, ... %min([1 0.3+theLevel]), ...
            'LineStyle', '-', 'LineWidth', 1.0);
        startPoint = startPoint + theLevelVerticesNum+1;
    end
end
