function visualizePSF(theOI, targetWavelength, psfRangeArcMin, varargin)
p = inputParser;
p.addParameter('axesHandle', [], @ishandle);
p.addParameter('withSuperimposedMosaic', [], @(x)(isa(x, 'coneMosaicHex')));
p.addParameter('figureTitle', '', @ischar);
% Parse input
p.parse(varargin{:});
axesHandle = p.Results.axesHandle;
theMosaic = p.Results.withSuperimposedMosaic;
figureTitle = p.Results.figureTitle;

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

% Extract support in arcmin
psfSupportMicrons = opticsGet(optics,'psf support','um');
if (isfield(optics, 'micronsPerDegree'))
    micronsPerDegree = optics.micronsPerDegree;
else
    focalLengthMeters = opticsGet(optics, 'focalLength');
    focalLengthMicrons = focalLengthMeters * 1e6;
    micronsPerDegree = focalLengthMicrons * tand(1);
end

xGridMinutes = 60*psfSupportMicrons{1}/micronsPerDegree;
yGridMinutes = 60*psfSupportMicrons{2}/micronsPerDegree;
xSupportMinutes = xGridMinutes(1,:);
ySupportMinutes = yGridMinutes(:,1);

% Extract slice through horizontal meridian
[~,idx] = min(abs(ySupportMinutes));
psfSlice = wavePSF(idx,:)/max(wavePSF(:));

if (isempty(axesHandle))
    figure(); clf;
    axesHandle = subplot('Position', [0.15 0.2 0.9 0.7]);
    fontSize = 20;
else
    fontSize = 12;
end
axes(axesHandle);

if (~isempty(theMosaic))
   theMosaic.visualizeGrid('axesHandle', axesHandle, ...
       'backgroundColor', 0.6*[1 1 1], ...
       'labelConeTypes', false);
   drawnow;
   hold on; 
   % transform minutes to meters
   xSupportMinutes = xSupportMinutes / 60 * theMosaic.micronsPerDegree * 1e-6;
   ySupportMinutes = ySupportMinutes / 60 * theMosaic.micronsPerDegree * 1e-6;
   psfRangeArcMin = psfRangeArcMin / 60 * theMosaic.micronsPerDegree * 1e-6;
   psfTicks = psfTicks / 60 * theMosaic.micronsPerDegree * 1e-6;
end

contourLevels = 0:0.1:1.0;
if (~isempty(theMosaic))
    contour(xSupportMinutes, ySupportMinutes, wavePSF/max(wavePSF(:)), contourLevels, ...
        'Color', 'c', 'LineWidth', 4);
    contour(xSupportMinutes, ySupportMinutes, wavePSF/max(wavePSF(:)), contourLevels, ...
        'Color', 'b', 'LineWidth', 1);
    hold on;
    plot(xSupportMinutes, psfRangeArcMin*(psfSlice-1), '-', 'Color', [1 1 0], 'LineWidth', 3.0);
    plot(xSupportMinutes, psfRangeArcMin*(psfSlice-1), '-', 'Color', [1 0.5 0], 'LineWidth', 1.0);
else
    contourf(xSupportMinutes, ySupportMinutes, wavePSF/max(wavePSF(:)), contourLevels, ...
        'Color', [0 0 0], 'LineWidth', 1.5);

    hold on;
    plot(xSupportMinutes, psfRangeArcMin*(psfSlice-1), 'c-', 'LineWidth', 3.0);
    plot(xSupportMinutes, psfRangeArcMin*(psfSlice-1), 'b-', 'LineWidth', 1.0);
end

axis 'image'; axis 'xy';
grid on; box on
set(gca, 'XLim', psfRangeArcMin*1.05*[-1 1], 'YLim', psfRangeArcMin*1.05*[-1 1], 'CLim', [0 1], ...
            'XTick', psfTicks, 'YTick', psfTicks, 'XTickLabel', psfTickLabels, 'YTickLabel', psfTickLabels);
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0]);
xlabel('\it space (arc min)');
ylabel('');
set(gca, 'FontSize', fontSize);
cmap = brewermap(1024, 'greys');
colormap(cmap);
if (isempty(figureTitle))
    title(gca, sprintf('%s\n(%2.0f nm) %dmm pupil', optics.name, targetWavelength, pupilDiameterMM));
else
    title(gca, figureTitle);
end
end
