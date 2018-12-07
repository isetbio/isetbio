function visualizePSF(theOI, targetWavelength, psfRangeArcMin, varargin)
p = inputParser;
p.addParameter('axesHandle', [], @ishandle);
% Parse input
p.parse(varargin{:});
axesHandle = p.Results.axesHandle;

psfRangeArcMin = 0.5*psfRangeArcMin;
psfTicksMin = 2*(-8:2:8);
if (psfRangeArcMin <= 10)
    psfTicks = psfTicksMin;
elseif (psfRangeArcMin <= 20)
    psfTicks = 2*psfTicksMin;
elseif (psfRangeArcMin <= 40)
    psfTicks = 4*psfTicksMin;
elseif (psfRangeArcMin <= 50)
    psfTicks = 5*psfTicksMin; 
elseif (psfRangeArcMin <= 100)
    psfTicks = 10*psfTicksMin; 
elseif (psfRangeArcMin <= 200)
    psfTicks = 20*psfTicksMin;
elseif (psfRangeArcMin <= 400)
    psfTicks = 40*psfTicksMin; 
end
optics = oiGet(theOI, 'optics');
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
    axesHandle = subplot(1,1,1);
    fontSize = 20;
else
    fontSize = 12;
end
axes(axesHandle);

contourLevels = 0:0.05:1.0;
contourf(xSupportMinutes, ySupportMinutes, wavePSF/max(wavePSF(:)), contourLevels);
hold on;
plot(xSupportMinutes, psfRangeArcMin*(2*psfSlice-1), 'c-', 'LineWidth', 4.0);
plot(xSupportMinutes, psfRangeArcMin*(2*psfSlice-1), 'b-', 'LineWidth', 1.0);
axis 'image'; axis 'xy';
grid on
set(gca, 'XLim', psfRangeArcMin*1.05*[-1 1], 'YLim', psfRangeArcMin*1.05*[-1 1], 'CLim', [0 1], ...
            'XTick', psfTicks, 'YTick', psfTicks);
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0], 'Color', [0 0 0]);
xlabel('\it space (arc min)');
set(gca, 'FontSize', fontSize);
cmap = brewermap(1024, 'greys');
colormap(cmap);
title(sprintf('%s PSF (%2.0f nm)', optics.name, targetWavelength));
end
