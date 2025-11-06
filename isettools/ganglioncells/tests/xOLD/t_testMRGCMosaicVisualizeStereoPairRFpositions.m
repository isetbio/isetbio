% Script to visualize large scale mRGC mosaics and their spatial density
% Usage:
%{
    t_testMRGCMosaicVisualizeStereoPairRFpositions
%}

eccentricityDegs = [0 0]; sizeDegs = [64 64];

% sampling for the density map
densityMapSampleSpacingDegs = 1.0;

% Compute density map for left eye
theLeftEyeDensity = computeMRGCRFmapDensityMap('left eye', eccentricityDegs, sizeDegs, densityMapSampleSpacingDegs;

% Compute density map for right eye
theRightEyeDensity = computeMRGCRFmapDensityMap('right eye', eccentricityDegs, sizeDegs, densityMapSampleSpacingDegs);

theLeftEyeDensity.mapLog = log10(theLeftEyeDensity.map);
theRightEyeDensity.mapLog = log10(theRightEyeDensity.map);

% density isocontour levels
densityRange = prctile(theLeftEyeDensity.mapLog(:), [1 100]);
zLevels = 0.5:0.125:4.0;

densityTicks = zLevels(5:8:end);
densityTickLabels = sprintf('1E%1.0f\n',densityTicks);

% LUT for contour plot
theLUT = brewermap(1024, '*spectral');
contourLineColor = [0 0 0];

hFig = figure(2); clf;
set(hFig, 'Position', [10 10 1250 610], 'Color', [1 1 1]);
ax = subplot(1,2,1);
contourf(ax, theLeftEyeDensity.spatialSupportXdegs, theLeftEyeDensity.spatialSupportYdegs, theLeftEyeDensity.mapLog, ...
    zLevels, 'LineColor', contourLineColor);
axis(ax, 'xy'); axis(ax, 'image');
set(ax, 'XLim', 32*[-1 1], 'YLim', 32*[-1 1], 'XTick', -30:5:30, 'YTick', -30:5:30);
set(ax, 'Color', theLUT(1,:), 'CLim', densityRange, 'FontSize', 16);
xtickangle(ax, 0);
xlabel(ax, 'temporal           VF x-eccentricity (degs)            nasal', 'FontAngle', 'italic');
ylabel(ax, 'inferior          VF y-eccentricity (degs)       superior', 'FontAngle', 'italic');
cb = colorbar(ax, 'NorthOutside', ...
    'Ticks', densityTicks, 'TickLabels', densityTickLabels);
cb.Title.String = 'mRGC RF density (RFs/deg^2)';
title(ax, 'left eye');
colormap(ax, theLUT);

ax = subplot(1,2,2);
contourf(ax, theRightEyeDensity.spatialSupportXdegs, theRightEyeDensity.spatialSupportYdegs, theRightEyeDensity.mapLog, ...
    zLevels, 'LineColor', contourLineColor);
axis(ax, 'xy'); axis(ax, 'image');
set(ax, 'XLim', 32*[-1 1], 'YLim', 32*[-1 1], 'XTick', -30:5:30, 'YTick', -30:5:30);
set(ax, 'Color', theLUT(1,:), 'CLim', densityRange, 'FontSize', 16);
xtickangle(ax, 0);
xlabel(ax, 'nasal               VF x-eccentricity (degs)         temporal', 'FontAngle', 'italic');
ylabel(ax, 'inferior          VF y-eccentricity (degs)       superior', 'FontAngle', 'italic');
cb = colorbar(ax, 'NorthOutside', ...
    'Ticks', densityTicks, 'TickLabels', densityTickLabels);
cb.Title.String = 'mRGC RF density (RFs/deg^2)';
title(ax, 'right eye');
colormap(ax, theLUT);

NicePlot.exportFigToPDF('MRGCRFdensityMaps.pdf', hFig, 300);



function densityMap = computeMRGCRFmapDensityMap(whichEye, eccentricityDegs, sizeDegs, deltaDegs)

    sourceLatticeSizeDegs = 64;
    degsToMMsConversionFunction = @RGCmodels.Watson.convert.rhoDegsToMMs;
    MMsToDegsConversionFunction = @RGCmodels.Watson.convert.rhoMMsToDegs;

    eccMicrons = 1e3 * degsToMMsConversionFunction(eccentricityDegs);
    eccXrightMicrons  =  1e3 * degsToMMsConversionFunction(eccentricityDegs(1)+0.5*sizeDegs(1));
    eccXleftMicrons   =  1e3 * degsToMMsConversionFunction(eccentricityDegs(1)-0.5*sizeDegs(1));
    eccYtopMicrons    =  1e3 * degsToMMsConversionFunction(eccentricityDegs(2)+0.5*sizeDegs(2));
    eccYbottomMicrons =  1e3 * degsToMMsConversionFunction(eccentricityDegs(2)-0.5*sizeDegs(2));
    sizeMicrons = [eccXrightMicrons-eccXleftMicrons eccYtopMicrons-eccYbottomMicrons];

    % Import mRGC RF positions
    [mRGCRFposMicrons, mRGCRFposDegs] = retinalattice.import.finalMRGCPositions(...
             sourceLatticeSizeDegs, ...
             eccMicrons, ...
             sizeMicrons, ...
             whichEye, ...
             MMsToDegsConversionFunction);
    xyMin = min(mRGCRFposDegs,[],1);
    xyMax = max(mRGCRFposDegs,[],1);

    mRGCRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(mRGCRFposDegs);

    xPosDegs = eccentricityDegs(1):deltaDegs:(eccentricityDegs(1)+0.5*sizeDegs(1));
    xPosDegs = unique([-fliplr(xPosDegs) xPosDegs]);
    yPosDegs = eccentricityDegs(2):deltaDegs:(eccentricityDegs(2)+0.5*sizeDegs(2));
    yPosDegs = unique([-fliplr(yPosDegs) yPosDegs]);

    sampledPositionsDegs{1} = xPosDegs;
    sampledPositionsDegs{2} = yPosDegs ;    
    densityMap.map = cMosaic.densityMap(mRGCRFposDegs, mRGCRFspacingsDegs, sampledPositionsDegs);
    densityMap.spatialSupportXdegs = xPosDegs ;
    densityMap.spatialSupportYdegs = yPosDegs ;
end



