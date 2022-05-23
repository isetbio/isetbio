function visualizeSpatialSpectralRadiance(spatialSizeDegs, radianceWaveSupport, radiancePhotons, visualizedWaves, visualizedXRangeDegs, visualizedYRangeDegs, figNo, radianceType)
% Visualize spatial spectral radiance slices at select wavelengths
%
% 7/24/18  npc  Wrote it
%

    [rows, cols, ~] = size(radiancePhotons);
    xSupportDegs = supportInDegs(cols, spatialSizeDegs(1));
    ySupportDegs = supportInDegs(rows, spatialSizeDegs(2)); 
    photonsRange = [min(radiancePhotons(:)) max(radiancePhotons(:))];
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 3, ...
       'colsNum', 3, ...
       'heightMargin', 0.03, ...
       'widthMargin',  0.01, ...
       'leftMargin',  0.01, ...
       'rightMargin', 0.01, ...
       'bottomMargin', 0.05, ...
       'topMargin',  0.01);

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1060 700]);
    
    for k = 1:numel(visualizedWaves)
        [~,idx] = min(abs(radianceWaveSupport-visualizedWaves(k)));
        photonsAtWaveband = squeeze(radiancePhotons(:,:,idx));
        subplot('Position', subplotPosVectors(floor((k-1)/3)+1, mod(k-1,3)+1).v);
        imagesc(xSupportDegs, ySupportDegs, ...
            photonsAtWaveband, photonsRange);
        colormap(gray(1024));
        axis('image'); axis('ij');
        if (k == 7)
            xlabel('space (degs)', 'FontWeight', 'bold');
        else
            set(gca, 'XTick', []);
        end
        set(gca, 'XLim', visualizedXRangeDegs, 'YLim', visualizedYRangeDegs);
        colorbar
        title(sprintf('%s photon rate @ %d nm', radianceType, radianceWaveSupport(idx)));
        set(gca, 'FontSize', 14, 'YTick', []);
    end
end