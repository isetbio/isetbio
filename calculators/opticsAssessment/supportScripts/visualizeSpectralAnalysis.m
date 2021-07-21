function hFig = visualizeSpectralAnalysis(figNo,  horizontalEcc, verticalEcc, ...
    psfSupportArcMin, psfImage, coneMosaicImage, ...
    spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, ...
    psfImageSpectrum, psfImageSpectrumRoatated,  ...
    coneMosaicImageSpectrum, coneSpectrumSlice, ...
    spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, effectivePSFSpectrumSliceX, ...
    spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, effectivePSFSpectrumSliceY, ...
    coneCutoffSF, psfXCutoffSF, psfYCutoffSF, mosaicNyquistFrequencyCPD)

    maxVisualizedSpatialFrequency = 200;
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', 3, ...
        'colsNum', 2, ...
        'heightMargin',  0.06, ...
        'widthMargin',    0.03, ...
        'leftMargin',     0.04, ...
        'rightMargin',    0.01, ...
        'bottomMargin',   0.04, ...
        'topMargin',      0.01);

    hFig = figure(1000+figNo); clf;
    set(hFig, 'Position', [10 10 870 1300], 'Color', [1 1 1]);

    % The PSF
    ax = subplot('Position', subplotPosVectors(1,1).v);
    imagesc(ax,psfSupportArcMin, psfSupportArcMin, psfImage);
    axis(ax,'image')
    set(ax, 'XLim', 20*[-1 1], 'YLim', 20*[-1 1], 'XTick', -20:10:20, 'YTick', -20:10:20, 'FontSize', 14);
    ylabel('space,y (arc min)');
    xlabel('space, x (arc min)');
    title(ax, sprintf('PSF @ (%2.0f,%2.0f)', horizontalEcc, verticalEcc));
    colormap(ax, brewermap(1024, '*spectral'));

    if (~isempty(coneMosaicImage))
        %The cone aperture
        ax = subplot('Position', subplotPosVectors(1,2).v);
        imagesc(ax, psfSupportArcMin, psfSupportArcMin, coneMosaicImage);
        axis(ax, 'image')
        xlabel('space, x (arc min)');
        set(ax, 'XLim', 20*[-1 1], 'YLim', 20*[-1 1], 'XTick', -20:10:20, 'YTick', -20:10:20, 'FontSize', 14);
        title(ax, 'cone aperture');
        colormap(ax, brewermap(1024, '*greys'));
    end
    
    % The amplitude spectrum of the PSF
    ax = subplot('Position', subplotPosVectors(2,1).v);
    imagesc(ax,spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegree, log10(psfImageSpectrum.^2));
    axis(ax, 'image');
    set(ax, 'XLim', maxVisualizedSpatialFrequency*[-1 1], 'YLim', maxVisualizedSpatialFrequency*[-1 1], ...
        'XTick', -maxVisualizedSpatialFrequency:50:maxVisualizedSpatialFrequency, ...
        'YTick', -maxVisualizedSpatialFrequency:50:maxVisualizedSpatialFrequency, 'CLim', [-3 0], 'FontSize', 14);
    ylabel(ax,'spatial frequency, y (cpd)');
    title(ax,'PSF power spectrum');
    colormap(ax,brewermap(1024, '*greys'));

    if (~isempty(coneMosaicImageSpectrum))
        % The amplitude spectrum of the cone mosaic
        ax = subplot('Position', subplotPosVectors(2,2).v);
        imagesc(ax,spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegree, log10(coneMosaicImageSpectrum.^2));
        axis(ax, 'image')
        set(ax, 'XLim', maxVisualizedSpatialFrequency*[-1 1], 'YLim', maxVisualizedSpatialFrequency*[-1 1], ...
            'XTick', -maxVisualizedSpatialFrequency:50:maxVisualizedSpatialFrequency, ...
            'YTick', -maxVisualizedSpatialFrequency:50:maxVisualizedSpatialFrequency, 'CLim', [-3 0], 'FontSize', 14);
        title(ax,'cone aperture power spectrum');
        colormap(ax,brewermap(1024, '*greys'));
    end
    
    % The amplitude spectrum of the rotated PSF
    ax = subplot('Position', subplotPosVectors(3,1).v);
    imagesc(ax,spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegree, log10(psfImageSpectrumRoatated.^2));
    axis(ax, 'image')
    set(ax, 'XLim', maxVisualizedSpatialFrequency*[-1 1], 'YLim', maxVisualizedSpatialFrequency*[-1 1], ...
        'XTick', -maxVisualizedSpatialFrequency:50:maxVisualizedSpatialFrequency, ...
        'YTick', -maxVisualizedSpatialFrequency:50:maxVisualizedSpatialFrequency, 'CLim', [-3 0], 'FontSize', 14);
    xlabel(ax,'spatial frequency, x (cpd)');
    ylabel(ax,'spatial frequency, y (cpd)');
    title(ax,'rotated PSF power spectrum');
    colormap(ax,brewermap(1024, '*greys'));

    
    if (spectralSupportCyclesPerDegreePositiveX(1) < 0.3)
        spectralSupportCyclesPerDegreePositiveX(1) = 0.3;
    end
    if (spectralSupportCyclesPerDegreePositiveY(1) < 0.3)
        spectralSupportCyclesPerDegreePositiveY(1) = 0.3;
    end
    
    % The slices of the power spectra
    ax = subplot('Position', subplotPosVectors(3,2).v);
    area(ax, spectralSupportCyclesPerDegreePositiveX, coneSpectrumSlice.^2, 0.001, ...
        'FaceAlpha', 0.5, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5); hold on;
    
    
    % The PSF slices
    plot(ax, spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX.^2, 'b-', ...
        'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 1.0);
    plot(ax, spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY.^2, 'r-', ...
        'MarkerSize', 8, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.0);

    % The effective PSF slices 
    plot(ax, spectralSupportCyclesPerDegreePositiveX, effectivePSFSpectrumSliceX.^2, 'b-', ...
        'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 2);
    plot(ax, spectralSupportCyclesPerDegreePositiveY, effectivePSFSpectrumSliceY.^2, 'r-', ...
        'MarkerSize', 8, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 2);
    
    scatter(ax, psfXCutoffSF, 3/100, 14*14, 'o', 'MarkerFaceColor', [0.5 1 1], 'MarkerEdgeColor', [0 0 1], 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);
    scatter(ax, psfYCutoffSF, 3/100, 14*14, 'o', 'MarkerFaceColor', [1 0.3 0.7], 'MarkerEdgeColor', [1 0 0],  'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);
    
    plot(ax, mosaicNyquistFrequencyCPD*[1 1], [0.001 1], 'k-', 'LineWidth', 3.0);
    plot(ax, mosaicNyquistFrequencyCPD*[1 1], [0.001 1], 'g--', 'LineWidth', 1.5);
    set(ax, 'XLim', [0.3 maxVisualizedSpatialFrequency], 'XTick', [0.1 0.3 1 3 10 30 100 300], ...
        'YLim', [0.001 1], 'XScale', 'log', 'YScale', 'log', 'YTick', [0.001 0.003 0.01 0.03 0.1 0.3 1], 'FontSize', 14);
    legend(ax,{'cone aperture', 'psf_x', 'psf_y'}, 'Location', 'SouthWest');
    %title(ax,'power spectra slices');
    axis(ax, 'square');
    grid(ax,'on'); box(ax, 'off');
    xlabel(ax,'spatial frequency, x (cpd)');
    ylabel(ax,'power spectrum (norm.)');


    drawnow;
end