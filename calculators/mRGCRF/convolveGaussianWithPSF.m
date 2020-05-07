function convolveGaussianWithPSF()
    
    for subjectID = 1:10
    eccXrange = [-10 -10];
    eccYrange = [0 0];
    deltaEcc = 1;
    
     % Compute the PSFs at the desired eccentricities
    [hEcc, vEcc, thePSFs, thePSFsupportDegs] = CronerKaplanRGCModel.psfAtEccentricity(subjectID, ...
        eccXrange, eccYrange, deltaEcc);
    
    gaussianRadiusDegs = 0.02;
    [X,Y] = meshgrid(thePSFsupportDegs,thePSFsupportDegs);
    gaussianSubregion = exp(-(X/gaussianRadiusDegs).^2).*exp(-(Y/gaussianRadiusDegs).^2);
    
    imagedSubregion = conv2(gaussianSubregion, squeeze(thePSFs(1, 1, 1,:,:)), 'same');
    
    
    % Visualize them
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 800 864]);
    
    cMap = brewermap(512, 'greys');
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 2, ...
           'heightMargin',  0.001, ...
           'widthMargin',    0.005, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.02, ...
           'topMargin',      0.01);
       
    ax = subplot('Position', subplotPosVectors(1,1).v);
    imagesc(ax, thePSFsupportDegs, thePSFsupportDegs, squeeze(thePSFs(1, 1, 1,:,:))); hold(ax, 'on');
    plot(ax, [0 0], [-1 1], 'r-'); plot(ax, [-1 1], [0 0], 'r-');
    set(ax, 'XLim', 0.5*max(thePSFsupportDegs)*[-1 1], 'YLim', 0.5*max(thePSFsupportDegs)*[-1 1]);
    axis(ax, 'square');  axis(ax, 'xy');
    set(ax, 'XTickLabel', {}, 'YTickLabel', {});
    colormap(ax,cMap);
    title(ax, sprintf('PSF at %2.1f, 2.1%f degs', eccXrange(1), eccYrange(1)));
    
    ax = subplot('Position', subplotPosVectors(1,2).v);
    contourf(ax, thePSFsupportDegs, thePSFsupportDegs, gaussianSubregion, 10);  hold(ax, 'on');
    plot(ax, [0 0], [-1 1], 'r-'); plot(ax, [-1 1], [0 0], 'r-');
    set(ax, 'XLim', 0.5*max(thePSFsupportDegs)*[-1 1], 'YLim', 0.5*max(thePSFsupportDegs)*[-1 1]);
    axis(ax, 'square');  axis(ax, 'xy');
    set(ax, 'XTickLabel', {}, 'YTickLabel', {});
    axis(ax, 'square');
    colormap(ax,cMap);
    title(ax, 'Gaussian subregion');
    
    ax = subplot('Position', subplotPosVectors(2,2).v);
    contourf(ax, thePSFsupportDegs, thePSFsupportDegs, imagedSubregion, 10);
    hold(ax, 'on');
    plot(ax, [0 0], [-1 1], 'r-'); plot(ax, [-1 1], [0 0], 'r-');
    set(ax, 'XLim', 0.5*max(thePSFsupportDegs)*[-1 1], 'YLim', 0.5*max(thePSFsupportDegs)*[-1 1]);
    axis(ax, 'square'); axis(ax, 'xy');
    set(ax, 'XTickLabel', {}, 'YTickLabel', {});
    axis(ax, 'square');
    colormap(ax,cMap);
    title(ax, 'imaged Gaussian subregion');
    
    ax = subplot('Position', subplotPosVectors(2,1).v);
    contourf(ax, thePSFsupportDegs, thePSFsupportDegs, imagedSubregion-gaussianSubregion, 10);
     hold(ax, 'on');
    plot(ax, [0 0], [-1 1], 'r-'); plot(ax, [-1 1], [0 0], 'r-');
    set(ax, 'XLim', 0.5*max(thePSFsupportDegs)*[-1 1], 'YLim', 0.5*max(thePSFsupportDegs)*[-1 1]);
    axis(ax, 'square');  axis(ax, 'xy');
    set(ax, 'XTickLabel', {}, 'YTickLabel', {});
    axis(ax, 'square');
    colormap(ax,cMap);
    title(ax, 'difference');
    
    drawnow;
    end
end

