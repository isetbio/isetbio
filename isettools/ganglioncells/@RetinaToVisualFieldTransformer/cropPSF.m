function cropPSF(obj,maxSpatialSupportDegs)

    figure(100); clf;
    subplot(2,2,1)
    imagesc(obj.thePSFData.supportX, obj.thePSFData.supportY, obj.thePSFData.data);
    XYLims = max(abs(obj.thePSFData.supportY(:)))*[-1 1];
    axis 'image';
    set(gca, 'XLim', XYLims, 'YLim', XYLims);
    title('psf, before cropping')

    subplot(2,2,2)
    imagesc(obj.theCircularPSFData.supportX, obj.theCircularPSFData.supportY, obj.theCircularPSFData.data);
    axis 'image';
    set(gca, 'XLim', XYLims, 'YLim', XYLims);
    title('circular psf, before cropping')

    % Reduce spatial support of the PSF to decrease compute time
    idx = find(abs(obj.thePSFData.supportX) < maxSpatialSupportDegs*60);
    idy = find(abs(obj.thePSFData.supportY) < maxSpatialSupportDegs*60);

    obj.thePSFData.supportX = obj.thePSFData.supportX(idx);
    obj.thePSFData.supportY = obj.thePSFData.supportY(idy);
    obj.thePSFData.data = obj.thePSFData.data(idy,idx);
    
    obj.theCircularPSFData.supportX = obj.theCircularPSFData.supportX(idx);
    obj.theCircularPSFData.supportY = obj.theCircularPSFData.supportY(idy);
    obj.theCircularPSFData.data = obj.theCircularPSFData.data(idy,idx);

    subplot(2,2,3)
    imagesc(obj.thePSFData.supportX, obj.thePSFData.supportY, obj.thePSFData.data);
    axis 'image';
    set(gca, 'XLim', XYLims, 'YLim', XYLims);

    title('psf, after cropping')

    subplot(2,2,4)
    imagesc(obj.theCircularPSFData.supportX, obj.theCircularPSFData.supportY, obj.theCircularPSFData.data);
    axis 'image';
    set(gca, 'XLim', XYLims, 'YLim', XYLims);

    title('circular psf, after cropping')

    colormap(gray)
end