function plotSpatialSamplingGrid(obj, ax, spatialSamplingGrid, plotTitle)
    
    xx = obj.mosaicCenterParams.positionDegs(1) + obj.mosaicCenterParams.sizeDegs(1)/2*[-1 -1 1  1 -1];
    yy = obj.mosaicCenterParams.positionDegs(2) + obj.mosaicCenterParams.sizeDegs(2)/2*[-1  1 1 -1 -1];
    plot(ax,xx,yy,'k-', 'LineWidth', 1);
    hold(ax, 'on');
    plot(ax, spatialSamplingGrid(:,1), spatialSamplingGrid(:,2), 'r+', 'MarkerSize', 12, 'LineWidth', 3.0);
    set(ax, 'FontSize', 16);
    axis(ax, 'equal')
    title(ax, plotTitle);
    box(ax, 'off');
    grid(ax, 'on');
    set(ax, 'XLim', obj.mosaicCenterParams.positionDegs(1) + obj.mosaicCenterParams.sizeDegs(1)*0.5*[-1 1] + [-0.1 0.1], ...
            'YLim', obj.mosaicCenterParams.positionDegs(2) + obj.mosaicCenterParams.sizeDegs(2)*0.5*[-1 1] + [-0.1 0.1]);
 
    drawnow;
end
