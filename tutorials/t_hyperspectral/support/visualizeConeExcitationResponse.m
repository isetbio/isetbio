function visualizeConeExcitationResponse(theConeMosaic, theConeExcitations, theOI, emPath, figNo)
% Co-visualize the mosaic outline on the retinal image and the mosaic activation
% by the retinal image. If the emPath has more than 1 position, the
% generated figure is dynamically updated to show the mosaic excitation
% during the duration of the eye movement path
%
% 7/24/18  npc  Wrote it
%

    % Find range of response
    isomerizationsRange = [min(theConeExcitations(:)) max(theConeExcitations(:))];
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [100 100 840 900]);
    
    % Determine outline of cone mosaic on retinal inage
    coneApertureMicrons = theConeMosaic.patternSampleSize(1)*1e6;
    emPathOffsetDegs = emPath(1,1,:)*coneApertureMicrons/theConeMosaic.micronsPerDegree;
    mosaicOutline.x =-emPathOffsetDegs(1) + 0.5*theConeMosaic.fov(1)*[-1 1 1 -1 -1];
    mosaicOutline.y = emPathOffsetDegs(2) + 0.5*theConeMosaic.fov(2)*[-1 -1 1 1 -1];
    
    % Render part of the retinal image together with the cone mosaic outline
    ax = subplot('Position', [0.01 0.56, 0.98, 0.44]);
    retinalImageRGB = oiGet(theOI, 'rgb');
    [rows, cols, ~] = size(retinalImageRGB);
    spatialSizeDegs = [oiGet(theOI, 'wAngular') oiGet(theOI, 'hAngular') ]; 
    xSupportDegs = supportInDegs(cols, spatialSizeDegs(1));
    ySupportDegs = supportInDegs(rows, spatialSizeDegs(2)); 
    imagesc(ax,xSupportDegs, ySupportDegs, oiGet(theOI, 'rgb'));
    hold(ax, 'on');
    mosaicProfilePlot = plot(ax,mosaicOutline.x, mosaicOutline.y, 'y-', 'LineWidth', 1.5);
    hold(ax, 'off');
    axis(ax, 'image');
    set(ax, 'XLim', [-1 3], 'YLim', [-0.5 1.5], 'YTick', [], 'FontSize', 16);
    xlabel(ax, 'space (degs)', 'FontWeight', 'bold');
    
    ax2 = subplot('Position', [0.01 0.06, 0.98, 0.44]);
    eyeMovementsNum = size(emPath,2);
    if (eyeMovementsNum > 1)
        for emIndex = 1:eyeMovementsNum
            cla(ax2);
            theConeMosaic.renderActivationMap(ax2, ...
                squeeze(theConeExcitations(1,:,:,emIndex)), ...
                'signalRange', isomerizationsRange, ...
                'visualizedConeAperture', 'geometricArea', ...
                'backgroundColor', [0 0 0]);
             axis(ax2, 'ij');
             ylabel(ax2, '');
             set(ax2, 'YTick', [], 'FontSize', 16, ...
                  'XLim', 0.5*theConeMosaic.fov(1)*[-1 1]*theConeMosaic.micronsPerDegree*1e-6, ...
                  'YLim', 0.5*theConeMosaic.fov(2)*[-1 1]*theConeMosaic.micronsPerDegree*1e-6);

             % Update mosaic outline
             emPathOffsetDegs = emPath(1,emIndex,:)*coneApertureMicrons/theConeMosaic.micronsPerDegree;
             mosaicOutline.x =-emPathOffsetDegs(1) + 0.5*theConeMosaic.fov(1)*[-1 1 1 -1 -1];
             mosaicOutline.y = emPathOffsetDegs(2) + 0.5*theConeMosaic.fov(2)*[-1 -1 1 1 -1]; 
             set(mosaicProfilePlot, 'XData', mosaicOutline.x, 'YData', mosaicOutline.y)
             drawnow;
        end
    else
        theConeMosaic.renderActivationMap(ax2, ...
            theConeExcitations, ...
            'signalRange', isomerizationsRange, ...
            'visualizedConeAperture', 'geometricArea', ...
            'backgroundColor', [0 0 0]);
         axis(ax2, 'ij');
         ylabel(ax2, '');
         set(ax2, 'YTick', [], 'FontSize', 16, ...
             'XLim', 0.5*theConeMosaic.fov(1)*[-1 1]*theConeMosaic.micronsPerDegree*1e-6, ...
             'YLim', 0.5*theConeMosaic.fov(2)*[-1 1]*theConeMosaic.micronsPerDegree*1e-6);
    end
end
