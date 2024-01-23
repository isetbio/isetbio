function visualizeDynamicMRGCMosaicResponse(figNo, theMRGCMosaic, visualizedMRGCindices, ...
    noiseFreeResponseSequence, noisyResponseInstancesSequence, responseTemporalSupportSeconds, ...
    exportVideo, videoFileName)

    %% Visualize the mRGCMosaic
    visualizedWidthDegs  = theMRGCMosaic.sizeDegs(1)*1.01;
    visualizedHeightDegs = theMRGCMosaic.sizeDegs(2)*1.01;
    domainVisualizationLimits(1:2) = theMRGCMosaic.eccentricityDegs(1) + visualizedWidthDegs  * 0.5*[-1 1];
    domainVisualizationLimits(3:4) = theMRGCMosaic.eccentricityDegs(2) + visualizedHeightDegs * 0.5*[-1 1];

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [1 1 1]);
    ax = subplot('Position', [0.035 0.03 0.96 0.96]);
    theMRGCMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'identifyInputCones', ~true, ...
        'identifyPooledCones',~true, ...
        'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
        'centerSubregionContourSamples', 24, ...
        'backgroundColor', [1 1 1], ...
        'labelRGCsWithIndices', visualizedMRGCindices, ...
        'labeledRGCsColor', [1 0 0], ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', struct('x', -10:10, 'y', -10:10), ...
        'plotTitle', '');


    %% Visualize the time course of activation of select mRGC responses
    hFig = figure(figNo+1); clf;
    set(hFig, 'Position', [10 10 750 1150], 'Color', [1 1 1]);

    
    visualizedMRGCXcoords = squeeze(theMRGCMosaic.rgcRFpositionsDegs(visualizedMRGCindices,1));
    trialNo = 1;
    spatioTemporalResponseOfSelectMRGCs = (squeeze(noiseFreeResponseSequence(trialNo,:,visualizedMRGCindices)))';
    activationRange = max(abs(noiseFreeResponseSequence(:)))*[-1 1];
    dt = responseTemporalSupportSeconds(2)-responseTemporalSupportSeconds(1);
    temporalSupportRange = [responseTemporalSupportSeconds(1)-dt/2 responseTemporalSupportSeconds(end)+dt/2]*1000;

    ax = subplot(2,1,1);
    imagesc(ax,responseTemporalSupportSeconds*1000, visualizedMRGCXcoords, spatioTemporalResponseOfSelectMRGCs);
    set(ax, 'XLim', temporalSupportRange);
    set(ax, 'CLim', activationRange);
    xlabel(ax,'time (msec)');
    set(ax, 'FontSize', 16, 'YDir', 'reverse');
    ylabel(ax,'mRGC cell x-position (degs)');
    c = colorbar(ax, 'northoutside');
    c.Label.String = 'mRGC activation';
    colormap(ax,brewermap(1024, '*greys'));

    ax = subplot(2,1,2);
    cLUT = brewermap(numel(visualizedMRGCXcoords), 'Spectral');
    for iCell = 1:numel(visualizedMRGCXcoords)
        lineColor = squeeze(cLUT(iCell,:));
        hold(ax, 'on')
        plot(ax, responseTemporalSupportSeconds*1000, spatioTemporalResponseOfSelectMRGCs(iCell,:), ...
            '-', 'LineWidth', 1.5, 'Color', lineColor);
        
    end

    colorbarTicks = prctile(visualizedMRGCXcoords(:), 0:20:100);
    c = colorbar(ax, 'northoutside');
    c.Label.String = 'mRGC x position';
    c.Ticks = 0:0.2:1;
    c.TickLabels = sprintf('%2.2f\n', colorbarTicks);
    colormap(ax,cLUT);

    set(ax, 'XLim', temporalSupportRange);
    set(ax, 'YLim', activationRange);
    xlabel(ax,'time (msec)');
    set(ax, 'FontSize', 16);
    ylabel(ax,'mRGC activation');
    drawnow


    if (exportVideo)
        %% Visualize the mean spatiotemporal response as a video
        hFig = figure(figNo+4); clf;
        set(hFig, 'Position', [10 10 1900 1050], 'Color', [1 1 1]);
        ax = subplot('Position', [0.035 0.03 0.96 0.96]);
    
        videoOBJ = VideoWriter(videoFileName, 'Uncompressed AVI');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();

        iTrial = 1;
        for iFrame = 1:size(noiseFreeResponseSequence,2)
            noiseFreeResponse = noiseFreeResponseSequence(iTrial,iFrame,:);
    
            theMRGCMosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', ax, ...
                'activation', noiseFreeResponse, ...
                'activationRange', max(abs(noiseFreeResponseSequence(:)))*[-1 1], ...
                'verticalActivationColorBarInside', true, ...
                'centerSubregionContourSamples', 24, ...
                'backgroundColor', [0 0 0], ...
                'labelRGCsWithIndices', visualizedMRGCindices, ...
                'labeledRGCsColor', [1 0 0], ...
                'domainVisualizationLimits', domainVisualizationLimits, ...
                'domainVisualizationTicks', struct('x', -10:10, 'y', -10:10), ...
                'plotTitle', sprintf('noise-free response (%2.3f sec)', responseTemporalSupportSeconds(iFrame)));
    
            drawnow;
            videoOBJ.writeVideo(getframe(hFig));
        end
        
    

        %% Visualize noisy spatiotemporal response instances
    
        visualizedTrials = size(noisyResponseInstancesSequence,1);

        for iTrial = 1:visualizedTrials
        for iFrame = 1:size(noisyResponseInstancesSequence,2)
            noiseFreeResponse = noisyResponseInstancesSequence(iTrial,iFrame,:);
    
            theMRGCMosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', ax, ...
                'activation', noiseFreeResponse, ...
                'activationRange', max(abs(noisyResponseInstancesSequence(:)))*[-1 1], ...
                'verticalActivationColorBarInside', true, ...
                'centerSubregionContourSamples', 24, ...
                'backgroundColor', [0 0 0], ...
                'labelRGCsWithIndices', visualizedMRGCindices, ...
                'labeledRGCsColor', [1 0 0], ...
                'domainVisualizationLimits', domainVisualizationLimits, ...
                'domainVisualizationTicks', struct('x', -10:10, 'y', -10:10), ...
                'plotTitle', sprintf('noisy response instance %d (%2.3f sec)', iTrial, responseTemporalSupportSeconds(iFrame)));
    
            drawnow;
            videoOBJ.writeVideo(getframe(hFig));
        end
        end
        
        videoOBJ.close();

        % Reformat video from AVI to mp4
        reformatVideo(videoFileName);
    end

end


function reformatVideo(videoFileName)
    if (exist('/opt/homebrew/bin/ffmpeg'))
        sysCommand = sprintf('/opt/homebrew/bin/ffmpeg -i %s.avi -c:v libx264 -crf 22 -pix_fmt yuv420p %s.mp4', ...
            videoFileName, videoFileName);
        system(sysCommand);
    
        sysCommand = sprintf('rm %s.avi ', videoFileName);
        system(sysCommand);
    end
end


function p = shadedAreaBetweenTwoLines(ax,x, y1, y2, ...
     faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)

    x = reshape(x, [1 numel(x)]);
    y1 = reshape(y1, [1 numel(x)]);
    y2 = reshape(y2, [1 numel(x)]);


    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    p = patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, 'edgeAlpha', 1.0,...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end
