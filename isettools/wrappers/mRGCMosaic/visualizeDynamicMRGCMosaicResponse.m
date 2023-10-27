function visualizeDynamicMRGCMosaicResponse(theMRGCMosaic, visualizedMRGCindices, ...
    noiseFreeResponseSequence, noisyResponseInstancesSequence, responseTemporalSupportSeconds, generateVideo)

    %% Visualize the mRGCMosaic
    visualizedWidthDegs  = theMRGCMosaic.sizeDegs(1)*1.01;
    visualizedHeightDegs = theMRGCMosaic.sizeDegs(2)*1.01;
    domainVisualizationLimits(1:2) = theMRGCMosaic.eccentricityDegs(1) + visualizedWidthDegs  * 0.5*[-1 1];
    domainVisualizationLimits(3:4) = theMRGCMosaic.eccentricityDegs(2) + visualizedHeightDegs * 0.5*[-1 1];

    hFig = figure(1); clf;
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


    %% Visualize the time course of activation of select mRGC responses in a 3D perspective plot
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1020 715], 'Color', [1 1 1]);

    visualizedMRGCXcoords = squeeze(theMRGCMosaic.rgcRFpositionsDegs(visualizedMRGCindices,1));

    subplotHeight = 0.1;
    subplotWidth = 0.4;
    dX = (1.0-subplotWidth)/(numel(visualizedMRGCXcoords)+5);
    dY = (1.0-subplotHeight)/(numel(visualizedMRGCXcoords)+5);
    faceColor = [1 0.85 0.85];
    edgeColor = 'none';
    faceAlpha = 0.75;
    lineWidth = 1.5;
    lineStyle = '-';
    for iRGC = 1:numel(visualizedMRGCXcoords)
        perspectiveScaleFactor = 0.7+0.3*iRGC/numel(visualizedMRGCXcoords);
        perspectiveScaleFactor2 = 0.5+1.6*iRGC/numel(visualizedMRGCXcoords);
        gain = 2 * perspectiveScaleFactor2;
        responseRange = gain*[min(noiseFreeResponseSequence(:)) max(noiseFreeResponseSequence(:))];

        xCoord = dX + iRGC * dX*0.2;
        yCoord = 1 - subplotHeight - iRGC*dY*perspectiveScaleFactor;
        ax = axes('Position', [xCoord yCoord subplotWidth*perspectiveScaleFactor2 subplotHeight*perspectiveScaleFactor2 ]);
        
        mm = 0; % min(gain*squeeze(noiseFreeResponseSequence(1,:,visualizedMRGCindices(iRGC))));
        shadedAreaBetweenTwoLines(ax, ...
            responseTemporalSupportSeconds*1000, ...
            gain*noiseFreeResponseSequence(1,:,visualizedMRGCindices(iRGC)), ...
            responseTemporalSupportSeconds*0+mm, ...
            faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);

        hold(ax, 'on')
        plot(ax, responseTemporalSupportSeconds*1000, responseTemporalSupportSeconds*0, 'k-');
        plot(ax, responseTemporalSupportSeconds*1000, gain*noiseFreeResponseSequence(1,:,visualizedMRGCindices(iRGC)), ...
            'r-', 'LineWidth', 1.5, 'MarkerFaceColor', [1 0.5 0.5]);
        

        set(ax, 'Color', 'none', 'XColor', 'none', 'YColor', [0.5 0.5 0.5], ...
            'XTick', [], 'YTick', [], 'YLim', responseRange, 'FontSize', 16);
        if (iRGC == numel(visualizedMRGCXcoords))
            set(ax, 'XColor', [0.5 0.5 0.5])
        end
        box(ax, 'off')
        
        if (iRGC == numel(visualizedMRGCXcoords))
            set(ax, 'XColor', [0 0 0], 'XTick', 0:100:20000);
            xlabel(ax, 'time (msec)')
        end
    end



    %% Visualize the peak response amplitudes of the select mRGCs
    hFig = figure(3); clf;
    set(hFig, 'Position', [10 10 1020 715], 'Color', [1 1 1]);

    maxResponseAmplitudes = max(abs(noiseFreeResponseSequence(1,:,visualizedMRGCindices)), [], 2);
size(maxResponseAmplitudes)

visualizedMRGCXcoords
    plot(visualizedMRGCXcoords, squeeze(maxResponseAmplitudes), 'ro');
    set(ax, 'XLim', domainVisualizationLimits(1:2), ...
            'YLim', [-1.0 1.0], 'FontSize', 16);
    grid(ax, 'on');
    xlabel('space (degs)')
    ylabel('mRGC response')

    pause


    if (1==2)
    %% Visualize the mean spatiotemporal response as a video
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [1 1 1]);
    ax = subplot('Position', [0.035 0.03 0.96 0.96]);

    videoFileName = 'FullMRGCMosaicNoiseFreeResponseToDriftingGrating';
    videoOBJ = VideoWriter(videoFileName, 'Uncompressed AVI');
    videoOBJ.FrameRate = 10;
    videoOBJ.open();

    for iFrame = 1:size(noiseFreeResponseSequence,2)
        noiseFreeResponse = noiseFreeResponseSequence(1,iFrame,:);

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
    videoOBJ.close();
    exportVideo(videoFileName);
    

    %% Visualize a noisy spatiotemporal response instance
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1900 1050], 'Color', [1 1 1]);
    ax = subplot('Position', [0.035 0.03 0.96 0.96]);

    videoFileName = 'FullMRGCMosaicNoisyResponseInstancesToDriftingGrating';
    videoOBJ = VideoWriter(videoFileName, 'Uncompressed AVI');
    videoOBJ.FrameRate = 10;
    videoOBJ.open();

    visualizedTrials = size(noisyResponseInstancesSequence,1);
    visualizedTrials = 1;

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
    exportVideo(videoFileName);


end

    


end


function exportVideo(videoFileName)
    sysCommand = sprintf('/opt/homebrew/bin/ffmpeg -i %s.avi -c:v libx264 -crf 22 -pix_fmt yuv420p %s.mp4', ...
        videoFileName, videoFileName);
    system(sysCommand);

    sysCommand = sprintf('rm %s.avi ', videoFileName);
    system(sysCommand);
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
