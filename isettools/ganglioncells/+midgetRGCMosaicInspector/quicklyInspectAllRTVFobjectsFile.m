function quicklyInspectAllRTVFobjectsFile()
    
    dropboxDir = midgetRGCMosaicInspector.localDropboxPath();
    [file,path] = uigetfile(fullfile(dropboxDir, '*.mat'), ...
                        'Select an RTVF file');

    if (file == 0)
        return;
    end

    progressBar = waitbar(0.2,'Loading all computed R2VFT objects. Please wait ...');
    pause(.1);

    fName = fullfile(path,file);
    load(fName, 'theRTFVTobjList', ...
        'theConesNumPooledByTheRFcenterGrid', ...
        'theOpticsPositionGrid');

    for iObj = 1:numel(theRTFVTobjList)
        if (isempty(theRTFVTobjList{iObj}))
            fprintf(2, '**** There is no RTVFobj at grid position %d ***** \n', iObj);
        end
    end

    % Plot the sampling grid for all RTVF objects, separately for the different # of
    % center cones examined
    iObj = 1;
    while (isempty(theRTFVTobjList{iObj})) && (iObj < numel(theRTFVTobjList))
        iObj = iObj + 1;
    end

    if (iObj > numel(theRTFVTobjList))
        return;
    end

    x = theRTFVTobjList{iObj}.theConeMosaic.coneRFpositionsDegs(:,1);
    y = theRTFVTobjList{iObj}.theConeMosaic.coneRFpositionsDegs(:,2);
    xLims = [min(x) max(x)];
    yLims = [min(y) max(y)];
    spatialTicksX = linspace(xLims(1), xLims(2), 7);
    spatialTicksY = linspace(yLims(1), yLims(2), 7);
    spatialTicksX = -2:0.5:2;
    spatialTicksY = -2:0.5:2;


    hFig = figure(9000); clf;
    set(hFig, 'Position', [10 10 1500 500], 'Color', [1 1 1]);

    theDifferentConesNumPooled = sort(unique(theConesNumPooledByTheRFcenterGrid));
    markerTypes = {'o', 'x', 's', 'd', 'h', '*', 'p'};

    for iConesNumPooledIndex = 1:numel(theDifferentConesNumPooled)
        conesNumPooled = theDifferentConesNumPooled(iConesNumPooledIndex);
        idx = find(theConesNumPooledByTheRFcenterGrid == conesNumPooled);
        ax = subplot(1, numel(theDifferentConesNumPooled), iConesNumPooledIndex);
        scatter(ax, theOpticsPositionGrid(idx,1), theOpticsPositionGrid(idx, 2), 200, ...
            'LineWidth', 1.5, ...
            'Marker',  markerTypes{iConesNumPooledIndex});

        axis(ax, 'equal');
        grid(ax, 'on');
        xlabel(ax, 'space, x (degs)')
        ylabel(ax, 'space, y (degs)')
        set(ax, 'XLim', xLims, 'YLim', yLims);
        set(ax, 'XTick', spatialTicksX, 'XTickLabel', sprintf('%2.2f\n', spatialTicksX));
        set(ax, 'YTick', spatialTicksY, 'YTickLabel', sprintf('%2.2f\n', spatialTicksY));
        set(ax, 'FontSize', 16);

        title(ax, sprintf('sampling grid for %d cone(s) RF centers', conesNumPooled));
    end
    NicePlot.exportFigToPDF('samplingGrids.pdf', hFig, 300);

    for iObj = 1:numel(theRTFVTobjList)
        waitbar(0.3+(iObj/numel(theRTFVTobjList))*0.7,progressBar, ...
            sprintf('Processing R2VF obj %d of %d', iObj, numel(theRTFVTobjList)));
        pause(0.1);
        if (~isempty(theRTFVTobjList{iObj}))
            midgetRGCMosaicInspector.peekIntoRTVFobj(...
                theRTFVTobjList{iObj}, ...
                iObj, ...
                theOpticsPositionGrid, ...
                theConesNumPooledByTheRFcenterGrid, ...
                iObj*1000);
        end
    end

    close(progressBar);
end