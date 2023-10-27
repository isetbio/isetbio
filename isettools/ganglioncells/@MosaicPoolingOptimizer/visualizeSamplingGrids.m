function visualizeSamplingGrids(obj)

    conesNumPooledByTheRFcenters = unique(obj.conesNumPooledByTheRFcenterGrid);

    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x1 medium');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    obj.theRGCMosaic.visualize('figureHandle', hFig, ...
                            'axesHandle',theAxes{1,1}, ...
                            'identifyInputCones', true, ...
                            'backgroundColor', [0.5 0.5 0.5], ...
                            'fontAngle', ff.axisFontAngle, ...
                            'fontSize', ff.fontSize, ...
                            'domainVisualizationTicks', struct(...
                                'x', 0-5:0.5:5, 'y', -5:0.5:5));

    set(theAxes{1,1}, 'TickDir', 'both');

    % Font size
    set(theAxes{1,1}, 'FontSize', ff.fontSize);

    % axis color and width
    set(theAxes{1,1}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
   
    if (~isempty(obj.exportedPDFFolder))
        pdfFileName = fullfile(obj.exportedPDFFolder,'rgcmosaic.pdf');
    else
        pdfFileName = 'rgcmosaic.pdf';
    end
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);


    for iConesNumPooled = 1:numel(conesNumPooledByTheRFcenters)

        theIdentifiedConeIndices = [];

        hFig = figure(iConesNumPooled+1); clf;
        ff = MSreadyPlot.figureFormat('1x1 medium');
        theAxes = MSreadyPlot.generateAxes(hFig,ff);
        set(hFig, 'Color', [1 1 1]);


        % Determine spatial grid coords for this # of center cones
        conesNumPooled = conesNumPooledByTheRFcenters(iConesNumPooled);
        gridNodesList = find(obj.conesNumPooledByTheRFcenterGrid == conesNumPooled);

        LcenterRGCs = obj.targetRGCindicesWithLconeMajorityCenter(gridNodesList);
        LcenterRGCsString = 'L-center RGCs: ';
        for i = 1:numel(LcenterRGCs)
            fprintf('L-center grid node %d (%d cone input RF center) is located at position: %f %f\n', ...
                gridNodesList(i),...
                conesNumPooled, ...
                obj.theRGCMosaic.rgcRFpositionsDegs(LcenterRGCs(i),1), ...
                obj.theRGCMosaic.rgcRFpositionsDegs(LcenterRGCs(i),2));

            if (i < numel(LcenterRGCs))
                LcenterRGCsString = sprintf('%s%d ,', LcenterRGCsString, LcenterRGCs(i));
            else
                LcenterRGCsString = sprintf('%s%d', LcenterRGCsString, LcenterRGCs(i));
            end

            inputConeIndicesForThisTargetRGC = find(squeeze(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,LcenterRGCs(i))) > 0.0001);
            theIdentifiedConeIndices = cat(1, theIdentifiedConeIndices, inputConeIndicesForThisTargetRGC);
        end

        McenterRGCs = obj.targetRGCindicesWithMconeMajorityCenter(gridNodesList);
        McenterRGCsString = 'M-center RGCs: ';
        for i = 1:numel(McenterRGCs)
            fprintf('M-center grid node %d (%d cone input RF center) is located at position: %f %f\n', ...
                gridNodesList(i),...
                conesNumPooled, ...
                obj.theRGCMosaic.rgcRFpositionsDegs(McenterRGCs(i),1), ...
                obj.theRGCMosaic.rgcRFpositionsDegs(McenterRGCs(i),2));

            if (i < numel(McenterRGCs))
                McenterRGCsString = sprintf('%s%d ,', McenterRGCsString, McenterRGCs(i));
            else
                McenterRGCsString = sprintf('%s%d', McenterRGCsString, McenterRGCs(i));
            end
            inputConeIndicesForThisTargetRGC = find(squeeze(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,McenterRGCs(i))) > 0.0001);
            theIdentifiedConeIndices = cat(1, theIdentifiedConeIndices, inputConeIndicesForThisTargetRGC);
        end

    
        obj.theRGCMosaic.inputConeMosaic.visualize(...
                            'figureHandle', hFig, ...
                            'axesHandle',theAxes{1,1}, ...
                            'labelConesWithIndices', theIdentifiedConeIndices, ...
                            'backgroundColor', [0.5 0.5 0.5], ...
                            'fontAngle', ff.axisFontAngle, ...
                            'fontSize', ff.fontSize, ...
                            'domainVisualizationTicks', struct(...
                                'x', 0-5:0.5:5, 'y', -5:0.5:5));

        %                    'plotTitle', sprintf('%d-cone center sanpling positions\n%s\n%s', conesNumPooled, LcenterRGCsString, McenterRGCsString), ...
        %                    'plotTitleColor', [1 1 0.5]);

        
        set(theAxes{1,1}, 'TickDir', 'both');

        % Font size
        set(theAxes{1,1}, 'FontSize', ff.fontSize);

        % axis color and width
        set(theAxes{1,1}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
   
        pdfFileName = sprintf('%d-coneCenterGrid', conesNumPooled);
        if (~isempty(obj.exportedPDFFolder))
            pdfFileName = fullfile(obj.exportedPDFFolder,pdfFileName);
        end
        NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

    end


end

function old(obj)
   hFig = figure(2); clf;
   set(hFig, 'Color', [1 1 1]);

   conesNumPooledByTheRFcenters = unique(obj.conesNumPooledByTheRFcenterGrid);

   for iConesNumPooled = 1:numel(conesNumPooledByTheRFcenters)
        % Determine spatial grid coords for this # of center cones
        conesNumPooled = conesNumPooledByTheRFcenters(iConesNumPooled);
        gridNodesList = find(obj.conesNumPooledByTheRFcenterGrid == conesNumPooled);

        theLconeMajoritySpatialGrid = obj.theRGCMosaic.rgcRFpositionsDegs(obj.targetRGCindicesWithLconeMajorityCenter(gridNodesList),:);
        theMconeMajoritySpatialGrid = obj.theRGCMosaic.rgcRFpositionsDegs(obj.targetRGCindicesWithMconeMajorityCenter(gridNodesList),:);

        plotSamplingGridOnTopOfMosaic = false;
        ax = subplot(2,numel(conesNumPooledByTheRFcenters), iConesNumPooled);
        plotMultiFocalSpatialSamplingGrid(obj, ax, ...
            theLconeMajoritySpatialGrid, ...
            gridNodesList, ...
            plotSamplingGridOnTopOfMosaic, ...
            [1 0 0], ...
            sprintf('spatial sampling grid for %d L-center cones', conesNumPooled));

        ax = subplot(2,numel(conesNumPooledByTheRFcenters), iConesNumPooled+numel(conesNumPooledByTheRFcenters));
        plotMultiFocalSpatialSamplingGrid(obj, ax, ...
            theMconeMajoritySpatialGrid, ...
            gridNodesList, ...
            plotSamplingGridOnTopOfMosaic, ...
            [0 0.5 0], ...
            sprintf('spatial sampling grid for %d M-center cones', conesNumPooled));
    end

end

function plotMultiFocalSpatialSamplingGrid(obj, ax, spatialSamplingGrid, gridNodesList, plotSamplingGridOnTopOfMosaic, color, plotTitle)

    if (plotSamplingGridOnTopOfMosaic)
        obj.theRGCMosaic.visualize(...
            'axesHandle', ax, ...
            'samplingGrid', spatialSamplingGrid, ...
            'samplingGridOutlineColor', [1 1 0],...
            'samplingGridFillColor', [0 0 1],...
            'plotTitle', plotTitle);
    else
        plotMultiFocalSpatialSamplingGridIndices(obj, ax, spatialSamplingGrid, gridNodesList, color, plotTitle);
    end

    
end

function plotMultiFocalSpatialSamplingGridIndices(obj, ax, spatialSamplingGrid, gridNodesList, fontColor, plotTitle)
    centerDegs = obj.theRGCMosaic.eccentricityDegs;
    sizeDegs = obj.theRGCMosaic.sizeDegs;

    xx = centerDegs(1) + sizeDegs(1)/2*[-1 -1 1  1 -1];
    yy = centerDegs(2) + sizeDegs(2)/2*[-1  1 1 -1 -1];
    plot(ax,xx,yy,'k-', 'LineWidth', 1);
    hold(ax, 'on');
    for iNode = 1:numel(gridNodesList)
        text(ax, spatialSamplingGrid(iNode,1), spatialSamplingGrid(iNode,2), sprintf('%d', gridNodesList(iNode)), ...
            'FontSize', 14, 'Color', fontColor);
    end
    set(ax, 'FontSize', 16);
    axis(ax, 'equal')
    title(ax, plotTitle);
    box(ax, 'off');
    grid(ax, 'on');
    set(ax, 'XLim', centerDegs(1) + sizeDegs(1)*0.5*[-1 1] + [-0.1 0.1], ...
            'YLim', centerDegs(2) + sizeDegs(2)*0.5*[-1 1] + [-0.1 0.1]);
 
    drawnow;
end
