function centerAndSurroundConePoolingLineWeightingFunctions(pdfExportSubDir, figNo, figPos, ...
	spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
	centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('axesToRenderIn', [], @(x)(isempty(x)||(isa(x, 'handle'))));
    p.addParameter('domainVisualizationLimits', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isnumeric(x))));
    p.addParameter('compositeInsteadOfComponent', false, @islogical);
    p.addParameter('flipXYaxes', false, @islogical);
    p.addParameter('noGrid', false, @islogical);

    p.parse(varargin{:});
    axesToRenderIn = p.Results.axesToRenderIn;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    compositeInsteadOfComponent = p.Results.compositeInsteadOfComponent;
    flipXYaxes = p.Results.flipXYaxes;
    noGrid = p.Results.noGrid;
    
	switch (whichMeridian)
		case 'horizontal'
            xyCenter = spatialSupportCenterDegs(1);
		    centerLineWeightingProfile = centerLineWeightingFunctions.xProfile;
		    surroundLineWeightingProfile = surroundLineWeightingFunctions.xProfile;
		    xLabelString = 'eccentricity, x (degs)';
		    pdfFileName = 'subregionConePoolingHorizontalLineWeightingFunctions.pdf';
		    
		case 'vertical'
            xyCenter = spatialSupportCenterDegs(2);
		    centerLineWeightingProfile = centerLineWeightingFunctions.yProfile;
		    surroundLineWeightingProfile = surroundLineWeightingFunctions.yProfile;
		    xLabelString = 'eccentricity, y (degs)';
    		pdfFileName = 'subregionConePoolingVerticalLineWeightingFunctions.pdf';
    	otherwise
    		error('Unknown subregion line weighting meridian: ''%s''. Expected either ''horizontal'' or ''vertical''.', whichMeridian);
    end

    spatialSupportRangeArcMin = 3 * spatialSupportTickSeparationArcMin;
    maxXY = round(spatialSupportRangeArcMin/2);
    spatialSupportDegs = (-maxXY:0.05:maxXY)/60;
    spatialSupportXYDegs(:,1) = spatialSupportCenterDegs(1) + spatialSupportDegs;
    spatialSupportXYDegs(:,2) = spatialSupportCenterDegs(2) + spatialSupportDegs;
    dx = (spatialSupportDegs(end)-spatialSupportDegs(1))*0.05;
    XLims = xyCenter + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];
    if (~isempty(domainVisualizationLimits))
        XLims = domainVisualizationLimits;
    end

    maxX = max(centerLineWeightingFunctions.xProfile.amplitude(:));
    maxY = max(centerLineWeightingFunctions.yProfile.amplitude(:));
    maxProfile = max([maxX maxY]);
    YTicks = maxProfile*(-0.25:0.25:1.0);
    YLims = maxProfile*[-0.4 1.0];

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure no left axis label');

    if (isempty(axesToRenderIn))
        % Initialize figure
        hFig = figure(figNo); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        set(hFig, 'Position', [figPos(1) figPos(2) ff.figureSize(1) ff.figureSize(2)]);
        ax = theAxes{1,1};
    else
        ax = axesToRenderIn;
    end

    % Plot
    idx = find(centerLineWeightingProfile.amplitude>0.0005*max(centerLineWeightingProfile.amplitude(:)));
    [~,iidx] = sort(centerLineWeightingProfile.spatialSupportDegs(idx), 'ascend');
    xx = centerLineWeightingProfile.spatialSupportDegs(idx(iidx));
    yy = centerLineWeightingProfile.amplitude(idx(iidx));

    if (flipXYaxes)
        xflip = yy;
        yflip = xx;
        xx = xflip;
        yy = yflip;
    end

    RGCMosaicAnalyzer.visualize.xyDataAsShadedArea(ax, xx, yy, 0, [0.9 0.1 0.3], 'none', 0.5, ff.lineWidth);
    hold(ax, 'on');
    plot(ax, xx, yy, 'r-', 'Color', [0.9 0.1 0.3]*0.5, 'LineWidth', ff.lineWidth);
    
    
    idx = find(surroundLineWeightingProfile.amplitude>0.001*max(surroundLineWeightingProfile.amplitude(:)));
    [~,iidx] = sort(surroundLineWeightingProfile.spatialSupportDegs(idx), 'ascend');
    xx = surroundLineWeightingProfile.spatialSupportDegs(idx(iidx));
    yy = -surroundLineWeightingProfile.amplitude(idx(iidx));

    if (flipXYaxes)
        xflip = yy;
        yflip = xx;
        xx = xflip;
        yy = yflip;
    end

    RGCMosaicAnalyzer.visualize.xyDataAsShadedArea(ax, xx, yy, 0,  [0.3 0.5 0.65], 'none', 0.5, ff.lineWidth);
    plot(ax, xx, yy, 'b-', 'Color', [0.3 0.5 0.65]*0.5, 'LineWidth', ff.lineWidth);
    

    hold(ax, 'off');
    axis(ax, 'square');

    switch (whichMeridian)
        case 'horizontal'
            xTicksDegs = spatialSupportCenterDegs(1) + (-3:3)*spatialSupportTickSeparationArcMin/60;
            if (~isempty(domainVisualizationTicks))
                xTicksDegs = domainVisualizationTicks;
            end

            if (flipXYaxes)
                xlabel(ax, 'y-integrated cone weights');
            else
                ylabel(ax, 'y-integrated cone weights');
            end

        case 'vertical'
            xTicksDegs = spatialSupportCenterDegs(2) + (-3:3)*spatialSupportTickSeparationArcMin/60;
            if (~isempty(domainVisualizationTicks))
                xTicksDegs = domainVisualizationTicks;
            end
            if (flipXYaxes)
                xlabel(ax, 'x-integrated cone weights');
            else
                ylabel(ax, 'x-integrated cone weights');
            end
    end

    if (spatialSupportTickSeparationArcMin/60 > 1-100*eps)
       xTickLabels = sprintf('%2.0f\n', xTicksDegs);
    elseif (spatialSupportTickSeparationArcMin/60>= 0.1-100*eps)
       xTickLabels = sprintf('%2.1f\n', xTicksDegs);
    elseif (spatialSupportTickSeparationArcMin/60 > 0.01-100*eps)
       xTickLabels = sprintf('%2.2f\n', xTicksDegs);
    else
       xTickLabels = sprintf('%2.3f\n', xTicksDegs);
    end

    if (flipXYaxes)
        set(ax, ...
            'XLim', YLims, 'YLim', XLims, ...
            'YTick',  YTicks, 'YTick', xTicksDegs, ...
            'YTickLabel', xTickLabels, ...
            'XTickLabel', {});
    else
        set(ax, ...
            'XLim', XLims, 'YLim', YLims, ...
            'XTick', xTicksDegs, 'YTick', YTicks, ...
            'XTickLabel', xTickLabels, ...
            'YTickLabel', {});
    end

    if (flipXYaxes)
        ylabel(ax, xLabelString);
        set(ax,'YAxisLocation','right', 'XDir', 'reverse')
    else
        xlabel(ax, xLabelString);
    end

    if (noGrid)
        ff.grid = 'off';
    end

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);


    if (isempty(axesToRenderIn))
        % Export figure
        % OLD Way
        %theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        %thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, pdfFileName);

        thePDFfileName = fullfile(pdfExportSubDir, pdfFileName);
        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    end
end