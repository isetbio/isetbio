function visualize(obj, varargin)
    
    % Parse optional input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('clearAxesBeforeDrawing', true, @islogical);
    p.addParameter('labelRetinalMeridians', false, @islogical);
    p.addParameter('component', 'RF centers', @(x)ismember(x, {'RF centers'}));
    p.addParameter('centerSubregionContourSamples', 10, @isscalar);
    p.addParameter('contourGenerationMethod', 'ellipseFitToPooledConeApertureImage', @(x)(ismember(x, {'ellipseFitBasedOnLocalSpacing', 'contourOfPooledConeApertureImage','ellipseFitToPooledConeApertureImage'})));
    p.addParameter('activation', []);
    p.addParameter('samplingGrid', [], @(x)(isempty(x) || (isnumeric(x) && (size(x,2) == 2)) ));
    p.addParameter('samplingGridOutlineColor', [1 0 0], @(x)(isempty(x)||(numel(x)==3)));
    p.addParameter('samplingGridFillColor', [1 0.6 0.6], @(x)(isempty(x)||(numel(x)==3)));
    p.addParameter('identifiedConeAperture', 'geometricArea', @(x)ismember(x, ...
        {'lightCollectingArea', 'geometricArea', 'coneSpacing', ...
        'lightCollectingAreaCharacteristicDiameter', 'lightCollectingArea2sigma', 'lightCollectingArea4sigma', 'lightCollectingArea5sigma', 'lightCollectingArea6sigma'}));
    p.addParameter('identifiedConeApertureThetaSamples', [], @(x)(isempty(x) || isscalar(x)));
   
    p.addParameter('identifyInputCones', false, @islogical);
    p.addParameter('identifyPooledCones', false, @islogical);
    p.addParameter('identifiedInputConeIndices', [], @isnumeric);
    p.addParameter('pooledConesLineWidth', 1.0, @isscalar);
    p.addParameter('plotRFoutlines', true, @islogical);
    p.addParameter('plottedRFoutlineLineWidth', 1.0, @isscalar);
    p.addParameter('activationRange', [],@(x)((isempty(x))||(numel(x)==2)));
    p.addParameter('activationColorMap', [], @(x)(isempty(x)||(size(x,2) == 3)));
    p.addParameter('horizontalActivationColorBar', false, @islogical);
    p.addParameter('verticalActivationColorBar', false, @islogical);
    p.addParameter('horizontalActivationColorBarInside', false, @islogical);
    p.addParameter('verticalActivationColorBarInside', false, @islogical);
    p.addParameter('fontSize', 16, @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('fontAngle', 'normal', @(x)(ismember(lower(x), {'normal', 'italic'})));
    p.addParameter('labelRGCsWithIndices', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('labeledRGCsColor', [], @(x)(isempty(x)||(size(x,2) == 3)));
    p.addParameter('labeledRGCsLineWidth', 1.5, @isscalar);
    p.addParameter('colorbarFontSize', 16, @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('colorBarTickLabelPostFix', '', @ischar);
    p.addParameter('colorbarTickLabelColor',  [], @(x)(isempty(x)||((isvector(x))&&(numel(x) == 3))));
    p.addParameter('backgroundColor', [], @(x)( (ischar(x)&&((strcmp(x,'none'))||(strcmp(x,'mean of color map'))) ) || isempty(x) || ((isvector(x))&&(numel(x) == 3))));
    p.addParameter('plotTitle', '', @(x)(isempty(x) || ischar(x) || islogical(x)));
    p.addParameter('plotTitleColor', [0 0 0], @isnumeric);
    p.addParameter('plotTitleFontSize', 16, @isscalar);
    p.addParameter('domainVisualizationLimits', [], @(x)((isempty(x))||(numel(x)==4)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isstruct(x)&&((isfield(x, 'x'))&&(isfield(x,'y'))))));
    p.addParameter('superimposedRect', [], @(x)(isempty(x)||isstruct(x)));
    p.addParameter('superimposedRectColor', [], @(x)( isempty(x) || ((isvector(x))&&(numel(x) == 3))));
    p.parse(varargin{:});

    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandle;
    clearAxesBeforeDrawing = p.Results.clearAxesBeforeDrawing;
    labelRetinalMeridians = p.Results.labelRetinalMeridians;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    visualizedComponent = p.Results.component;
    centerSubregionContourSamples = p.Results.centerSubregionContourSamples;
    contourGenerationMethod = p.Results.contourGenerationMethod;
    identifiedConeAperture = p.Results.identifiedConeAperture;
    identifiedConeApertureThetaSamples = p.Results.identifiedConeApertureThetaSamples;
    identifiedInputConeIndices = p.Results.identifiedInputConeIndices;
    identifyInputCones = p.Results.identifyInputCones;
    identifyPooledCones = p.Results.identifyPooledCones;
    pooledConesLineWidth = p.Results.pooledConesLineWidth;
    plotRFoutlines = p.Results.plotRFoutlines;
    plottedRFoutlineLineWidth = p.Results.plottedRFoutlineLineWidth;
    activation = p.Results.activation;
    samplingGrid = p.Results.samplingGrid;
    samplingGridOutlineColor = p.Results.samplingGridOutlineColor;
    samplingGridFillColor = p.Results.samplingGridFillColor;
    activationRange = p.Results.activationRange;
    activationColorMap = p.Results.activationColorMap;
    colorBarTickLabelPostFix = p.Results.colorBarTickLabelPostFix;
    verticalColorBar = p.Results.verticalActivationColorBar;
    horizontalColorBar = p.Results.horizontalActivationColorBar;
    verticalColorBarInside = p.Results.verticalActivationColorBarInside;
    horizontalColorBarInside = p.Results.horizontalActivationColorBarInside;
    colorbarTickLabelColor = p.Results.colorbarTickLabelColor;
    backgroundColor = p.Results.backgroundColor;
    fontSize = p.Results.fontSize;
    fontAngle = p.Results.fontAngle;
    labelRGCsWithIndices = p.Results.labelRGCsWithIndices;
    labeledRGCsColor = p.Results.labeledRGCsColor;
    labeledRGCsLineWidth = p.Results.labeledRGCsLineWidth;
    colorbarFontSize = p.Results.colorbarFontSize;
    plotTitle = p.Results.plotTitle;
    plotTitleColor = p.Results.plotTitleColor;
    plotTitleFontSize = p.Results.plotTitleFontSize;
    superimposedRect = p.Results.superimposedRect;
    superimposedRectColor = p.Results.superimposedRectColor;

    % Generate the visualization cache
    xSupport = [];
    ySupport = []; 
    obj.generateVisualizationCache(xSupport, ySupport, centerSubregionContourSamples, contourGenerationMethod);

    % Determine X,Y limits
    if (isempty(domainVisualizationLimits))
        if (isfield(obj.visualizationCache, 'surroundConesXrange'))
            xRange = obj.visualizationCache.surroundConesXrange;
        else
            xRange(1) = min(squeeze(obj.inputConeMosaic.coneRFpositionsDegs(:,1)));
            xRange(2) = max(squeeze(obj.inputConeMosaic.coneRFpositionsDegs(:,1)));
        end

        if (xRange(2) == xRange(1))
            xRange = xRange(1) + 0.02*[-1 1];
        end
        if (isfield(obj.visualizationCache, 'surroundConesYrange'))
            yRange = obj.visualizationCache.surroundConesYrange;
        else
            yRange(1) = min(squeeze(obj.inputConeMosaic.coneRFpositionsDegs(:,2)));
            yRange(2) = max(squeeze(obj.inputConeMosaic.coneRFpositionsDegs(:,2)));
        end

        if (yRange(2) == yRange(1))
            yRange = yRange(1) + 0.02*[-1 1];
        end
        xx = xRange(2)-xRange(1);
        yy = yRange(2)-yRange(1);
        XLims(1) = xRange(1)-xx*0.02;
        XLims(2) = xRange(2)+xx*0.02;
        YLims(1) = yRange(1)-yy*0.02;
        YLims(2) = yRange(2)+yy*0.02;
        domainVisualizationLimits(1) = XLims(1);
        domainVisualizationLimits(2) = XLims(2);
        domainVisualizationLimits(3) = YLims(1);
        domainVisualizationLimits(4) = YLims(2);
    else
        XLims(1) = domainVisualizationLimits(1);
        XLims(2) = domainVisualizationLimits(2);
        YLims(1) = domainVisualizationLimits(3);
        YLims(2) = domainVisualizationLimits(4);
    end
    
    if (isempty(domainVisualizationTicks))
        xo = (XLims(1)+XLims(2))/2;
        xx = XLims(2)-XLims(1);
        yo = (YLims(1)+YLims(2))/2;
        yy = YLims(2)-YLims(1);
        ticksX = [XLims(1) xo XLims(2)];
        ticksY = [YLims(1) yo YLims(2)];
        
        if (xx > 10)
            domainVisualizationTicks.x = round(ticksX);
        elseif (xx > 5)
            domainVisualizationTicks.x = round(ticksX*10)/10;
        elseif (xx > 1)
            domainVisualizationTicks.x = round(ticksX*100)/100;
        else
            domainVisualizationTicks.x = round(ticksX*1000)/1000;
        end
        if (yy > 10)
            domainVisualizationTicks.y = round(ticksY);
        elseif (yy > 5)
            domainVisualizationTicks.y = round(ticksY*10)/10;
        elseif (yy > 1)
            domainVisualizationTicks.y = round(ticksY*100)/100;
        else
            domainVisualizationTicks.y = round(ticksY*1000)/1000;
        end 
    end
    
    switch (visualizedComponent)
        case 'RF centers'
            [hFig, ax] = visualizeRFcenters(obj, hFig, ax, clearAxesBeforeDrawing, ...
                labelRetinalMeridians, ...
                domainVisualizationTicks, domainVisualizationLimits, ...
                identifiedConeAperture, identifiedConeApertureThetaSamples, ...
                identifyInputCones, identifyPooledCones, identifiedInputConeIndices, ...
                pooledConesLineWidth, labelRGCsWithIndices, ...
                labeledRGCsColor, labeledRGCsLineWidth, ...
                plotRFoutlines, plottedRFoutlineLineWidth, superimposedRect, superimposedRectColor, ...
                activation, activationRange, activationColorMap, ...
                colorBarTickLabelPostFix, colorbarTickLabelColor, ...
                verticalColorBar, horizontalColorBar, colorbarFontSize, ...
                verticalColorBarInside, horizontalColorBarInside, ...
                backgroundColor, fontSize, fontAngle, ...
                plotTitle, plotTitleColor, plotTitleFontSize);

        otherwise
            error('Uknown visualized component: ''%s''.', visualizedComponent);
    end


    % Superimpose a sampling grid
    if (~isempty(samplingGrid))
       plot(ax, samplingGrid(:,1), samplingGrid(:,2), '^', ...
           'Color', samplingGridFillColor , 'LineWidth', 5.0, 'MarkerSize', 20);
       plot(ax, samplingGrid(:,1), samplingGrid(:,2), '^', ...
           'Color', samplingGridOutlineColor , 'LineWidth', 2.0, 'MarkerSize', 16);
    end

end


function [hFig, ax] = visualizeRFcenters(obj,hFig, ax, clearAxesBeforeDrawing, ...
        labelRetinalMeridians, ...
        domainVisualizationTicks, domainVisualizationLimits, ...
        identifiedConeAperture, identifiedConeApertureThetaSamples, ...
        identifyInputCones, identifyPooledCones, identifiedInputConeIndices, pooledConesLineWidth, labelRGCsWithIndices, ...
        labeledRGCsColor, labeledRGCsLineWidth, ...
        plotRFoutlines, plottedRFoutlineLineWidth, superimposedRect, superimposedRectColor, ...
        activation, activationRange, activationColorMap, ...
        colorBarTickLabelPostFix, colorbarTickLabelColor, ...
        verticalColorBar, horizontalColorBar, colorbarFontSize, ...
        verticalColorBarInside, horizontalColorBarInside, ...
        backgroundColor, fontSize, fontAngle, plotTitle, plotTitleColor,  plotTitleFontSize)

    
    if (isempty(ax))
        if (isempty(hFig))
            hFig = figure(); clf;
            if (isempty(obj.name))
                set(hFig, 'Color', [1 1 1], 'Position', [10 10 1300 900]);
            else
                set(hFig, 'Color', [1 1 1], 'Position', [10 10 1300 900], 'Name', obj.name);
            end
        end
        ax = subplot('Position', [0.05 0.06 0.91 0.91]);
    else
        if (clearAxesBeforeDrawing)
            cla(ax);
        end
    end

    if (~isempty(activation))
        % Visualize RF centers, color-coded with their activation level
        if (isempty(activationColorMap))
            linearRamp = linspace(0,1,1024);
            linearRamp = linearRamp(:);
            cMap = [linearRamp*0 linearRamp linearRamp*0];
        else
            cMap = activationColorMap;
        end
        
        if (ischar(backgroundColor) && strcmp(backgroundColor, 'mean of color map'))
            midRow = round(size(cMap,1)/2);
            backgroundColor = squeeze(cMap(midRow,:));
        elseif (isempty(backgroundColor))
            backgroundColor = squeeze(cMap(1,:));
        end
        
        if (isempty(activationRange))
            activationRange = [min(activation(:)) max(activation(:))];
        end

        activation = (activation - activationRange(1))/(activationRange(2)-activationRange(1));
        activation(activation<0) = 0;
        activation(activation>1) = 1;

        currentFacesNum = 0;
        for iRGC = 1:obj.rgcsNum
            newVerticesNum = obj.visualizationCache.rfCenterPatchData.verticesNumForRGC(iRGC);
            idx = currentFacesNum + (1:newVerticesNum);
            obj.visualizationCache.rfCenterPatchData.faceVertexCData(idx,:) = activation(iRGC);
            currentFacesNum = currentFacesNum + newVerticesNum;
        end
    else
        if (isempty(backgroundColor))
            backgroundColor = [1 1 1];
        end

        % All mRGC centers in gray
        cMap = [0 0 0; 0.5 0.5 0.5; 0 0 0];
        S.FaceVertexCData = obj.visualizationCache.rfCenterPatchData.faceVertexCData*0+0.5;
    end

    if (plotRFoutlines) || (~isempty(activation))
        % Plot the RFs
        S.Vertices = obj.visualizationCache.rfCenterPatchData.vertices;
        S.Faces = obj.visualizationCache.rfCenterPatchData.faces;
        S.FaceVertexCData = obj.visualizationCache.rfCenterPatchData.faceVertexCData;
    
        S.FaceColor = 'flat';
        if (~isempty(activation))
            S.EdgeColor = [0 0 0]';
        else
            S.FaceColor = [0.95 0.95 0.95]*0.6;
            S.EdgeColor = [0 0 0];
        end
        if (isempty(activation))
            S.FaceAlpha = 0.4;
        else
            S.FaceAlpha = 0.9;
        end

        S.EdgeAlpha = 1.0;
        S.LineWidth = plottedRFoutlineLineWidth;
        patch(S, 'Parent', ax)
    end


    if (identifyPooledCones)
        hold(ax, 'on')
        if (identifyInputCones)
            lConeInputLineColor = [1 0.1 0.5];
            mConeInputLineColor = [0.1 1 0.5];
            lineSegmentWidth = pooledConesLineWidth;
        else
            lConeInputLineColor = [0 0 0];
            mConeInputLineColor = [0 0 0];
            lineSegmentWidth = pooledConesLineWidth;
            % Put a single dot in all mRGC RF centers with a single input
            if (~isempty(obj.visualizationCache.rfCenterSingleConeInputDotPositions))
                plot(ax, obj.visualizationCache.rfCenterSingleConeInputDotPositions(:,1), ...
                         obj.visualizationCache.rfCenterSingleConeInputDotPositions(:,2), 'k.');
            end
        end
        
        % Render line segments from centroid to pulled cones
        renderPooledConesLineSegments(obj, ax, lConeInputLineColor, mConeInputLineColor, lineSegmentWidth);
    end

    % Identify input cones
    if (identifyInputCones)
        hold(ax, 'on')
        if (isempty(identifiedInputConeIndices))
            obj.inputConeMosaic.visualize(...
                'figureHandle', hFig, 'axesHandle', ax, ...
                'clearAxesBeforeDrawing', false, ...
                'visualizedConeAperture', identifiedConeAperture, ...
                'conesAlpha', 0.7, ...
                'conesEdgeAlpha', 0.7, ...
                'visualizedConeApertureThetaSamples', identifiedConeApertureThetaSamples, ...
                'labelRetinalMeridians', labelRetinalMeridians, ...
                'domainVisualizationTicks', domainVisualizationTicks, ...
                'domainVisualizationLimits', domainVisualizationLimits, ...
                'backgroundColor', backgroundColor);
        else
            obj.inputConeMosaic.visualize(...
                'figureHandle', hFig, 'axesHandle', ax, ...
                'clearAxesBeforeDrawing', false, ...
                'visualizedConeAperture', identifiedConeAperture, ...
                'labelConesWithIndices', identifiedInputConeIndices, ...
                'conesAlpha', 0.7, ...
                'conesEdgeAlpha', 0.7, ...
                'visualizedConeApertureThetaSamples', identifiedConeApertureThetaSamples, ...
                'labelRetinalMeridians', labelRetinalMeridians, ...
                'domainVisualizationTicks', domainVisualizationTicks, ...
                'domainVisualizationLimits', domainVisualizationLimits, ...
                'backgroundColor', backgroundColor);
        end
    else
        hold(ax, 'on')
        obj.inputConeMosaic.visualize(...
            'figureHandle', hFig, 'axesHandle', ax, ...
            'clearAxesBeforeDrawing', false, ...
            'labelCones', false, ...
            'visualizedConeApertureThetaSamples', identifiedConeApertureThetaSamples, ...
            'labelRetinalMeridians', labelRetinalMeridians, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'backgroundColor', backgroundColor);
    end

    if (~isempty(labelRGCsWithIndices))
        if (plotRFoutlines) || (~isempty(activation))
            hold(ax, 'on')
            for iRGC = 1:numel(labelRGCsWithIndices)
                theRGCindex = labelRGCsWithIndices(iRGC);
                if (iscell(obj.visualizationCache.rfCenterContourData{theRGCindex}))
                    S = obj.visualizationCache.rfCenterContourData{theRGCindex}{1};
                else
                    S = obj.visualizationCache.rfCenterContourData{theRGCindex};
                end

                S.FaceVertexCData = 0.5;
                S.FaceColor = 'flat';
                if (isempty(labeledRGCsColor))
                    S.EdgeColor = [1 1 0];
                else
                    S.EdgeColor = labeledRGCsColor;
                end

                S.FaceAlpha = 0.0;
                S.LineWidth = labeledRGCsLineWidth;
                
                S.LineStyle = '-';
                patch(S, 'Parent', ax);
            end
        end
    end


    if (~isempty(superimposedRect))
        hold(ax, 'on');
        if (isempty(superimposedRectColor))
            superimposedRectColor = [1 0 0];
        end

        x1 = superimposedRect.center(1) - 0.5*superimposedRect.xRange;
        x2 = superimposedRect.center(1) + 0.5*superimposedRect.xRange;
        y1 = superimposedRect.center(2) - 0.5*superimposedRect.yRange;
        y2 = superimposedRect.center(2) + 0.5*superimposedRect.yRange;
        xx = [x1 x1 x2 x2 x1];
        yy = [y1 y2 y2 y1 y1];
        if (all(superimposedRectColor(1:2) == [1 1]))
            plot(ax, xx,yy, '-', 'LineWidth', 3, 'Color', [0 0 0]);
        else
            plot(ax, xx,yy, '-', 'LineWidth', 3, 'Color', [1 1 0]);
        end

        plot(ax, xx,yy, '--', 'LineWidth', 3, 'Color', superimposedRectColor);
    end

    % Finalize plot
    set(ax, 'FontSize', fontSize, 'FontAngle', fontAngle);

    minTickIncrement = min([min(abs(diff(domainVisualizationTicks.x))) min(abs(diff(domainVisualizationTicks.y)))]);
    if (minTickIncrement >= 2)
       set(ax, 'XTickLabel', sprintf('%1.0f\n', domainVisualizationTicks.x), ...
               'YTickLabel', sprintf('%1.0f\n', domainVisualizationTicks.y));
    elseif (minTickIncrement >= 1)
       set(ax, 'XTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.x), ...
               'YTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.y));
    elseif (minTickIncrement >= 0.1)
       set(ax, 'XTickLabel', sprintf('%1.2f\n', domainVisualizationTicks.x), ...
               'YTickLabel', sprintf('%1.2f\n', domainVisualizationTicks.y));
    else
        set(ax, 'XTickLabel', sprintf('%1.3f\n', domainVisualizationTicks.x), ...
                'YTickLabel', sprintf('%1.3f\n', domainVisualizationTicks.y));
    end

    if (~identifyInputCones)
        colormap(ax, cMap);
    end


    if (isempty(backgroundColor))
        set(ax, 'CLim', [0 1], 'Color', 'none');
    else
        set(ax, 'CLim', [0 1], 'Color', backgroundColor);
    end

    if (isempty(colorbarTickLabelColor))
        colorbarTickLabelColor = [1 0.5 0];
    end
    
    box(ax, 'on')

    % Colorbar and colorbar ticks
    if (~isempty(activation))
        if (verticalColorBar) || (horizontalColorBar) || (verticalColorBarInside) || (horizontalColorBarInside)
            colorBarTicks = [0.00 0.25 0.5 0.75 1.0];
            colorBarTickLabels = cell(1, numel(colorBarTicks));
            colorBarTickLevels = activationRange(1) + (activationRange(2)-activationRange(1)) * colorBarTicks;
            
            for k = 1:numel(colorBarTicks)
                if (max(abs(colorBarTickLevels)) >= 10)
                    colorBarTickLabels{k} = sprintf('%2.0f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                elseif (max(abs(colorBarTickLevels)) >= 1)
                    colorBarTickLabels{k} = sprintf('%2.1f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                elseif (max(abs(colorBarTickLevels)) > 0.1)
                    colorBarTickLabels{k} = sprintf('%2.2f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                else
                    colorBarTickLabels{k} = sprintf('%2.3f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                end
            end
            
            if (isempty(colorbarFontSize))
                colorbarFontSize = fontSize/2;
            end
    
  
            if (verticalColorBar)
                colorbar(ax, 'eastOutside', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor);
            elseif (verticalColorBarInside)
                colorbar(ax, 'east', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor,  'FontWeight', 'Bold', 'FontSize', colorbarFontSize, 'FontName', 'Spot mono');
            elseif (horizontalColorBar)
                colorbar(ax,'northOutside', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor);
            elseif (horizontalColorBarInside)
                colorbar(ax,'north', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor,  'FontWeight', 'Bold', 'FontSize', colorbarFontSize, 'FontName', 'Spot mono');
            end
        else
            colorbar(ax, 'off');
        end
    else
        colorbar(ax, 'off');
    end

    if (plotTitle)    
       title(ax, plotTitle, 'Color', plotTitleColor, 'FontSize', plotTitleFontSize);
    end


    fprintf('\nDrawning mRGCmosaic patch. Please wait ...');
    tic
    drawnow;
    fprintf('Done in %2.2f seconds\n', toc);

end


function renderPooledConesLineSegments(obj,ax, lConeInputLineColor, mConeInputLineColor, lineSegmentWidth)
   
    if (~isempty(obj.visualizationCache.rfCenterConeConnectionLineSegments))
        % Plot the connections from the RF center to the input L-cones
        idx = find(obj.visualizationCache.rfCenterConeConnectionLineSegments.coneTypes == cMosaic.LCONE_ID);
    
        plot(ax, ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
            'Color', [0 0 0],...
            'LineWidth', lineSegmentWidth*2);  

        plot(ax, ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
            'Color', lConeInputLineColor,...
            'LineWidth', lineSegmentWidth); 

    
        idx = find(obj.visualizationCache.rfCenterConeConnectionLineSegments.coneTypes == cMosaic.MCONE_ID);
        plot(ax, ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
            'Color', [0 0 0],...
            'LineWidth', lineSegmentWidth*2);
        plot(ax, ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Xpos(:, idx), ...
            obj.visualizationCache.rfCenterConeConnectionLineSegments.Ypos(:,idx), ...
            'Color', mConeInputLineColor,...
            'LineWidth', lineSegmentWidth);
    end

end




