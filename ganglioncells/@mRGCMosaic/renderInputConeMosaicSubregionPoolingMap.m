function subregionLineWeightingFunctions = renderInputConeMosaicSubregionPoolingMap(theAxes, ...
            theInputConeMosaic, theContourData, ...
            theSubregionConeIndices, theSubregionConeWeights, flatTopSaturationLevel, ...
            spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
            scaleBarDegs, gridless, noXLabel, noYLabel, varargin)

     % Parse optional input
    p = inputParser;
    p.addParameter('computeLineWeightingFunctionsAndReturn', false, @islogical);
    p.addParameter('clearAxes', true, @islogical);
    p.addParameter('doNotLabelScaleBar', false, @islogical);
    p.addParameter('renderConeMap', true, @islogical);
    p.addParameter('spatialSupportSamplesNum', 300, @isscalar);
    p.addParameter('plotTitle', '', @ischar);
    p.parse(varargin{:});

    computeLineWeightingFunctionsAndReturn = p.Results.computeLineWeightingFunctionsAndReturn;
    clearAxes = p.Results.clearAxes;
    doNotLabelScaleBar = p.Results.doNotLabelScaleBar;
    renderConeMap = p.Results.renderConeMap;
    spatialSupportSamplesNum = p.Results.spatialSupportSamplesNum;
    plotTitle = p.Results.plotTitle;
    

    if (isempty(spatialSupportTickSeparationArcMin))
        spatialSupportTickSeparationArcMin = 3;
    end
    
    spatialSupportRangeArcMin = 4 * spatialSupportTickSeparationArcMin;
    maxXY = round(spatialSupportRangeArcMin/2);
    spatialSupportSamplingArcMin = spatialSupportRangeArcMin/spatialSupportSamplesNum;
    spatialSupportDegs = (-maxXY:spatialSupportSamplingArcMin:maxXY)/60;
    spatialSupportXYDegs(:,1) = spatialSupportCenterDegs(1) + spatialSupportDegs;
    spatialSupportXYDegs(:,2) = spatialSupportCenterDegs(2) + spatialSupportDegs;

    dx = (spatialSupportDegs(end)-spatialSupportDegs(1))*0.05;
    XLims = spatialSupportCenterDegs(1) + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];
    YLims = spatialSupportCenterDegs(2 )+ [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];


    % Generate image of pooled cones
    [subregionConeMap, subregionConeMapFlatTop] = mRGCMosaic.retinalSubregionConeMapFromPooledConeInputs(...
      theInputConeMosaic, theSubregionConeIndices, theSubregionConeWeights, spatialSupportXYDegs, ...
      flatTopSaturationLevel);
        
    % Return the subregion line weighting functions along X and Y
    subregionLineWeightingFunctions.xProfile = struct(...
        'spatialSupportDegs', spatialSupportXYDegs(:,1), ...
        'amplitude', sum(subregionConeMap,1));

    subregionLineWeightingFunctions.yProfile = struct(...
        'spatialSupportDegs', spatialSupportXYDegs(:,2), ...
        'amplitude', sum(subregionConeMap,2));


    if (computeLineWeightingFunctionsAndReturn)
        return;
    end

    if (clearAxes)
        cla(theAxes, 'reset');
    end
    if (renderConeMap)
        imagesc(theAxes, spatialSupportXYDegs(:,1), spatialSupportXYDegs(:,2), subregionConeMapFlatTop);
    end

    hold(theAxes, 'on');

    if (~isempty(theContourData))
        S.Vertices = theContourData.vertices;
        S.Faces = theContourData.faces;
        S.FaceVertexCData = [0.5 0.5 0.5];
        S.FaceColor = 'flat';
        S.EdgeColor = [0 0 0];
        S.FaceAlpha = 0.0;
        S.LineWidth = 2.0;
        patch(S, 'Parent', theAxes);
    end

    % The input cones
    for iInputCone = 1:numel(theSubregionConeIndices)
        switch theInputConeMosaic.coneTypes(theSubregionConeIndices(iInputCone))
            case cMosaic.LCONE_ID
                coneColor = [1 0 0];
            case cMosaic.MCONE_ID
                coneColor = [0 0.75 0];
            case cMosaic.SCONE_ID
                coneColor = [0 0 1];
        end
        plot(theAxes, ...
            theInputConeMosaic.coneRFpositionsDegs(theSubregionConeIndices(iInputCone),1), ...
            theInputConeMosaic.coneRFpositionsDegs(theSubregionConeIndices(iInputCone),2), ...
            'o', 'MarkerSize', 12, 'MarkerFaceColor', coneColor, 'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.0);
    end

    axis(theAxes, 'image'); axis(theAxes, 'xy');

    if (gridless)
        grid(theAxes, 'off');
    else
        grid(theAxes, 'on');
    end

    
    % Add a scale bar for comparison with physiology

    if ((~isempty(scaleBarDegs)) && (scaleBarDegs>0.00))
        xOffset = XLims(1)+0.05*(XLims(2)-XLims(1));
        yOffset = YLims(1)+0.05*(YLims(2)-YLims(1));
        plot(theAxes, xOffset+[0 scaleBarDegs], yOffset*[1 1], 'r-', 'LineWidth', 2);

        if (~doNotLabelScaleBar)
            if (scaleBarDegs>=1.0)
                text(theAxes, xOffset+0*scaleBarDegs, yOffset, sprintf(' %2.1f degs', scaleBarDegs), 'FontSize', 20, 'Color', [1 0 0]);
            elseif (scaleBarDegs>=0.1)
                text(theAxes, xOffset+0*scaleBarDegs, yOffset, sprintf(' %2.2f degs', scaleBarDegs), 'FontSize', 20, 'Color', [1 0 0]);
            else
                text(theAxes, xOffset+0*scaleBarDegs, yOffset, sprintf(' %2.3f degs', scaleBarDegs), 'FontSize', 20, 'Color', [1 0 0]);
            end
        end
    end

    xTicksDegs = spatialSupportCenterDegs(1) + (-3:3)*spatialSupportTickSeparationArcMin/60;
    yTicksDegs = spatialSupportCenterDegs(2) + (-3:3)*spatialSupportTickSeparationArcMin/60;

    if (spatialSupportTickSeparationArcMin/60 > 1-100*eps)
       xTickLabels = sprintf('%2.0f\n', xTicksDegs);
       yTickLabels = sprintf('%2.0f\n', yTicksDegs);
    elseif (spatialSupportTickSeparationArcMin/60>= 0.1-100*eps)
       xTickLabels = sprintf('%2.1f\n', xTicksDegs);
       yTickLabels = sprintf('%2.1f\n', yTicksDegs);
    elseif (spatialSupportTickSeparationArcMin/60 > 0.01-100*eps)
       xTickLabels = sprintf('%2.2f\n', xTicksDegs);
       yTickLabels = sprintf('%2.2f\n', yTicksDegs);
    else
       xTickLabels = sprintf('%2.3f\n', xTicksDegs);
       yTickLabels = sprintf('%2.3f\n', yTicksDegs);
    end

    set(theAxes, ...
        'CLim', [0 max(subregionConeMapFlatTop(:))], ...
        'XLim', XLims, 'YLim', YLims, ...
        'XTick', xTicksDegs, 'YTick', yTicksDegs, ...
        'XTickLabel', xTickLabels, ...
        'YTickLabel', yTickLabels);

    if (~noXLabel)
        xlabel(theAxes, 'eccentricity, x (degs)');
    end
    if (~noYLabel)
        ylabel(theAxes, 'eccentricity, y (deg)');
    end

    title(theAxes, plotTitle);

    xtickangle(theAxes, 0);
    colormap(theAxes, brewermap(1024, 'greys'));
    drawnow;
end
