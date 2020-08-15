function plotDeconvolutionModel(deconvolutionModel)
    
    visualizeCenterDeconvolutionModel(deconvolutionModel.center);
    visualizeSurroundDeconvolutionModel(deconvolutionModel.surround);
    
end

function visualizeSurroundDeconvolutionModel(deconvolutionModel)
    hFig = figure(2);
    set(hFig, 'Name', sprintf('RF surround deconvolution model for subject %d and ''%s'' quadrant.', ...
        deconvolutionModel.subjectID, deconvolutionModel.quadrant));
    
    rowsNum = 2;
    colsNum = size(deconvolutionModel.eccDegs,1);
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.04, ...
            'bottomMargin', 0.04, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.01, ...
            'heightMargin', 0.08, ...
            'topMargin', 0.03, ...
            'rowsNum', rowsNum, ...
            'colsNum', colsNum);
        
end

function visualizeCenterDeconvolutionModel(deconvolutionModel)

    w = WatsonRGCModel();
    retinalNasalEccentricitiesDegs = logspace(log10(0.1), log10(30), 100);
    coneRFSpacingsDegs = w.coneRFSpacingAndDensityAlongMeridian(retinalNasalEccentricitiesDegs, ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    coneApertureDiameterDegs = WatsonRGCModel.coneApertureToDiameterRatio * coneRFSpacingsDegs;
    coneRadiiDegs = 0.5*coneApertureDiameterDegs;
    
    ck = CronerKaplanRGCModel('generateAllFigures', false);
    ck.setupPlotLab(0, 25, 15);
    
    hFig = figure(1);
    set(hFig, 'Name', sprintf('RF center deconvolution model for subject %d and ''%s'' quadrant.', ...
        deconvolutionModel.subjectID, deconvolutionModel.quadrant));
    
    rowsNum = 2;
    colsNum = size(deconvolutionModel.eccDegs,1);
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.04, ...
            'bottomMargin', 0.04, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.01, ...
            'heightMargin', 0.08, ...
            'topMargin', 0.03, ...
            'rowsNum', rowsNum, ...
            'colsNum', colsNum);
        
    
   
    for eccIndex = 1:size(deconvolutionModel.eccDegs,1)
        [~,idx]  = min(abs(retinalNasalEccentricitiesDegs-deconvolutionModel.eccDegs(eccIndex,1)));
        coneRadiusDegsAtClosestEccentricity = coneRadiiDegs(idx);
        
        % Visual characteristic radius as a function of number of cone in RF center
        ax = theAxesGrid{1, eccIndex};
        
        scatter(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), deconvolutionModel.visualCharacteristicRadiusMin(eccIndex,:), ...
            'r^', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.8 0.3 0.3]);
        hold(ax,'on');
        
        scatter(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), deconvolutionModel.visualCharacteristicRadiusMax(eccIndex,:), ...
            'rv', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.8 0.3 0.3]);
        
        % The mean (min+max)/2 characteristic radius of the VISUAL center
        plot(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), ...
              0.5*(deconvolutionModel.visualCharacteristicRadiusMax(eccIndex,:)+deconvolutionModel.visualCharacteristicRadiusMin(eccIndex,:)), ...
            'r-', 'LineWidth', 2);
        
        % The mean (min+max)/2 characteristic radius of the RETINAL center
        plot(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), ...
            0.5*(deconvolutionModel.retinalCharacteristicRadiusMin(eccIndex,:) +  deconvolutionModel.retinalCharacteristicRadiusMax(eccIndex,:)), ...
            'b--', 'LineWidth', 1.5);
        
        % The characteristic radius of a single cone
        plot(ax,retinalNasalEccentricitiesDegs, retinalNasalEccentricitiesDegs*0 + coneRadiusDegsAtClosestEccentricity, ...
            'k--', 'LineWidth', 1.5);
         
        set(ax, 'XLim', [0.9 31], 'YLim', [0.001 0.1], ...
            'YTick', [0.003 0.01 0.03 0.1 1], 'YTickLabel', {'0.003', '0.01', '0.03', '0.1', '0.3', '1.0'}, ...
            'XTick', [1 2 3 4 6 10 15 30 100], ...
            'XScale', 'log', 'YScale', 'log');
        if (eccIndex == 1)
            ylabel(ax, ['characteristic radius ({\color{red} visual, \color{blue} retinal})']);
        else
            set(ax, 'YTickLabel', {})
        end
        title(ax,sprintf('ecc: (%2.1f, %2.1f)', deconvolutionModel.eccDegs(eccIndex,1), deconvolutionModel.eccDegs(eccIndex,2)));
          
        
        % Visual gain attenuation as a function of number of cone in RF center
        ax = theAxesGrid{2, eccIndex};
        plot(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), deconvolutionModel.visualGain(eccIndex,:), ...
            'r-', 'MarkerFaceColor', [1 0.5 0.5]);
        hold(ax, 'on');
        plot(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), deconvolutionModel.retinalGain(eccIndex,:), ...
            'b-', 'LineWidth', 1.5);
        set(ax, 'XLim', [0.9 31], 'YLim', [0 1.04], 'XScale', 'log', ...
            'XTick', [1 2 3 4 6 10 15 30 100])
        xlabel(ax,'# of center cones');
        if (eccIndex == 1)
            ylabel(ax, ['peak gain ({\color{red} visual, \color{blue} retinal})']);
        else
            set(ax, 'YTickLabel', {})
        end
        drawnow
   end
               
end







function old()
    eccTested = deconvolutionModel.tabulatedEccentricities;
    subjectIDs = deconvolutionModel.opticsParams.PolansWavefrontAberrationSubjectIDsToAverage;
    quadrants = deconvolutionModel.opticsParams.quadrantsToAverage;
    medianVisualRadiusFit = logspace(log10(0.01), log10(1), 100);
    
    plotlabOBJ = setupPlotLab([26 14], 'both');
    
    % Display the visual radius data
    hFig = figure(100); clf;
    plotRows = 4; plotCols = 7;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', plotRows, ...
            'colsNum', plotCols, ...
            'leftMargin', 0.04, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.01, ...
            'heightMargin', 0.02, ...
            'bottomMargin', 0.05, ...
            'topMargin', 0.01);
        
    for eccIndex = 1:numel(eccTested)
        eccDegs = eccTested(eccIndex);
        row = floor((eccIndex-1)/plotCols)+1;
        col = mod(eccIndex-1, plotCols)+1;
        
        if (row <=  plotRows) && (col <= plotCols)
            yLims = [0.003 0.6];
            showIdentityLine = true;
            visualRadius = deconvolutionModel.nonAveragedVisualRadius{eccIndex};
            medianVisualRadius = (median(visualRadius,1, 'omitnan'))';
            retinalPoolingRadiiFit = deconvolutionModel.modelFunctionRadius(deconvolutionModel.fittedParamsRadius(eccIndex,:), medianVisualRadiusFit);
            idx = find(deconvolutionModel.retinalPoolingRadii >= deconvolutionModel.coneApertureRadii(eccIndex));
            retinalPoolingRadii = deconvolutionModel.retinalPoolingRadii(1,idx);
        
            plotData(theAxesGrid{row,col}, eccDegs, retinalPoolingRadii, visualRadius,  medianVisualRadius, ...
                retinalPoolingRadiiFit, medianVisualRadiusFit, subjectIDs, quadrants, ...
                row==plotRows, col==1, 0, yLims, showIdentityLine, deconvolutionModel.coneApertureRadii(eccIndex), 'visual pooling radius');
        end
    end
    
    
    % Display the visual gain data
    hFig = figure(101); clf;
    plotRows = 4; plotCols = 7;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', plotRows, ...
            'colsNum', plotCols, ...
            'leftMargin', 0.04, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.01, ...
            'heightMargin', 0.02, ...
            'bottomMargin', 0.05, ...
            'topMargin', 0.01);
        
    for eccIndex = 1:numel(eccTested)
        eccDegs = eccTested(eccIndex);
        row = floor((eccIndex-1)/plotCols)+1;
        col = mod(eccIndex-1, plotCols)+1;
        
        if (row <=  plotRows) && (col <= plotCols)
            yLims = [0.03 1];
            showIdentityLine = ~true;
            visualGain = deconvolutionModel.nonAveragedVisualGain{eccIndex};
            medianVisualGain = (median(visualGain,1, 'omitnan'))';
            retinalPoolingRadiiFit = logspace(log10(0.001), log10(1), 100);
            medianVisualGainFit = deconvolutionModel.modelFunctionGain(deconvolutionModel.fittedParamsGain(eccIndex,:), retinalPoolingRadiiFit);
            idx = find(deconvolutionModel.retinalPoolingRadii >= deconvolutionModel.coneApertureRadii(eccIndex));
            retinalPoolingRadii = deconvolutionModel.retinalPoolingRadii(1,idx);
            
            plotData(theAxesGrid{row,col}, eccDegs, retinalPoolingRadii, visualGain,  medianVisualGain, ...
                retinalPoolingRadiiFit, medianVisualGainFit, subjectIDs, quadrants, ...
                row==plotRows, col==1, 0, yLims, showIdentityLine, deconvolutionModel.coneApertureRadii(eccIndex), 'peak sensitivity attenuation');
        end
    end
    
        
end

function plotData(theAxes, ecc, xQuantity, visualQuantity, visualQuantityMedian,  ...
    fittedXQuantity, fittedVisualQuantityMedian, subjectIDs, quadrants, ...
    showXLabel, showYLabel, imposedRefractionErrorDiopters, yLims, showIdentityLine, identityLineBreakPoint, theXLabel)
    
    hold(theAxes, 'on');
    if (showIdentityLine)
        line(theAxes,  [identityLineBreakPoint 1], identityLineBreakPoint*[1 1], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.0);
        line(theAxes, identityLineBreakPoint*[1 1], [identityLineBreakPoint 1], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.0);
        line(theAxes, [identityLineBreakPoint 1], [identityLineBreakPoint 1], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.0);
    else
        line(theAxes, identityLineBreakPoint*[1 1], [0.003 1], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.0);
    end

    % Individual subjects/quadrants
    for k = 1:size(visualQuantity,1)
        line(theAxes, xQuantity, visualQuantity(k,:), ...
                'Color', 0.8*[0.6 0.6 1], 'LineWidth', 1.5); 
    end

 
    plot(theAxes, fittedXQuantity, fittedVisualQuantityMedian, 'r-', 'LineWidth',2.0);
    plot(theAxes, xQuantity, visualQuantityMedian, 'o', ...
        'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgecolor', [1 0 0]);
    
    grid(theAxes, 'off');
    axis(theAxes, 'square');
    set(theAxes, 'XScale', 'log', 'YScale', 'log');
    set(theAxes, 'XTick', [0.01 0.03 0.1 0.3 0.6], 'YTick', [0.01 0.03 0.1 0.3 0.6 1]);
    set(theAxes, 'XLim', [0.003 0.6], 'YLim', yLims);
    
    if (imposedRefractionErrorDiopters == 0)
        title(theAxes, sprintf('%2.1f degs', ecc));
    else
        title(theAxes, sprintf('%2.1f degs)', ecc));
    end
    
    if (showXLabel)
        xlabel(theAxes, 'retinal pooling radius');
    else
        set(theAxes, 'XTickLabel', {});
    end
    if (showYLabel)
        ylabel(theAxes, theXLabel);
    else
        set(theAxes, 'YTickLabel', {});
    end

end

function plotlabOBJ = setupPlotLab(figSize, tickDir)
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', [0 0 0], ...
            'lineColor', [0.5 0.5 0.5], ...
            'scatterMarkerEdgeColor', [0.3 0.3 0.9], ...
            'lineMarkerSize', 12, ...
            'axesBox', 'off', ...
            'axesTickDir', tickDir, ...
            'renderer', 'painters', ...
            'axesTickLength', [0.02 0.02], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', figSize(1), ...
            'figureHeightInches', figSize(2));
end

function plotlabOBJ = setupPlotLabForFittedModel(figSize, tickDir, colorOrder)
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', colorOrder, ...
            'lineColor', [0.5 0.5 0.5], ...
            'scatterMarkerEdgeColor', [0.3 0.3 0.9], ...
            'lineMarkerSize', 11, ...
            'axesBox', 'off', ...
            'axesTickDir', tickDir, ...
            'renderer', 'painters', ...
            'axesTickLength', [0.02 0.02], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', figSize(1), ...
            'figureHeightInches', figSize(2));
end
 