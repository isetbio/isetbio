function generateDeconvolutionModel(obj, deconvolutionOpticsParams, modelPrefix)
    
    % Validate the deconvolutionOpticsParams
    obj.validateDeconvolutionOpticsParams(deconvolutionOpticsParams);
    quadrantsToAverage = deconvolutionOpticsParams.quadrantsToAverage;
    subjectsToAverage = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToAverage;
    
    eccTested = [0 0.25 0.5 1 1.5 2.0 2.5 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
    eccTested = [0 0.25 0.5 1 1.5 2.0 2.5 3 4];
    rootDir = obj.psfDeconvolutionDir;
  
    defocusMode = 'subjectDefault';
    if (strcmp(defocusMode, 'subjectDefault'))
        imposedRefractionErrorDiopters = 0; 
    else
        imposedRefractionErrorDiopters =  0.01; 
    end
    
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
        
    w = WatsonRGCModel();
    coneRFSpacingsDegs  = w.coneRFSpacingAndDensityAlongMeridian(abs(eccTested), ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    coneApertureRadii = 0.8*0.5*coneRFSpacingsDegs;
    
    for eccIndex = 1:numel(eccTested)
        eccDegs = eccTested(eccIndex);
        dataFileName = fullfile(rootDir, sprintf('ecc_-%2.1f_deconvolutions_refractionError_%2.2fD.mat', eccDegs, imposedRefractionErrorDiopters));
        load(dataFileName, 'retinalPoolingRadii', 'visualRadius', 'subjectIDs', 'quadrants');
        
        if (eccDegs == 0)
            w = WatsonRGCModel();
            coneRFSpacingDeg  = w.coneRFSpacingAndDensityAlongMeridian(0, ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio',false);
            % min retinal radius = 0.5 * cone spacing at 0 deg  
            minRetinalRadius = 0.8*0.5*coneRFSpacingDeg;
            % Find closest retinal pooling radius
            [dd,idx] = sort(abs(minRetinalRadius-retinalPoolingRadii));
            minRetinalRadius = retinalPoolingRadii(idx(1));
            d1 = dd(2)/(dd(1)+dd(2));
            d2 = dd(1)/(dd(1)+dd(2));
            visualRadii = d1 * squeeze(visualRadius(idx(1),:,:)) + d2 * squeeze(visualRadius(idx(2),:,:));
            minVisualRadiusDegs = median(median(visualRadii));
        end
        
        % Get data for the quadrant of interest
        visualRadius = CronerKaplanRGCModel.quadrantData(visualRadius, quadrantsToAverage, quadrants, subjectsToAverage, subjectIDs);
        
        % Only include points for retinal pooling radii > = cone aperture
        idx = find(retinalPoolingRadii >= coneApertureRadii(eccIndex));
        retinalPoolingRadii = retinalPoolingRadii(1,idx);
        visualRadius = visualRadius(:,idx);
 
        row = floor((eccIndex-1)/plotCols)+1;
        col = mod(eccIndex-1, plotCols)+1;

        % median over included subjects/ecc quadrants
        medianVisualRadius = (median(visualRadius,1, 'omitnan'))';
        
        
        % Extend the visual radius data to 1.5
        medianVisualRadiusExtended = [medianVisualRadius; [0.7 0.8 1 1.5]'];
        retinalPoolingRadiiExtended = [retinalPoolingRadii [0.7 0.8 1 1.5]]';
        
        % Fit the visualRadius with a saturating function
        modelFunctionRadius = @(p,x)(p(1) - (p(2)-p(1))*exp(-p(3)*x));
        initialParams = [5 10 0.2];
        [fittedParams, fittedParamsSE] = nonLinearFitData(medianVisualRadiusExtended, retinalPoolingRadiiExtended, modelFunctionRadius, initialParams);
        fittedParamsRadius(eccIndex,:) = fittedParams;
        medianVisualRadiusFit = logspace(log10(0.01), log10(1), 100);
        retinalPoolingRadiiFit = modelFunctionRadius(fittedParams, medianVisualRadiusFit);
      
        if (row <=  plotRows) && (col <= plotCols)
            yLims = [0.003 0.6];
            showIdentityLine = true;
            plotData(theAxesGrid{row,col}, eccDegs, retinalPoolingRadii, visualRadius,  medianVisualRadius, ...
                retinalPoolingRadiiFit, medianVisualRadiusFit, subjectIDs, quadrants, ...
                row==plotRows, col==1, imposedRefractionErrorDiopters, yLims, showIdentityLine, coneApertureRadii(eccIndex), 'visual pooling radius');
        end
    end
           
    if (strcmp(defocusMode, 'subjectDefault'))
        pdfFileName = sprintf('%s_defocusDefault_radius_', modelPrefix);
    else
        pdfFileName = sprintf('%s_defocus%2.2fD_radius_', modelPrefix, imposedRefractionErrorDiopters);
    end
    plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, pwd());
    
    
    % Display the visual gain data
    hFig = figure(101); clf;
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
        dataFileName = fullfile(rootDir, sprintf('ecc_-%2.1f_deconvolutions_refractionError_%2.2fD.mat', eccDegs, imposedRefractionErrorDiopters));
        load(dataFileName, 'retinalPoolingRadii',  'visualGain', 'subjectIDs', 'quadrants');
        retinalPoolingRadiiOriginal = retinalPoolingRadii;
        
        % Get data for the quadrant of interest
        visualGain = CronerKaplanRGCModel.quadrantData(visualGain, quadrantsToAverage, quadrants, subjectsToAverage, subjectIDs);
        
        % Only include points for retinal pooling radii > = cone aperture
        idx = find(retinalPoolingRadii >= 0.5*coneApertureRadii(eccIndex));
        retinalPoolingRadii = retinalPoolingRadii(1,idx);
        visualGain = visualGain(:,idx);
        
        row = floor((eccIndex-1)/plotCols)+1;
        col = mod(eccIndex-1, plotCols)+1;
        
        % Median over all subjects/ecc quadrants
        medianVisualGain = (median(visualGain,1, 'omitnan'))';
        
        % Extend the visual gain data 
        medianVisualGainExtended = [medianVisualGain; [1 1 1]'];
        retinalPoolingRadiiExtended = [retinalPoolingRadii [0.8 1 2]]';
        
        % Fit the visualGain with a saturating function
        modelFunctionGain = @(p,x)((p(1)*x.^p(3))./(x.^p(3)+p(2)));
        initialParams = [0.9 0.002 2];
        [fittedParams, fittedParamsSE] = nonLinearFitData(retinalPoolingRadiiExtended, medianVisualGainExtended, modelFunctionGain, initialParams);
        fittedParamsGain(eccIndex,:) = fittedParams;
        retinalPoolingRadiiFit = logspace(log10(0.001), log10(1), 100);
        medianVisualGainFit = modelFunctionGain(fittedParams, retinalPoolingRadiiFit);
        
        if (row <=  plotRows) && (col <= plotCols)
            yLims = [0.03 1];
            showIdentityLine = false;
            plotData(theAxesGrid{row,col}, eccDegs, retinalPoolingRadii, visualGain,  medianVisualGain, ...
                retinalPoolingRadiiFit, medianVisualGainFit, subjectIDs, quadrants, ...
                row==plotRows, col==1, imposedRefractionErrorDiopters, yLims, showIdentityLine, ...
                coneApertureRadii(eccIndex), 'peak sensitivity attenuation');
        end
    end
    
    if (strcmp(defocusMode, 'subjectDefault'))
        pdfFileName = sprintf('%s_defocusDefault_gain_', modelPrefix);
    else
        pdfFileName = sprintf('%s_defocus%2.2fD_gain_', modelPrefix, imposedRefractionErrorDiopters);
    end
            
    plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, pwd());
    
    matFileName = strrep(pdfFileName, '_gain', '');
    retinalPoolingRadii = retinalPoolingRadiiOriginal;
    modelFileName = fullfile(rootDir, sprintf('%sFittedModel.mat', matFileName));
    save(modelFileName, 'retinalPoolingRadii', ...
        'eccTested', 'fittedParamsRadius', 'fittedParamsGain', ...
        'modelFunctionRadius', 'modelFunctionGain', 'minVisualRadiusDegs');
   
    visualizeFittedModel(modelFileName);
end

function visualizeFittedModel(modelFileName)

    load(modelFileName, 'eccTested', 'fittedParamsRadius', 'fittedParamsGain', 'modelFunctionRadius', 'modelFunctionGain');
   
    cMap =  (brewermap(numel(eccTested), 'spectral'))/1.3;
    plotlabOBJ = setupPlotLabForFittedModel([15 8], 'both', cMap);
    
    hFig = figure(1234); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 1, ...
            'colsNum', 2, ...
            'leftMargin', 0.08, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.1, ...
            'heightMargin', 0.08, ...
            'bottomMargin', 0.1, ...
            'topMargin', 0.06);
        
    w = WatsonRGCModel();
    coneRFSpacingsDegs  = w.coneRFSpacingAndDensityAlongMeridian(abs(eccTested), ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', true);
    coneApertureRadii = 0.8*0.5*coneRFSpacingsDegs;
        
    ax = theAxesGrid{1,1};
    hold(ax, 'on');
    visualRadiusFit = logspace(log10(0.01), log10(1), 64);
    for eccIndex = 1:numel(eccTested)
        retinalPoolingRadius = modelFunctionRadius(fittedParamsRadius(eccIndex,:), visualRadiusFit);
        idx = find(retinalPoolingRadius >= 0.5*coneApertureRadii(eccIndex));
        
        line(ax, [0.001 1], [0.001 1], 'LineStyle','--');
        line(ax,retinalPoolingRadius(idx), visualRadiusFit(idx), 'Color', cMap(eccIndex,:));
        scatter(ax,retinalPoolingRadius(idx), visualRadiusFit(idx),49);
        
    end
    axis(ax, 'square');
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.003, 1], 'YLim', [0.003 1], ...
        'XTick', [0.01 0.03 0.1 0.3 1], 'YTick', [0.01 0.03 0.1 0.3 1]);
    xlabel(ax, 'visual radius (degs)');
    ylabel(ax, 'retinal radius (degs)');
    
    ax = theAxesGrid{1,2};
    hold(ax, 'on');
    legends = {};
    retinalPoolingRadiusFit = logspace(log10(0.003), log10(1), 64);
    for eccIndex = 1:numel(eccTested)
        idx = find(retinalPoolingRadiusFit >= 0.5*coneApertureRadii(eccIndex));
        retinalPoolingRadiusFitForThisEcc = retinalPoolingRadiusFit(idx);
        gainAttenuation = modelFunctionGain(fittedParamsGain(eccIndex,:), retinalPoolingRadiusFitForThisEcc);
        line(ax,retinalPoolingRadiusFitForThisEcc, gainAttenuation, 'Color', cMap(eccIndex,:));
        legends{numel(legends)+1} = sprintf('ecc = %2.1f^o\n', eccTested(eccIndex));
    end
    for eccIndex = 1:numel(eccTested)
         idx = find(retinalPoolingRadiusFit >= 0.5*coneApertureRadii(eccIndex));
        retinalPoolingRadiusFitForThisEcc = retinalPoolingRadiusFit(idx);
        gainAttenuation = modelFunctionGain(fittedParamsGain(eccIndex,:), retinalPoolingRadiusFitForThisEcc);
        scatter(ax,retinalPoolingRadiusFitForThisEcc, gainAttenuation, 49);
    end
    axis(ax, 'square');
    %legend(ax, legends, 'Location', 'SouthEast');
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.003, 1], 'YLim', [0 1], ...
        'XTick',  [0.001 0.003 0.01 0.03 0.1 0.3 1], 'YTick', [0.01 0.03 0.1 0.3 1]);
    xlabel(ax, 'retinal radius (degs)');
    ylabel(ax, 'peak sensitivity attenuation');
    
    plotlabOBJ.exportFig(hFig, 'pdf', 'modelSummary', pwd());
end


function [fittedParams, fittedParamsSE] = nonLinearFitData(x,y, modelFunction, initialParams)
    opts.RobustWgtFun = 'talwar';
    [fittedParams,~,~,varCovarianceMatrix,~] = nlinfit(x,y,modelFunction,initialParams,opts);
    % standard error of the mean
    fittedParamsSE = sqrt(diag(varCovarianceMatrix));
    fittedParamsSE = fittedParamsSE';
end

function renderPSF(ax, thePSFSupportDegs, thePSF, xLims, cMap, theTitle, showXLabel, peakV)
    imagesc(ax, thePSFSupportDegs, thePSFSupportDegs, thePSF); hold(ax, 'on');
    zLevels = logspace(log10(0.05), log10(1), 6)*peakV;
    [X,Y] = meshgrid(thePSFSupportDegs, thePSFSupportDegs);
    contour(ax, X,Y, thePSF, zLevels, 'LineColor', [0.5 0.5 0.5]);
    midRow = round(size(thePSF,1)/2)+1;
    profile = thePSF(midRow,:);
    deltaX = xLims(2)-xLims(1);
    area(ax,thePSFSupportDegs, xLims(1) + deltaX*profile/peakV, xLims(1), 'FaceColor', [1 0.5 0.5], 'EdgeColor', [1 0 0]);
    %plot(ax, [0 0], [-1 1], 'k-'); plot(ax, [-1 1], [0 0], 'k-');
    axis(ax, 'square');  axis(ax, 'xy');
    set(ax, 'XLim', xLims, 'YLim', xLims, 'CLim', [0 peakV]);
    set(ax, 'XTickLabel', {}, 'YTickLabel', {});
    set(ax, 'XTick', -1:0.05:1, 'YTick', -1:0.05:1, 'XTickLabel', {}, 'YTickLabel', {});
    if (showXLabel)
        set(ax, 'XTickLabel', [-1:0.05:1]*60);
        xlabel(ax,'arc min');
    end
    colormap(ax,cMap);
    title(ax, theTitle);
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
 