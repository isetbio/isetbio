function analyzeConvolutionResults
    
    defocusMode = 'subjectDefault';
    goodSubjects = [1 3 5 6 7 8 9 10] %[5 6 7 8 9 10]; - 7 -9
    
    groupsToAnalyze = [1 2 3];
    goodSubjectIndices = [];
    for k = 1:numel(groupsToAnalyze)
        goodSubjectIndices = cat(2, goodSubjectIndices, goodSubjects + (k-1)*10);
    end

    
    if (strcmp(defocusMode, 'subjectDefault'))
        imposedRefractionErrorDiopters = 0; 
    else
        imposedRefractionErrorDiopters =  0.01; 
    end
    
    eccTested = [0 0.5 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
    plotRows = 3;
    plotCols = 7;
    
    
    % Display the employed PSFs
    showPSFs = ~true;
    if (showPSFs)
        plotThePSFs(goodSubjects, eccTested, [3 1 2]);
    end
 

    plotlabOBJ = setupPlotLab([26 14], 'both');
    
    % Display the visual radius data
    hFig = figure(100); clf;
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
    coneRFSpacingsDegs  = w.coneRFSpacingAndDensityAlongMeridian(eccTested, ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    coneApertureRadii = 0.5*coneRFSpacingsDegs;
    
    for eccIndex = 1:numel(eccTested)
        eccDegs = eccTested(eccIndex);
        dataFileName = sprintf('cellEcc_%2.1f_cellRefractionError_%2.2fD_VisualGain.mat', -eccDegs, imposedRefractionErrorDiopters);
        load(dataFileName, 'retinalPoolingRadii', 'visualRadius');
      
        % Only include points for retinal pooling radii > = cone aperture
        idx = find(retinalPoolingRadii >= 0.5*coneApertureRadii(eccIndex));
        retinalPoolingRadii = retinalPoolingRadii(idx);
        visualRadius = visualRadius(idx, :);
        visualRadius = visualRadius(:, goodSubjectIndices);

        row = floor((eccIndex-1)/plotCols)+1;
        col = mod(eccIndex-1, plotCols)+1;

    
        medianVisualRadius = (median(visualRadius,2, 'omitnan'))';

        % Fit the visualRadius with a saturating function
        modelFunctionRadius = @(p,x)(p(1) - (p(2)-p(1))*exp(-p(3)*x));
        initialParams = [5 10 0.2];
        
        [fittedParams, fittedParamsSE] = nonLinearFitData([medianVisualRadius 0.7 0.8 1 1.5], [retinalPoolingRadii 0.7 0.8 1 1.5], modelFunctionRadius, initialParams);
        fittedParamsRadius(eccIndex,:) = fittedParams;
        medianVisualRadiusFit = logspace(log10(0.01), log10(1), 100);
        retinalPoolingRadiiFit = modelFunctionRadius(fittedParams, medianVisualRadiusFit);
      
        if (row <=  plotRows) && (col <= plotCols)
            yLims = [0.003 0.6];
            showIdentityLine = true;
            plotData(theAxesGrid{row,col}, eccDegs, retinalPoolingRadii, visualRadius,  medianVisualRadius, ...
                retinalPoolingRadiiFit, medianVisualRadiusFit, ...
                row==plotRows, col==1, imposedRefractionErrorDiopters, yLims, showIdentityLine, coneApertureRadii(eccIndex), 'visual pooling radius');

        end
    end
    
    if (strcmp(defocusMode, 'subjectDefault'))
        pdfFileName = sprintf('defocusDefault_radius');
    else
        pdfFileName = sprintf('defocus%2.2fD_radius', imposedRefractionErrorDiopters);
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
        dataFileName = sprintf('cellEcc_%2.1f_cellRefractionError_%2.2fD_VisualGain.mat', -eccDegs, imposedRefractionErrorDiopters);
        load(dataFileName, 'retinalPoolingRadii',  'visualGain');
        retinalPoolingRadiiOriginal = retinalPoolingRadii;
        
        % Only include points for retinal pooling radii > = cone aperture
        idx = find(retinalPoolingRadii >= 0.5*coneApertureRadii(eccIndex));
        retinalPoolingRadii = retinalPoolingRadii(idx);
        visualGain = visualGain(idx, :);
        visualGain = visualGain(:, goodSubjectIndices);
        
        
        row = floor((eccIndex-1)/plotCols)+1;
        col = mod(eccIndex-1, plotCols)+1;
        
        medianVisualGain = (median(visualGain,2, 'omitnan'))';
        
        % Fit the visualGain with a saturating function
        modelFunctionGain = @(p,x)(p(1) - (p(2)-p(1))*exp(-p(3)*x));
        modelFunctionGain = @(p,x)((p(1)*x.^p(3))./(x.^p(3)+p(2)));
        
        initialParams = [0.9 0.002 2];
        [fittedParams, fittedParamsSE] = nonLinearFitData([retinalPoolingRadii 0.8 1 2], [medianVisualGain 1 1 1], modelFunctionGain, initialParams);
        fittedParamsGain(eccIndex,:) = fittedParams;
        retinalPoolingRadiiFit = logspace(log10(0.001), log10(1), 100);
        medianVisualGainFit = modelFunctionGain(fittedParams, retinalPoolingRadiiFit);
        
        if (row <=  plotRows) && (col <= plotCols)
            yLims = [0.03 1];
            showIdentityLine = false;
            plotData(theAxesGrid{row,col}, eccDegs, retinalPoolingRadii, visualGain,  medianVisualGain, ...
                retinalPoolingRadiiFit, medianVisualGainFit, ...
                row==plotRows, col==1, imposedRefractionErrorDiopters, yLims, showIdentityLine, coneApertureRadii(eccIndex), 'peak sensitivity attenuation');
        end
        
    end
    
    if (strcmp(defocusMode, 'subjectDefault'))
        pdfFileName = sprintf('defocusDefault_gain');
    else
        pdfFileName = sprintf('defocus%2.2fD_gain', imposedRefractionErrorDiopters);
    end
            
    plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, pwd());
    
    matFileName = strrep(pdfFileName, '_gain', '');
    retinalPoolingRadii = retinalPoolingRadiiOriginal;
    save(sprintf('%sFittedModel.mat', matFileName), 'retinalPoolingRadii', 'eccTested', 'fittedParamsRadius', 'fittedParamsGain', 'modelFunctionRadius', 'modelFunctionGain');
   
    pause
    visualizeFittedModel(defocusMode);
    
end

function visualizeFittedModel(defocusMode)

    if (strcmp(defocusMode, 'subjectDefault'))
        matFileName = 'defocusDefaultFittedModel.mat';
    else
        matFileName = 'defocusNoneFittedModel.mat';
    end

    load(matFileName, 'eccTested', 'fittedParamsRadius', 'fittedParamsGain', 'modelFunctionRadius', 'modelFunctionGain');
    

   
    cMap =  (brewermap(numel(eccTested), 'spectral'))/1.3;
    plotlabOBJ = setupPlotLabForFittedModel([15 8], 'both', cMap);
    
    hFig = figure(1); clf;
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
    coneRFSpacingsDegs  = w.coneRFSpacingAndDensityAlongMeridian(eccTested, ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', ~false);
    coneApertureRadii = 0.5*coneRFSpacingsDegs;
        
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
    fittedXQuantity, fittedVisualQuantityMedian, ...
    showXLabel, showYLabel, imposedRefractionErrorDiopters, yLims, showIdentityLine, identityLineBreakPoint, theXLabel)
    
    hold(theAxes, 'on');
    if (showIdentityLine)
        line(theAxes,  [identityLineBreakPoint 1], identityLineBreakPoint*[1 1], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.0);
        line(theAxes, identityLineBreakPoint*[1 1], [identityLineBreakPoint 1], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.0);
        line(theAxes, [identityLineBreakPoint 1], [identityLineBreakPoint 1], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.0);
    else
        line(theAxes, identityLineBreakPoint*[1 1], [0.003 1], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.0);
    end

    usedSubjectsNum = size(visualQuantity,2)/3;
    groupColors = [0.4 0.4 0.4; 1 0.4 0.4; 0.4 0.4 1];
    for group = 3:-1:1
        color = groupColors(group,:);
        for subIndex = 1:usedSubjectsNum
            line(theAxes, xQuantity, visualQuantity(:,(group-1)*usedSubjectsNum+subIndex), ...
                'Color', 0.8*[0.6 0.6 1], 'LineWidth', 1.5); 
        end
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

function   plotThePSFs(goodSubjects, eccTested, groupIndices)

    plotlabOBJ = setupPlotLab([28 14], 'in');
    cMap = brewermap(512, 'greys');
    
    for subk = 1:numel(goodSubjects)
        
        subIndex = goodSubjects(subk);
        hFig = figure(subIndex); clf;
        theAxesGrid = plotlab.axesGrid(hFig, ...
                'rowsNum', 3, ...
                'colsNum', numel(eccTested), ...
                'leftMargin', 0.08, ...
                'rightMargin', 0.01, ...
                'widthMargin', 0.03, ...
                'heightMargin', 0.08, ...
                'bottomMargin', 0.06, ...
                'topMargin', 0.06);
            
       
        eccTested = sort(eccTested);
        
        for kkk = 1:numel(groupIndices)
            groupIndex = groupIndices(kkk);
            for eccIndex = 1:numel(eccTested)
                eccDegs = eccTested(eccIndex);
                switch (groupIndex)
                case 1
                    eccDegs = [eccDegs 0];
                case 2
                    eccDegs = [0 eccDegs];
                case 3
                    eccDegs = [0 -eccDegs];
                end
        
                eccXrange = eccDegs(1)*[1 1];
                eccYrange = eccDegs(2)*[1 1];
                deltaEcc = 0.5;
                % Compute the subject PSF at the desired eccentricity
                [hEcc, vEcc, thePSFs, thePSFsupportDegs] = CronerKaplanRGCModel.psfAtEccentricity(subIndex, ...
                    imposedRefractionErrorDiopters, eccXrange, eccYrange, deltaEcc);
        
                thePSF(groupIndex,eccIndex,:,:) = squeeze(thePSFs(1, 1, 1,:,:));
            end
        end
        peakV = max(thePSF(:))*0.9;
        
        for kkk = 1:numel(groupIndices)
            groupIndex = groupIndices(kkk);
            for eccIndex = 1:numel(eccTested)
                
                row = kkk;
                col = eccIndex;
                eccDegs = eccTested(eccIndex);
                switch (groupIndex)
                case 1
                    eccDegs = [eccDegs 0];
                case 2
                    eccDegs = [0 eccDegs];
                case 3
                    eccDegs = [0 -eccDegs];
                end
       
                xLims = [-0.1 0.1];
                
                theTitle = sprintf('ecc (degs) = (%2.1f , %2.1f)', eccDegs(1), eccDegs(2));
                renderPSF(theAxesGrid{row,col}, thePSFsupportDegs, squeeze(thePSF(groupIndex,eccIndex,:,:)), xLims, cMap, theTitle, (kkk == 3) && (col== 1), peakV);
                drawnow;
            end
            
            if (strcmp(defocusMode, 'subjectDefault'))
                pdfFileName = sprintf('defocusDefault_psf_subj_%d', subIndex);
            else
                pdfFileName = sprintf('defocus%2.2fD_psf_subj_%d', imposedRefractionErrorDiopters, subIndex);
            end
            
            plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, pwd());
        end
    end
 end % plotPSFs
    