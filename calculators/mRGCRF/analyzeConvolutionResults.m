function analyzeConvolutionResults
    
    defocusMode = 'subjectDefault';
    if (strcmp(defocusMode, 'subjectDefault'))
        imposedRefractionErrorDiopters = 0; 
    else
        imposedRefractionErrorDiopters =  0.01; 
    end
    
    eccTested = -[0 1 2.5 5 7.5 10 15 20 25];
    plotRows = 2;
    plotCols = 5;
    
    plotlabOBJ = setupPlotLab([28 14], 'in');
    cMap = brewermap(512, 'greys');
    
    % Display the employed PSFs
    goodSubjects = [5 6 8 10]; %[5 6 7 8 9 10]; - 7 -9
    goodSubjectIndices = [goodSubjects goodSubjects+10 goodSubjects+20];
     
    plotPSFs = false;
    if (plotPSFs)
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
            
        groupIndices = [3 1 2];
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

    plotlabOBJ = setupPlotLab([20 10], 'both');
    
    % Display the visual radius data
    hFig = figure(100); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', plotRows, ...
            'colsNum', plotCols, ...
            'leftMargin', 0.08, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.03, ...
            'heightMargin', 0.08, ...
            'bottomMargin', 0.06, ...
            'topMargin', 0.06);
        
    for eccIndex = 1:numel(eccTested)
        eccDegs = eccTested(eccIndex);
        dataFileName = sprintf('cellEcc_%2.1f_cellRefractionError_%2.2fD_VisualGain.mat', eccDegs, imposedRefractionErrorDiopters);
        load(dataFileName, 'retinalPoolingRadii', 'visualRadius');
      
        visualRadius = visualRadius(:,goodSubjectIndices);
        row = floor((eccIndex-1)/plotCols)+1;
        col = mod(eccIndex-1, plotCols)+1;

    
        medianVisualRadius = (median(visualRadius,2, 'omitnan'))';

        % Fit the visualRadius with a saturating function
         modelFunction = @(p,x)(p(1) - (p(2)-p(1))*exp(-p(3)*x));
        initialParams = [5 10 0.2];
        
        [fittedParams, fittedParamsSE] = nonLinearFitData([medianVisualRadius 1], [retinalPoolingRadii 1], modelFunction, initialParams);
        fittedParamsRadius = fittedParams
        medianVisualRadiusFit = logspace(log10(0.01), log10(1), 100);
        retinalPoolingRadiiFit = modelFunction(fittedParams, medianVisualRadiusFit);
      
        
        yLims = [0.01 0.6];
        showIdentityLine = true;
        plotData(theAxesGrid{row,col}, eccDegs, retinalPoolingRadii, visualRadius,  medianVisualRadius, ...
            retinalPoolingRadiiFit, medianVisualRadiusFit, ...
            row==plotRows, col==1, imposedRefractionErrorDiopters, yLims, showIdentityLine, 'visual pooling radius');
           
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
            'leftMargin', 0.08, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.03, ...
            'heightMargin', 0.08, ...
            'bottomMargin', 0.06, ...
            'topMargin', 0.06);
        
    for eccIndex = 1:numel(eccTested)
        eccDegs = eccTested(eccIndex);
        dataFileName = sprintf('cellEcc_%2.1f_cellRefractionError_%2.2fD_VisualGain.mat', eccDegs, imposedRefractionErrorDiopters);
        load(dataFileName, 'retinalPoolingRadii',  'visualGain');
        visualGain = visualGain(:,goodSubjectIndices);
        
        row = floor((eccIndex-1)/plotCols)+1;
        col = mod(eccIndex-1, plotCols)+1;
        
        medianVisualGain = (median(visualGain,2, 'omitnan'))';
        
        % Fit the visualGain with a saturating function
         modelFunction = @(p,x)(p(1) - (p(2)-p(1))*exp(-p(3)*x));
        initialParams = [0.9 2 26];
        [fittedParams, fittedParamsSE] = nonLinearFitData([retinalPoolingRadii 1], [medianVisualGain 1], modelFunction, initialParams);
        fittedParamsGain = fittedParams
        retinalPoolingRadiiFit = logspace(log10(0.01), log10(1), 100);
        medianVisualGainFit = modelFunction(fittedParams, retinalPoolingRadiiFit);
        
        yLims = [0 1];
        showIdentityLine = false;
        plotData(theAxesGrid{row,col}, eccDegs, retinalPoolingRadii, visualGain,  medianVisualGain, ...
            retinalPoolingRadiiFit, medianVisualGainFit, ...
            row==plotRows, col==1, imposedRefractionErrorDiopters, yLims, showIdentityLine, 'peak sensitivity attenuation');
    end
    
    if (strcmp(defocusMode, 'subjectDefault'))
        pdfFileName = sprintf('defocusDefault_gain');
    else
        pdfFileName = sprintf('defocus%2.2fD_gain', imposedRefractionErrorDiopters);
    end
            
    plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, pwd());
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
    showXLabel, showYLabel, imposedRefractionErrorDiopters, yLims, showIdentityLine, theXLabel)
    
    hold(theAxes, 'on');
    if (showIdentityLine)
        line(theAxes, [0.01 1], [0.01 1], 'Color', [0 0 0], 'LineStyle', ':');
    end

    usedSubjectsNum = size(visualQuantity,2)/3;
    groupColors = [0.4 0.4 0.4; 1 0.4 0.4; 0.4 0.4 1];
    for group = 3:-1:1
        color = groupColors(group,:);
        for subIndex = 1:usedSubjectsNum
            line(theAxes, xQuantity, visualQuantity(:,(group-1)*usedSubjectsNum+subIndex), 'Color', color); 
        end
    end
    plot(theAxes, xQuantity, visualQuantityMedian, 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgecolor', [0 1 1]);
    plot(theAxes, fittedXQuantity, fittedVisualQuantityMedian, 'b-', 'LineWidth',4.0);
    plot(theAxes, fittedXQuantity, fittedVisualQuantityMedian, 'c-', 'LineWidth', 2.0);
    
    axis(theAxes, 'square');
    set(theAxes, 'XScale', 'log', 'YScale', 'log');
    set(theAxes, 'XTick', [ 0.01 0.03 0.1 0.3 0.6], 'YTick', [0.01 0.03 0.1 0.3 0.6 1]);
    set(theAxes, 'XLim', [0.01 0.6], 'YLim', yLims);
    
    if (imposedRefractionErrorDiopters == 0)
        title(theAxes, sprintf('cell eccentricity: %2.1f degs\n(subject defocus)', ecc));
    else
        title(theAxes, sprintf('cell eccentricity: %2.1f degs\n(defocus error: %2.2fD)', ecc, imposedRefractionErrorDiopters));
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
