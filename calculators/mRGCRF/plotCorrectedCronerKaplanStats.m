function plotCorrectedCronerKaplanStats

    ck = CronerKaplanRGCModel();

    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', [1 0 0; 0 0 1], ...
            'axesBox', 'off', ...
            'axesTickDir', 'in', ...
            'renderer', 'painters', ...
            'axesTickLength', [0.01 0.01], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 14, ...
            'figureHeightInches', 14);

        
    eccDegs = logspace(log10(0.1), log10(25), 16);
    centerVisualRadiusDegs = ck.centerRadiusFunction(ck.centerRadiusParams, eccDegs);
    minRFradius = 0.013 % 0.6*min(ck.centerData('size').radiusDegs)
    centerVisualRadiusDegs(centerVisualRadiusDegs<minRFradius) = minRFradius;
    surroundVisualRadiusDegs = ck.surroundRadiusFunction(ck.surroundRadiusParams, eccDegs);
    centerVisualPeakSensitivity = ck.centerPeakSensitivityFunction(ck.centerPeakSensitivityParams, centerVisualRadiusDegs);
    surroundVisualPeakSensitivity = ck.surroundPeakSensitivityFunction(ck.surroundPeakSensitivityParams, surroundVisualRadiusDegs);
    w = WatsonRGCModel();
    coneRFSpacingDeg  = w.coneRFSpacingAndDensityAlongMeridian(eccDegs, ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio',false);
        
    
    defocusMode = 'subjectDefault';
    if (strcmp(defocusMode, 'subjectDefault'))
        matFileName = 'defocusDefaultFittedModel.mat';
    else
        matFileName = 'defocusNoneFittedModel.mat';
    end

    load(matFileName, 'retinalPoolingRadii', 'eccTested', 'fittedParamsRadius', 'fittedParamsGain', 'modelFunctionRadius', 'modelFunctionGain');

    
        
    hFig = figure(2); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', 2, ...
            'colsNum', 2, ...
            'leftMargin', 0.08, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.1, ...
            'heightMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'topMargin', 0.04);
        
    plotDataForSubregion(theAxesGrid{1,1}, theAxesGrid{2,1}, eccDegs, centerVisualRadiusDegs, centerVisualPeakSensitivity, ...
        eccTested, fittedParamsRadius, fittedParamsGain, modelFunctionRadius, modelFunctionGain, coneRFSpacingDeg,  [1 0 0]);


    plotDataForSubregion(theAxesGrid{1,2}, theAxesGrid{2,2}, eccDegs, surroundVisualRadiusDegs, surroundVisualPeakSensitivity, ...
        eccTested, fittedParamsRadius, fittedParamsGain, modelFunctionRadius, modelFunctionGain, coneRFSpacingDeg, [0 0 1]);
    
    plotlabOBJ.exportFig(hFig, 'pdf', 'ModelCorrections', pwd());
    
    pause
    

    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'lineColor', [0.5 0.5 0.5], ...
            'scatterMarkerEdgeColor', [0.3 0.3 0.9], ...
            'lineMarkerSize', 11, ...
            'axesBox', 'off', ...
            'axesTickDir', 'in', ...
            'renderer', 'painters', ...
            'axesTickLength', [0.02 0.02], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 28, ...
            'figureHeightInches', 16.5);
     
    
    eccDegs = [0.25:0.05:1 1.1:0.1:2 3 4 5 6 7 8 9 10];
    centerVisualRadiusDegs = ck.centerRadiusFunction(ck.centerRadiusParams, eccDegs);
    minRFradius = 0.013 % 0.6*min(ck.centerData('size').radiusDegs)
    centerVisualRadiusDegs(centerVisualRadiusDegs<minRFradius) = minRFradius;
    surroundVisualRadiusDegs = ck.surroundRadiusFunction(ck.surroundRadiusParams, eccDegs);
    centerVisualPeakSensitivity = ck.centerPeakSensitivityFunction(ck.centerPeakSensitivityParams, centerVisualRadiusDegs);
    surroundVisualPeakSensitivity = ck.surroundPeakSensitivityFunction(ck.surroundPeakSensitivityParams, surroundVisualRadiusDegs);
    w = WatsonRGCModel();
    coneRFSpacingDeg  = w.coneRFSpacingAndDensityAlongMeridian(eccDegs, ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio',false);
        
    plotMeanRFs(1, eccDegs, centerVisualRadiusDegs, centerVisualPeakSensitivity, ...
        surroundVisualRadiusDegs, surroundVisualPeakSensitivity, ...
        eccTested, fittedParamsRadius, fittedParamsGain, modelFunctionRadius, modelFunctionGain, ...
        coneRFSpacingDeg, 'normalized', plotlabOBJ );
%     pause
    
%     plotMeanRFs(2, eccDegs, centerVisualRadiusDegs, centerVisualPeakSensitivity, ...
%         surroundVisualRadiusDegs, surroundVisualPeakSensitivity, ...
%         eccTested, fittedParamsRadius, fittedParamsGain, modelFunctionRadius, modelFunctionGain, ...
%         coneRFSpacingDeg, 'actual', plotlabOBJ );
    
    plotMean2DRFs(3, eccDegs, centerVisualRadiusDegs, centerVisualPeakSensitivity, ...
        surroundVisualRadiusDegs, surroundVisualPeakSensitivity, ...
        eccTested, fittedParamsRadius, fittedParamsGain, modelFunctionRadius, modelFunctionGain, ...
        coneRFSpacingDeg, 'actual', plotlabOBJ );
end

function plotMean2DRFs(figNo, eccDegs, centerVisualRadiusDegs, centerVisualPeakSensitivity, ...
        surroundVisualRadiusDegs, surroundVisualPeakSensitivity, ...
        fittedEcc, fittedParamsRadius, fittedParamsGain, modelFunctionRadius, modelFunctionGain, ...
        coneRFSpacingDeg, scalingMode,plotlabOBJ  )
    
    videoOBJ = VideoWriter('rfs2DVideo', 'MPEG-4'); % H264 format
        videoOBJ.FrameRate = 5;
        videoOBJ.Quality = 100;
        videoOBJ.open();
        
    for k = 1:numel(eccDegs)
        
        hFig = figure(figNo*100+k); clf;
        [~,eccIndex] = min(abs(eccDegs(k)-abs(fittedEcc)));
        
        centerRetinalRadiusDegs(k) = modelFunctionRadius(fittedParamsRadius(eccIndex,:), centerVisualRadiusDegs(k));
        
        surroundRetinalRadiusDegs(k) = modelFunctionRadius(fittedParamsRadius(eccIndex,:), surroundVisualRadiusDegs(k)); 
        centerRetinalPeakSensitivity(k) = centerVisualPeakSensitivity(k) / modelFunctionGain(fittedParamsGain(eccIndex,:),  centerRetinalRadiusDegs(k));
        surroundRetinalPeakSensitivity(k) = surroundVisualPeakSensitivity(k) / modelFunctionGain(fittedParamsGain(eccIndex,:),  surroundRetinalRadiusDegs(k));
        ax = subplot('Position', [0.01 0.01 0.98 0.98]);
        
        render2DDogRF(ax,eccDegs(k), coneRFSpacingDeg(k) , centerVisualRadiusDegs(k), centerVisualPeakSensitivity(k), surroundVisualRadiusDegs(k), surroundVisualPeakSensitivity(k), ...
                       centerRetinalRadiusDegs(k), centerRetinalPeakSensitivity(k), surroundRetinalRadiusDegs(k), surroundRetinalPeakSensitivity(k), scalingMode);
    
        plotlabOBJ.exportFig(hFig, 'pdf', sprintf('rf2D_%2.1fdegs', eccDegs(k)), pwd());
        videoOBJ.writeVideo(getframe(hFig));
    end
    videoOBJ.close();
end


function plotMeanRFs(figNo, eccDegs, centerVisualRadiusDegs, centerVisualPeakSensitivity, ...
        surroundVisualRadiusDegs, surroundVisualPeakSensitivity, ...
        fittedEcc, fittedParamsRadius, fittedParamsGain, modelFunctionRadius, modelFunctionGain, ...
        coneRFSpacingDeg, scalingMode,plotlabOBJ  )
    
    videoOBJ = VideoWriter('rfsVideo', 'MPEG-4'); % H264 format
        videoOBJ.FrameRate = 5;
        videoOBJ.Quality = 100;
        videoOBJ.open();
        
    for k = 1:numel(eccDegs)
        
        hFig = figure(figNo*100+k); clf;
        [~,eccIndex] = min(abs(eccDegs(k)-abs(fittedEcc)));
        
        centerRetinalRadiusDegs(k) = modelFunctionRadius(fittedParamsRadius(eccIndex,:), centerVisualRadiusDegs(k));
        
        surroundRetinalRadiusDegs(k) = modelFunctionRadius(fittedParamsRadius(eccIndex,:), surroundVisualRadiusDegs(k)); 
        centerRetinalPeakSensitivity(k) = centerVisualPeakSensitivity(k) / modelFunctionGain(fittedParamsGain(eccIndex,:),  centerRetinalRadiusDegs(k));
        surroundRetinalPeakSensitivity(k) = surroundVisualPeakSensitivity(k) / modelFunctionGain(fittedParamsGain(eccIndex,:),  surroundRetinalRadiusDegs(k));
        ax = subplot('Position', [0.01 0.01 0.98 0.98]);
        
        renderDogRF(ax,eccDegs(k), coneRFSpacingDeg(k) , centerVisualRadiusDegs(k), centerVisualPeakSensitivity(k), surroundVisualRadiusDegs(k), surroundVisualPeakSensitivity(k), ...
                       centerRetinalRadiusDegs(k), centerRetinalPeakSensitivity(k), surroundRetinalRadiusDegs(k), surroundRetinalPeakSensitivity(k), scalingMode);
    
        plotlabOBJ.exportFig(hFig, 'pdf', sprintf('rf_%2.1fdegs_%s', eccDegs(k), scalingMode), pwd());
        videoOBJ.writeVideo(getframe(hFig));
    end
    videoOBJ.close();
end

    
function renderDogRF(ax, eccDegs, coneRFSpacingDeg , centerVisualRadiusDegs, centerVisualPeakSensitivity, surroundVisualRadiusDegs, surroundVisualPeakSensitivity, ...
                         centerRetinalRadiusDegs, centerRetinalPeakSensitivity, surroundRetinalRadiusDegs, surroundRetinalPeakSensitivity, scalingMode)
    x = linspace(0, centerVisualRadiusDegs*3, 100);
    xC = [fliplr(-x(2:end)) x];
    visualRFcenter = centerVisualPeakSensitivity * exp(-(xC/centerVisualRadiusDegs).^2);
    retinalRFcenter = centerRetinalPeakSensitivity * exp(-(xC/centerRetinalRadiusDegs).^2);
    x = linspace(0, surroundVisualRadiusDegs*2, 100);
    xS = [fliplr(-x(2:end)) x];
    visualRFsurround = -surroundVisualPeakSensitivity * exp(-(xS/surroundVisualRadiusDegs).^2);
    retinalRFsurround = -surroundRetinalPeakSensitivity * exp(-(xS/surroundRetinalRadiusDegs).^2);
    
    coneRect = [-1 -1 1 1 -1; 0 1 1 0 0];
    coneRect(1,:) = coneRect(1,:) * coneRFSpacingDeg/2;
    xLims = [-0.8 0.8];
    hold(ax, 'on');
    cs1 = (centerVisualRadiusDegs/surroundVisualRadiusDegs)^2 * centerVisualPeakSensitivity/surroundVisualPeakSensitivity;
    cs2 = (centerRetinalRadiusDegs/surroundRetinalRadiusDegs)^2 * centerRetinalPeakSensitivity/surroundRetinalPeakSensitivity;
    
    if (strcmp(scalingMode, 'normalized'))
        yLim = [-70 1300];
        yTicks = -0.02:0.02:1;
        text(-0.2, 40, sprintf('%2.2f^o', eccDegs), 'FontSize', 20);
        text(0.2, 40, sprintf('S/C int. sens. ratio: %2.2f (visual)/ %2.2f (retinal)', cs1, cs2), 'FontSize', 20);
        coneRect(2,:) = coneRect(2,:) * max(retinalRFcenter)*exp(-1);
    else
        yLim = [-70 13000];
        yTicks = -100:50:5000;
        coneRect(2,:) = coneRect(2,:) * max(retinalRFcenter)*exp(-1);
        text(-0.2, 400, sprintf('%2.2f^o', eccDegs), 'FontSize', 20);
        text(0.2, 400, sprintf('S/C int. sens. ratio: %2.2f (visual) / %2.2f (retinal)', 1/cs1, 1/cs2), 'FontSize', 20);
    end
    
    line(ax, xC, visualRFcenter, 'Color', [0 0 0], 'LineStyle', '-');
    area(ax, xC, retinalRFcenter, 'FaceColor', 'r');
    line(ax, xS, visualRFsurround, 'Color', [0 0 0], 'LineStyle', '-');
    area(ax, xS, retinalRFsurround, 'FaceColor', 'b');
    patch(ax, coneRect(1,:), coneRect(2,:), [0 0.5 0.7], 'LineWidth', 1, 'FaceAlpha', 0.2);
    grid(ax, 'off');
    xTicksArcMin = -60:0.5:60;
   
    set(ax, 'XLim', xLims, 'XTick', xTicksArcMin/60,  'YLim', yLim, 'YTick', yTicks,...
        'XColor', 'none', 'YColor', 'none');
end

function render2DDogRF(ax, eccDegs, coneRFSpacingDeg , centerVisualRadiusDegs, centerVisualPeakSensitivity, surroundVisualRadiusDegs, surroundVisualPeakSensitivity, ...
                         centerRetinalRadiusDegs, centerRetinalPeakSensitivity, surroundRetinalRadiusDegs, surroundRetinalPeakSensitivity, scalingMode)
    
                     
    xLims = [-0.5 0.5];
    hold(ax, 'on');
    xSv = surroundVisualRadiusDegs * cosd(0:360);
    ySv = surroundVisualRadiusDegs * sind(0:360);
    xS = surroundRetinalRadiusDegs * cosd(0:360);
    yS = surroundRetinalRadiusDegs * sind(0:360);
    patch(ax,xS, yS, [0.5 0.5 1], 'LineWidth', 1, 'FaceAlpha', 0.2);
    
    xCv = centerVisualRadiusDegs * cosd(0:360);
    yCv = centerVisualRadiusDegs * sind(0:360);
    xC = centerRetinalRadiusDegs * cosd(0:360);
    yC = centerRetinalRadiusDegs * sind(0:360);
    
    patch(ax,xC, yC, [1 0.5 0.5], 'LineWidth', 1, 'FaceAlpha', 0.2);
    
    line(ax, xCv, yCv, 'LineStyle', '--', 'LineColor', 'k');
    line(ax, xSv, ySv, 'LineStyle', '--', 'LineColor', 'k');
    grid(ax, 'off');
    axis(ax, 'square');
    box(ax, 'on');
    set(ax, 'XLim', xLims, 'YLim', xLims, 'XTick', [], 'YTick', []);
end


function plotDataForSubregion(theAxes1, theAxes2, eccDegs, CronnerKaplanModelForSubregionVisualRadiusDegs, CronnerKaplanModelForSubregionVisualPeakSensitivity, ...
    eccTested, fittedParamsRadius, fittedParamsGain, modelFunctionRadius, modelFunctionGain, coneRFSpacingDeg, color)
    
    % Compute correction for subregion radius
    for iecc = 1:numel(eccDegs)
        % find nearest ecc in eccTested
        [~,eccIndex] = min(abs(eccDegs(iecc)-abs(eccTested)));
        retinalRadiusDegs(iecc) = modelFunctionRadius(fittedParamsRadius(eccIndex,:), CronnerKaplanModelForSubregionVisualRadiusDegs(iecc));
    end
    radiusCorrection = retinalRadiusDegs ./ CronnerKaplanModelForSubregionVisualRadiusDegs;
    
    
    
        
    theAxes = theAxes1;
    line(theAxes, eccDegs, coneRFSpacingDeg*0.5, 'LineStyle', '--', 'Color', 'k');
    hold(theAxes, 'on');
    scatter(theAxes, eccDegs, CronnerKaplanModelForSubregionVisualRadiusDegs, 'o', 'MarkerFaceColor',  color, 'MarkerEdgeColor', color);
    scatter(theAxes, eccDegs,  CronnerKaplanModelForSubregionVisualRadiusDegs.*radiusCorrection, 'd', 'MarkerFaceColor', color/2, 'MarkerEdgeColor', color);
    
    xlabel(theAxes, 'eccentricity (degs)');
    ylabel(theAxes, 'radius (degs)')
    set(theAxes, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.1 100], 'YLim', [0.003 1]);
    set(theAxes, 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], 'YTick', [0.0001 0.003 0.01 0.03 0.1 0.3 1 3 10]);

   
    
    % Now the visual gain attenuation
    for iecc = 1:numel(eccDegs)
        [~,eccIndex] = min(abs(eccDegs(iecc)-abs(eccTested)));
        gainAttenuation(iecc) = modelFunctionGain(fittedParamsGain(eccIndex,:), retinalRadiusDegs(iecc));
    end
    gainCorrection = 1./gainAttenuation;
    
    theAxes = theAxes2;
    hold(theAxes, 'on');
    scatter(theAxes, eccDegs, CronnerKaplanModelForSubregionVisualPeakSensitivity, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color);
    scatter(theAxes, eccDegs, CronnerKaplanModelForSubregionVisualPeakSensitivity.*gainCorrection, 'd', 'MarkerFaceColor', color/2, 'MarkerEdgeColor', color);
   
    xlabel(theAxes, 'ecc (degs)');
    ylabel(theAxes, 'peak sensitivity');
    set(theAxes, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.1 100], 'YLim', [0.3 30000]);
    set(theAxes, 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], 'YTick',  [0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000 10000 30000 100000]);
    

end

