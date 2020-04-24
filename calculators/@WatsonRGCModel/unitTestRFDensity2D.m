function unitTestRFDensity2D()

    eccMinDegs = 1/60;
    eccMaxDegs = 2.5;
    eccSamplesNum = 32;
    eccDegs = logspace(log10(eccMinDegs), log10(eccMaxDegs), eccSamplesNum);
    
    obj = WatsonRGCModel();
    theView = 'right eye visual field';
    theView = 'left eye retina';
    
    [coneDensity2DMap, coneMeridianDensities, densitySupport, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = obj.compute2DConeRFDensity(eccDegs, theView);
    
    figNo = 1;
    doIt(figNo, obj, supportUnits, eccMinDegs, eccMaxDegs, densityUnits, ...
        densitySupport, coneDensity2DMap, coneMeridianDensities, horizontalMeridianLabel, verticalMeridianLabel, theView, densityLabel);
    
    
    [mRGCDensity2DMap, mRGCMeridianDensities, densitySupport, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = obj.compute2DmRGCRFDensity(eccDegs, theView, 'subtype', 'ON');
    
    densityLabel = strrep(densityLabel, 'mRGCRF', 'ON mRGC RF');
    figNo = 2;
    doIt(figNo, obj, supportUnits, eccMinDegs, eccMaxDegs, densityUnits, ...
        densitySupport, mRGCDensity2DMap, mRGCMeridianDensities, horizontalMeridianLabel, verticalMeridianLabel, theView, densityLabel);
    
    
    % cone to mRGC ratio map
    [conesToMRGCratio2DMap, spatialSupport, horizontalMeridianLabel, verticalMeridianLabel, ratioLabel, ...
            meridianConeToMRGratio, supportUnits] = obj.compute2DConeToMRGCRFRatio(eccDegs,  theView);
        
    figNo = 3;
    doIt(figNo, obj, supportUnits, eccMinDegs, eccMaxDegs, 'unitless', ...
        spatialSupport, conesToMRGCratio2DMap, meridianConeToMRGratio, horizontalMeridianLabel, verticalMeridianLabel, theView, ratioLabel);
    
end

function doIt(figNo, obj, supportUnits, eccMinDegs, eccMaxDegs, densityUnits, densitySupport, density2DMap, meridianDensities, horizontalMeridianLabel, verticalMeridianLabel, theView, densityLabel)
    if (strcmp(supportUnits, obj.visualDegsEccUnits))
        xLims = eccMaxDegs*[-1 1];
        xLimsPos = [eccMinDegs eccMaxDegs];
        xTicks = -100:2:100;
        xTickLabels = sprintf('%2.0f\n', xTicks);
        xTicksPos = [1/60 5/60 15/60 30/60 1 2 4 8 16];
        xTickPosLabels = {'1''',  '5''', '0.25', '0.5', '1', '2', '4', '8', '16'};
    else
        eccMinMM = obj.rhoDegsToMMs(eccMinDegs);
        eccMaxMM = obj.rhoDegsToMMs(eccMaxDegs);
        xLims = eccMaxMM*[-1 1];
        xLimsPos = [eccMinMM eccMaxMM];
        xTicks = -10:0.5:10;
        xTickLabels = {};
        for k = 1:numel(xTicks)
            if (mod(k,2) == 1)
                xTickLabels{k} = sprintf('%2.0f', xTicks(k));
            else
                xTickLabels{k} = '';
            end
        end
        xTicksPos = [0.01 0.03 0.1 0.3 1 3 10];
        for k = 1:numel(xTicksPos)
            if (xTicksPos(k) < 0.1)
                xTickPosLabels{k} = sprintf('%2.2f', xTicksPos(k));
            elseif (xTicksPos(k) < 1)
                xTickPosLabels{k} = sprintf('%2.1f', xTicksPos(k));
            else
                xTickPosLabels{k} = sprintf('%2.0f', xTicksPos(k));
            end
        end

    end
    

    if (strcmp(densityUnits, obj.visualDegsDensityUnits))
        densityLims = [50 20000];
        densityLevels = logspace(log10(densityLims(1)), log10(densityLims(2)), 15);
        scaling = 'log';
    elseif (strcmp(densityUnits, obj.retinalMMDensityUnits))
        densityLims = [50 20000];
        densityLims(1) = densityLims(1) / obj.alpha(0);
        densityLims(2) = densityLims(2) / obj.alpha(eccMaxDegs);
        densityLevels = logspace(log10(densityLims(1)), log10(densityLims(2)), 15);
        scaling = 'log';
    elseif (strcmp(densityUnits, 'unitless'))
        densityLims = [0 10];
        densityLevels = [1.05 1.25 1.5:0.5:10];
        scaling = 'linear';
    end
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1460 1600], 'Color', [1 1 1]);
    subplot(2,2,1);
    renderContourPlot(densitySupport, density2DMap, densityLevels, scaling, horizontalMeridianLabel, verticalMeridianLabel, xLims, xTicks, obj.figurePrefs.fontSize);
    title(sprintf('%s - %s', densityLabel, theView));
    
    subplot(2,2,3);
    midPoint = (size(density2DMap,1)-1)/2+1;
    renderProfilePlot(...
        squeeze(densitySupport(1,:)), squeeze(density2DMap(midPoint,:)), ...
        meridianDensities.ecc, meridianDensities.nasal, 'nasal', ...
        meridianDensities.ecc, meridianDensities.temporal, 'temporal', ...
        'log', horizontalMeridianLabel, densityLabel, xLims, xTicks, xTickLabels, densityLims, densityLevels, obj.figurePrefs.fontSize);
    
    subplot(2,2,2);
    midPoint = (size(density2DMap,1)-1)/2+1;
    renderProfilePlot(...
        squeeze(densitySupport(1,:)), squeeze(density2DMap(:,midPoint)), ....
        meridianDensities.ecc, meridianDensities.inferior, 'inferior', ...
        meridianDensities.ecc, meridianDensities.superior, 'superior', ...
        'log', verticalMeridianLabel, densityLabel, xLims, xTicks, xTickLabels, densityLims, densityLevels, obj.figurePrefs.fontSize);

    subplot(2,2,4);
    hold on
    theLegends = cell(1, numel(obj.enumeratedMeridianNames));
    for k = 1:numel(obj.enumeratedMeridianNames)
        meridianName = obj.enumeratedMeridianNames{k};
        theLegends{k} = sprintf('%s (%s)', meridianName, theView);
        switch (meridianName)
            case 'temporal meridian'
                plot(meridianDensities.ecc, log10(meridianDensities.temporal), '-', 'Color', obj.meridianColors(meridianName), 'LineWidth', 1.5);
            case 'nasal meridian'
                plot(meridianDensities.ecc, log10(meridianDensities.nasal), '-', 'Color', obj.meridianColors(meridianName), 'LineWidth', 1.5);
            case 'inferior meridian'
                plot(meridianDensities.ecc, log10(meridianDensities.inferior), '-', 'Color', obj.meridianColors(meridianName), 'LineWidth', 1.5);
            case 'superior meridian'
                plot(meridianDensities.ecc, log10(meridianDensities.superior), '-', 'Color', obj.meridianColors(meridianName), 'LineWidth', 1.5);
        end
    end
    legend(theLegends);
    axis 'square';
    xlabel(sprintf('eccentricity (%s)', supportUnits));
    ylabel(densityLabel);

    if (max(densityLevels)<20)
        yTickLabels = sprintf('%2.1f\n',densityLevels);
    else
        yTickLabels = sprintf('%2.0f\n',densityLevels);
    end
    
    set(gca, 'XLim', xLimsPos, 'YLim', log10(densityLims), ...
        'XTick', xTicksPos, ...
        'XTickLabel', xTickPosLabels, ...
        'YTick', log10(densityLevels), ...
        'YTickLabel', yTickLabels, ...
        'XScale', 'log', 'YScale', 'linear', ...
        'FontSize', obj.figurePrefs.fontSize);

    % Finalize figure
    grid(gca, 'on');
    
end

function renderProfilePlot(spatialSupport, profile, ...
    posSpatialSupport1, profile1, profile1Legend, ...
    posSpatialSupport2, profile2, profile2Legend, ...
    scaling, xLabelString, yLabelString, xLims, xTicks, xTickLabels, yLims, yTicks, fontSize)
    
    if (strcmp(scaling, 'log'))
        plot(spatialSupport, log10(profile), 'ko', 'lineWidth', 1.5, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 10); hold on;
        plot(posSpatialSupport1, log10(profile1), 'r-','lineWidth', 1.5);
        plot(posSpatialSupport2, log10(profile2), 'b-', 'lineWidth', 1.5);
        
        yLims = log10(yLims);
        if (max(yLims) < log10(20))
            yTickLabels = sprintf('%2.1f\n',yTicks);
        else
            yTickLabels = sprintf('%2.0f\n',yTicks);
        end
        yTicks = log10(yTicks);
    else
        plot(spatialSupport, profile, 'ko', 'lineWidth', 1.5, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 10);
        hold on;
        plot(posSpatialSupport1, profile1, 'r-', 'lineWidth', 1.5);
        plot(posSpatialSupport2, profile2, 'b-', 'lineWidth', 1.5);
        yTickLabels = sprintf('%2.0f\n',yTicks);
    end
    ylabel(yLabelString);
    
    legend({'full', profile1Legend, profile2Legend});
    axis 'square';
    
    xlabel(xLabelString);
    
    % Only positive
    set(gca, 'XLim', xLims, 'YLim', yLims, ...
        'XTick', xTicks, ...
        'XTickLabel', xTickLabels, ...
        'YTick', yTicks, ...
        'YTickLabel', yTickLabels, ...
        'FontSize', fontSize);

    % Finalize figure
    grid(gca, 'on');
end


function renderContourPlot(spatialSupport, Z, zLevels, scaling, xLabelString, yLabelString, xLims, xTicks, fontSize)
    X = squeeze(spatialSupport(1,:));
    Y = squeeze(spatialSupport(2,:));
    
    contourLabelSpacing = 4000;
    if (strcmp(scaling, 'log'))
        [cH, hH] = contourf(X,Y, log10(Z), log10(zLevels), 'LineColor', [0 0 0], ...
             'LineWidth', 2.0, 'ShowText', 'on', 'LabelSpacing', contourLabelSpacing);
    else
        [cH, hH] = contourf(X,Y, Z, zLevels, 'LineColor', [0 0 0], ...
             'LineWidth', 2.0, 'ShowText', 'on', 'LabelSpacing', contourLabelSpacing);
    end

    hold on;
    clabel(cH,hH,'FontWeight','normal', 'FontSize', 12, 'Color', [1 0 0], 'BackgroundColor', 'none')

    
    cMap = brewermap(1024, 'greys');
    axis 'square';
    colormap(cMap);
    
    xlabel(xLabelString); ylabel(yLabelString);
    
    
    set(gca, 'XLim', xLims, 'YLim', xLims,...
        'XTick', xTicks, ...
        'YTick', xTicks, ...
        'FontSize', fontSize);
    %colorbar('Ticks', log10(zLevels), 'TickLabels', sprintf('%2.0f\n', zLevels), 'Location', 'North');
    % Finalize figure
    grid(gca, 'on');
    axis(gca, 'square');
    drawnow;
    
    
    if (max(zLevels)>=1000) 
        % Label contours as %2.1k
        labels=hH.TextPrims;
        for idx = numel(labels):-1:1
            theZLevel=str2num(labels(idx).String);
            if (strcmp(scaling, 'log'))
                theZLevel = 10^theZLevel;
            end
            hH.TextPrims(idx).String=sprintf('%2.0fk', theZLevel/1000);
        end
    end
    
end
