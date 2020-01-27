function computeConnectivityMatrix

    % Instantiate a WatsonRGCModel object.
    WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);
    
    nSamples = 50; maxEcc = 30;
    xPos = logspace(log10(0.05), log10(maxEcc), nSamples);
    xPos = [-fliplr(xPos) xPos]; yPos = xPos;
    [X,Y] = meshgrid(xPos, yPos);
    eccXYposDegs = [X(:) Y(:)];
    
   

    % Points along the nasal meridian
    nasalIndices = find((eccXYposDegs(:,1) >= 0) & (eccXYposDegs(:,2) == min(abs(yPos))));

    % Points along the temporal meridian
    temporalIndices = find((eccXYposDegs(:,1) <= 0) & (eccXYposDegs(:,2) == min(abs(yPos))));
    
    % Points along the superior meridian
    superiorIndices = find((eccXYposDegs(:,2) >= 0) & (eccXYposDegs(:,1) == min(abs(xPos))));
    
    % Points along the inferior meridian
    inferiorIndices = find((eccXYposDegs(:,2) <= 0) & (eccXYposDegs(:,1) == min(abs(xPos))));
    
    
    % Meridian colors
    quadrantColors = containers.Map();
    quadrantColors('nasal')    = [0.0 0.7 0.0; 0.0 1.00 0.0];  % green
    quadrantColors('temporal') = [1.0 0.0 0.0; 1.0 0.50 0.50]; % red
    quadrantColors('superior') = [0.0 0.0 1.0; 0.5 0.5 1.0];   % blue
    quadrantColors('inferior') = [0.0 0.0 0.0; 0.5 0.5 0.5];   % gray
    
    % Meridian axes names
    xEccentricityAxisName = sprintf('space (deg)\n   <- nasal  |  temporal ->');
    yEccentricityAxisName = sprintf('space (deg)\n<- inferior  |  superior ->');
    
    
    conesPerRGCRFRatioLevels = [1:1:21]; 
    midgetRGCRFPerConeRatioLevels = [0 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
    
    logConeDensityLevels = [1:0.1:4];
    cMapConeDensity = brewermap(numel(logConeDensityLevels), '*greys');
    cMapConesPerRGCRatio = brewermap(numel(conesPerRGCRFRatioLevels), 'YlGnBu');
    cMapRGCRFPerConeRatio = brewermap(numel(conesPerRGCRFRatioLevels), '*YlGnBu');
    
    [midgetRGCRFToConeRatios, coneRFDensities, mRGCRFDensities] ...
        = WatsonRGCCalc.ratioOfMidgetRGCsToCones(eccXYposDegs, 'left');
    conesPerMidgetRGCRF = 1 ./ midgetRGCRFToConeRatios;
    nasalRatios = 1./midgetRGCRFToConeRatios(nasalIndices);
    temporalRatios = 1./midgetRGCRFToConeRatios(temporalIndices);
    superiorRatios = 1./midgetRGCRFToConeRatios(superiorIndices);
    inferiorRatios = 1./midgetRGCRFToConeRatios(inferiorIndices);
    
    
    
    figure(100);
    clf;
    
    subplot(3,4,1);
    makeContourPlot(log10(coneRFDensities), xPos, yPos, logConeDensityLevels, logConeDensityLevels(1:2:end), cMapConeDensity, ...
        'cone density (log(cones)/deg^2 left eye)', ...
        xEccentricityAxisName, yEccentricityAxisName);
    
    subplot(3,4,2);
    makeContourPlot(log10(mRGCRFDensities), xPos, yPos, logConeDensityLevels, logConeDensityLevels(1:2:end), cMapConeDensity, ...
        'mRGC RF density (log(mRGC)/deg^2 left eye)', ...
        xEccentricityAxisName, yEccentricityAxisName);
    
    
    subplot(3,4,3);
    makeContourPlot(midgetRGCRFToConeRatios, xPos, yPos, midgetRGCRFPerConeRatioLevels, midgetRGCRFPerConeRatioLevels(1:1:end), cMapRGCRFPerConeRatio, ...
        'mRGC RF / cones (left eye)', ...
        xEccentricityAxisName, yEccentricityAxisName);
    
    subplot(3,4,4);
    makeContourPlot(conesPerMidgetRGCRF, xPos, yPos, conesPerRGCRFRatioLevels, conesPerRGCRFRatioLevels(1:2:end), cMapConesPerRGCRatio, ...
        'cones / mRGC RF center (left eye)', ...
        xEccentricityAxisName, yEccentricityAxisName);
    
    subplot(3,4,8);
    makeLinePlot(eccXYposDegs(nasalIndices,1), nasalRatios, ...
                 eccXYposDegs(temporalIndices,1), temporalRatios, ...
                 quadrantColors, {'nasal' 'temporal'}, ...
                 maxEcc, -30:5:30, [1 20], 1:1:20, xEccentricityAxisName, 'cones / mRGC RF center');
             
    subplot(3,4,12);
    makeLinePlot(eccXYposDegs(superiorIndices,2), superiorRatios, ...
                 eccXYposDegs(inferiorIndices,2), inferiorRatios, ...
                 quadrantColors, {'superior' 'inferior'}, ...
                 maxEcc, -30:5:30, [1 20], 1:1:20, yEccentricityAxisName, 'cones / mRGC RF center');
             
             
    subplot(3,4,7);
    makeLinePlot(eccXYposDegs(nasalIndices,1), nasalRatios, ...
                 eccXYposDegs(temporalIndices,1), temporalRatios, ...
                 quadrantColors, {'nasal' 'temporal'}, ...
                 5, -5:1:5, [1 5.0], 1:0.5:5.0, xEccentricityAxisName, 'cones / mRGC RF center');
             
    subplot(3,4,11);
    makeLinePlot(eccXYposDegs(superiorIndices,2), superiorRatios, ...
                 eccXYposDegs(inferiorIndices,2), inferiorRatios, ...
                 quadrantColors, {'superior' 'inferior'}, ...
                 5, -5:1:5, [1 5.0], 1:0.5:5.0, yEccentricityAxisName, 'cones / mRGC RF center');
             
             
end


function makeContourPlot(ratioMap, xPos, yPos, zLevels, zLevelsOnCbar, cMap, titleName, xLabelName, yLabelName)
    rows = numel(yPos);
    cols = numel(xPos);
    [X,Y] = meshgrid(xPos, yPos);
    contourf(X,Y, reshape(ratioMap, [rows, cols]), zLevels);
    set(gca, 'XLim', [min(xPos) max(xPos)], 'YLim', [min(yPos) max(yPos)], 'FontSize', 14);
    axis 'image';
    colormap(gca,cMap)
    hL = colorbar;
    set(hL, 'Ticks', zLevelsOnCbar, 'Limits', [min(zLevels) max(zLevels)]);
    title(titleName);
    xlabel(xLabelName, 'FontName', 'Menlo', 'FontSize', 12, 'FontAngle', 'italic');
    ylabel(yLabelName, 'FontName', 'Menlo', 'FontSize', 12, 'FontAngle', 'italic');
end

function makeLinePlot(eccXYposDegs1, ratios1, eccXYposDegs2, ratios2, ...
                 quadrantColors, dataNames, ...
                 maxEcc, eccTicks, yRange, yTicks, xLabelName, yLabelName)
             
    color1 = quadrantColors(dataNames{1});
    color2 = quadrantColors(dataNames{2});
    
    plot(eccXYposDegs1, ratios1, ...
        'ko-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'Color', color1(1,:), ...
        'MarkerFaceColor', color1(2,:), ...
        'MarkerEdgeColor', color1(1,:));
    hold on;
    
    plot(eccXYposDegs2, ratios2, ...
        'ko-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'Color', color2(1,:), ...
        'MarkerFaceColor', color2(2,:), ...
        'MarkerEdgeColor', color2(1,:));
    
    xlabel(xLabelName, 'FontName', 'Menlo', 'FontAngle', 'italic');
    ylabel(yLabelName, 'FontName', 'Menlo', 'FontAngle', 'italic');
    hL = legend(dataNames, 'Location', 'North');
    axis 'square'
    set(gca, 'XLim', maxEcc*[-1 1], 'YLim', yRange, ...
         'XScale', 'linear', 'YScale', 'linear', ...
         'XTick', eccTicks, 'YTick', yTicks, ...
         'FontSize', 14);
    grid('on');
end



function oldConnectivityMatrix
    eccDegs = 0.0:0.05:30;
    meridians = {...
        'nasal meridian' ...
        'superior meridian' ...
        'temporal meridian' ...
        'inferior meridian' ...
    };
    
    for mIndex = 1:numel(meridians)
        s = coneDensityReadData(...
            'eccentricity',eccDegs, ...
            'eccentricityUnits', 'deg', ...
            'angle',(mIndex-1)*90 + eccDegs*0, ...
            'angleUnits', 'deg', ...
            'whichEye','left', ...
            'useParfor', true);
        coneDensityFromISETbio(mIndex,:) = s;
    end
    
    % Instantiate a WatsonRGCModel object.
    WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);
    
    for mIndex = 1:numel(meridians)
        meridian = meridians{mIndex};
        
        % The density of the ON or OFF midget RGC mosaic is half the density of
        % the total midgetRGC mosaic assuming equal numerosities of the ON and the OFF submosaics
        midgetRGCRFDensity(mIndex,:) = 0.5*WatsonRGCCalc.midgetRGCRFDensity(eccDegs, meridian, 'RFs per mm2');
        [coneSpacingDegs, coneDensity(mIndex,:)] = WatsonRGCCalc.coneRFSpacingAndDensity(eccDegs, meridian, 'Cones per mm2');

        midgetRGCRFSinglePolarityDensity = midgetRGCRFDensity;
        averageNumberOfConesPerMidgetRGCcenter(mIndex,:) = coneDensity(mIndex,:) ./ midgetRGCRFSinglePolarityDensity(mIndex,:);
        averageNumberOfISETBioConesPerMidgetRGCcenter(mIndex,:) = coneDensityFromISETbio(mIndex,:) ./ midgetRGCRFSinglePolarityDensity(mIndex,:);
    end
    
    figure(1); clf;
    mIndex = 1;
    subplot(2,3,1);
    plot(eccDegs, squeeze(midgetRGCRFDensity(mIndex,:)), 'r.-', 'LineWidth', 1.5);
    hold on;
    plot(eccDegs, squeeze(coneDensity(mIndex,:)), 'r.', 'LineWidth', 1.5);
    plot(eccDegs, squeeze(coneDensityFromISETbio(mIndex,:)), 'ro', 'LineWidth', 1.0);
    set(gca, 'XScale', 'log', 'YScale', 'log', 'YLim', [500 500*1000], 'YTick', [300 1000 3000 10000 30000 100000 300000]);
    grid on;
    ylabel('density (count/mm^2)');
    xlabel('eccentricity (degs)');
    legend({'RGC', 'cones', 'cones (isetbio)'});
    title(meridians{mIndex});
    
    subplot(2,3,2);
    mIndex = 2;
    plot(eccDegs, squeeze(midgetRGCRFDensity(mIndex,:)), 'b.-', 'LineWidth', 1.5);
    hold on;
    plot(eccDegs, squeeze(coneDensity(mIndex,:)), 'b.', 'LineWidth', 1.5);
    plot(eccDegs, squeeze(coneDensityFromISETbio(mIndex,:)), 'bo', 'LineWidth', 1.0);
    set(gca, 'XScale', 'log', 'YScale', 'log', 'YLim', [500 500*1000], 'YTick', [300 1000 3000 10000 30000 100000 300000]);
    grid on;
    ylabel('density (count/mm^2)');
    xlabel('eccentricity (degs)');
    legend({'RGC', 'cones', 'cones (isetbio)'});
    title(meridians{mIndex});
    
    subplot(2,3,4);
    mIndex = 3;
    plot(eccDegs, squeeze(midgetRGCRFDensity(mIndex,:)), 'm.-', 'LineWidth', 1.5);
    hold on;
    plot(eccDegs, squeeze(coneDensity(mIndex,:)), 'm.', 'LineWidth', 1.5);
    plot(eccDegs, squeeze(coneDensityFromISETbio(mIndex,:)), 'mo', 'LineWidth', 1.0);
    set(gca, 'XScale', 'log', 'YScale', 'log', 'YLim', [500 500*1000], 'YTick', [300 1000 3000 10000 30000 100000 300000]);
    grid on;
    ylabel('density (count/mm^2)');
    xlabel('eccentricity (degs)');
    legend({'RGC', 'cones', 'cones (isetbio)'});
    title(meridians{mIndex});
    
    subplot(2,3,5);
    mIndex = 4;
    plot(eccDegs, squeeze(midgetRGCRFDensity(mIndex,:)), 'k.-', 'LineWidth', 1.5);
    hold on;
    plot(eccDegs, squeeze(coneDensity(mIndex,:)), 'k.', 'LineWidth', 1.5);
    plot(eccDegs, squeeze(coneDensityFromISETbio(mIndex,:)), 'ko', 'LineWidth', 1.0);
    set(gca, 'XScale', 'log', 'YScale', 'log', 'YLim', [500 500*1000], 'YTick', [300 1000 3000 10000 30000 100000 300000]);
    grid on;
    ylabel('density (count/mm^2)');
    xlabel('eccentricity (degs)');
    legend({'RGC', 'cones', 'cones (isetbio)'});
    title(meridians{mIndex});
    
    
    subplot(2,3,3);
    plot(eccDegs, squeeze(averageNumberOfConesPerMidgetRGCcenter(1,:)), 'ro-', 'LineWidth', 1.5);
    hold on;
    plot(eccDegs, squeeze(averageNumberOfConesPerMidgetRGCcenter(2,:)), 'bo-', 'LineWidth', 1.5);
    plot(eccDegs, squeeze(averageNumberOfConesPerMidgetRGCcenter(3,:)), 'mo-', 'LineWidth', 1.5);
    plot(eccDegs, squeeze(averageNumberOfConesPerMidgetRGCcenter(4,:)), 'ko-', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log', 'YScale', 'log', 'YLim', [0.9 30], 'YTick', [1 1.5 2 3 5 7 10 15 20 30], 'XTick', [0.1 0.3 1 3 10 30]);
    grid on;
    ylabel('cones per midget RGC');
    xlabel('eccentricity (degs)');
    legend(meridians);
    
    subplot(2,3,6);
    plot(eccDegs, squeeze(averageNumberOfISETBioConesPerMidgetRGCcenter(1,:)), 'ro-', 'LineWidth', 1.5);
    hold on;
    plot(eccDegs, squeeze(averageNumberOfISETBioConesPerMidgetRGCcenter(2,:)), 'bo-', 'LineWidth', 1.5);
    plot(eccDegs, squeeze(averageNumberOfISETBioConesPerMidgetRGCcenter(3,:)), 'mo-', 'LineWidth', 1.5);
    plot(eccDegs, squeeze(averageNumberOfISETBioConesPerMidgetRGCcenter(4,:)), 'ko-', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log', 'YScale', 'log', 'YLim', [0.9 30], 'YTick', [1 1.5 2 3 5 7 10 15 20 30], 'XTick', [0.1 0.3 1 3 10 30]);
    grid on;
    ylabel('ISETBio cones per midget RGC');
    xlabel('eccentricity (degs)');
    legend(meridians);
    
    
end

