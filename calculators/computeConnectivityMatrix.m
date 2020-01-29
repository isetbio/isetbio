function computeConnectivityMatrix

    % Instantiate a WatsonRGCModel object.
    WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);
    
    % Sample the retinal space
    samplesPerMeridian = 50; maxEcc = 30; spacing = 'log';
    [eccXYposDegs, quadrants] = sampleRetinalSpace(WatsonRGCCalc,samplesPerMeridian, maxEcc, spacing);
    
    % Compute midgetRGCRF to cone ratios at the above eccentricities
    whichEye = 'left';
    [midgetRGCRFToConeRatios, coneRFDensities, mRGCRFDensities] ...
        = WatsonRGCCalc.ratioOfMidgetRGCsToCones(eccXYposDegs, whichEye);
    
    % Compute cones per midget RGC
    conesPerMidgetRGCRF = 1 ./ midgetRGCRFToConeRatios;
    
    
    % Compute midget spacing using equation (9) in the Watson (2014) paper.
    singlePolarityMidgetRGCRFDensities = 0.5*mRGCRFDensities; % only ON or OFF
    singlePolarityMidgetRGCRFSpacings = sqrt(2.0./(sqrt(3.0)*singlePolarityMidgetRGCRFDensities));
     
    % Compute cone spacing from density
    coneSpacings = sqrt(2.0./(sqrt(3.0)*coneRFDensities));
    
    % Plot cones per mRGC along meridians
    figureNo = 1;
    plotMeridianStats(figureNo, conesPerMidgetRGCRF,eccXYposDegs, quadrants, maxEcc, [1 20], 1:1:20, 'cones / mRGC RF center');
    
    % Plot mRGC RF spacing along meridians
    figureNo = 2;
    plotMeridianStats(figureNo, singlePolarityMidgetRGCRFSpacings,eccXYposDegs, quadrants, 16, [0 0.16], 0:0.04:0.16, 'mRGC RF spacing');
    
    % Plot cone spacing along meridians
    figureNo = 3;
    plotMeridianStats(figureNo, coneSpacings,eccXYposDegs, quadrants, 16, [0 0.16], 0:0.04:0.16, 'cone spacing');
    
end

function plotMeridianStats(figureNo, stats,eccXYposDegs, quadrants, maxEcc, statsRange, statsTicks, statsName)
    nasalStats = stats(quadrants('nasal').meridianIndices);
    nasalEccDegs = eccXYposDegs(quadrants('nasal').meridianIndices,1);
    
    temporalStats = stats(quadrants('temporal').meridianIndices);
    temporalEccDegs = eccXYposDegs(quadrants('temporal').meridianIndices,1);
    
    superiorStats = stats(quadrants('superior').meridianIndices);
    superiorEccDegs = eccXYposDegs(quadrants('superior').meridianIndices,2);
    
    inferiorStats = stats(quadrants('inferior').meridianIndices);
    inferiorEccDegs = eccXYposDegs(quadrants('inferior').meridianIndices,2);
    
    % Meridian axes names
    xEccentricityAxisName = sprintf('space (deg)\n   <- nasal  |  temporal ->');
    yEccentricityAxisName = sprintf('space (deg)\n<- inferior  |  superior ->');
    
    figure(figureNo); clf;
    subplot(1,2,1);
    makeLinePlot(nasalEccDegs, nasalStats, ...
                 temporalEccDegs, temporalStats, ...
                 quadrants('nasal').colors, quadrants('temporal').colors, ...
                 {'nasal' 'temporal'}, ...
                 maxEcc, -30:5:30, statsRange, statsTicks, xEccentricityAxisName, statsName);
             
    subplot(1,2,2);
    makeLinePlot(superiorEccDegs, superiorStats, ...
                 inferiorEccDegs, inferiorStats, ...
                 quadrants('superior').colors, quadrants('inferior').colors, ...
                 {'superior' 'inferior'}, ...
                 maxEcc, -30:5:30, statsRange, statsTicks, yEccentricityAxisName, statsName);
             
             
end


    
%     
%     
%     conesPerRGCRFRatioLevels = [1:1:21]; 
%     midgetRGCRFPerConeRatioLevels = [0 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
%     
%     logConeDensityLevels = [1:0.1:4];
%     cMapConeDensity = brewermap(numel(logConeDensityLevels), '*greys');
%     cMapConesPerRGCRatio = brewermap(numel(conesPerRGCRFRatioLevels), 'YlGnBu');
%     cMapRGCRFPerConeRatio = brewermap(numel(conesPerRGCRFRatioLevels), '*YlGnBu');
%     
%     
%     
%     
%     
%     figure(100);
%     clf;
%     
%     subplot(3,4,1);
%     makeContourPlot(log10(coneRFDensities), xPos, yPos, logConeDensityLevels, logConeDensityLevels(1:2:end), cMapConeDensity, ...
%         'cone density (log(cones)/deg^2 left eye)', ...
%         xEccentricityAxisName, yEccentricityAxisName);
%     
%     subplot(3,4,2);
%     makeContourPlot(log10(mRGCRFDensities), xPos, yPos, logConeDensityLevels, logConeDensityLevels(1:2:end), cMapConeDensity, ...
%         'mRGC RF density (log(mRGC)/deg^2 left eye)', ...
%         xEccentricityAxisName, yEccentricityAxisName);
%     
%     
%     subplot(3,4,3);
%     makeContourPlot(midgetRGCRFToConeRatios, xPos, yPos, midgetRGCRFPerConeRatioLevels, midgetRGCRFPerConeRatioLevels(1:1:end), cMapRGCRFPerConeRatio, ...
%         'mRGC RF / cones (left eye)', ...
%         xEccentricityAxisName, yEccentricityAxisName);
%     
%     subplot(3,4,4);
%     makeContourPlot(conesPerMidgetRGCRF, xPos, yPos, conesPerRGCRFRatioLevels, conesPerRGCRFRatioLevels(1:2:end), cMapConesPerRGCRatio, ...
%         'cones / mRGC RF center (left eye)', ...
%         xEccentricityAxisName, yEccentricityAxisName);
%     
%     subplot(3,4,8);
%     makeLinePlot(eccXYposDegs(nasalIndices,1), nasalRatios, ...
%                  eccXYposDegs(temporalIndices,1), temporalRatios, ...
%                  quadrantColors, {'nasal' 'temporal'}, ...
%                  maxEcc, -30:5:30, [1 20], 1:1:20, xEccentricityAxisName, 'cones / mRGC RF center');
%              
%     subplot(3,4,12);
%     makeLinePlot(eccXYposDegs(superiorIndices,2), superiorRatios, ...
%                  eccXYposDegs(inferiorIndices,2), inferiorRatios, ...
%                  quadrantColors, {'superior' 'inferior'}, ...
%                  maxEcc, -30:5:30, [1 20], 1:1:20, yEccentricityAxisName, 'cones / mRGC RF center');
%              
%              
%     subplot(3,4,7);
%     makeLinePlot(eccXYposDegs(nasalIndices,1), nasalRatios, ...
%                  eccXYposDegs(temporalIndices,1), temporalRatios, ...
%                  quadrantColors, {'nasal' 'temporal'}, ...
%                  5, -5:1:5, [1 5.0], 1:0.5:5.0, xEccentricityAxisName, 'cones / mRGC RF center');
%              
%     subplot(3,4,11);
%     makeLinePlot(eccXYposDegs(superiorIndices,2), superiorRatios, ...
%                  eccXYposDegs(inferiorIndices,2), inferiorRatios, ...
%                  quadrantColors, {'superior' 'inferior'}, ...
%                  5, -5:1:5, [1 5.0], 1:0.5:5.0, yEccentricityAxisName, 'cones / mRGC RF center');
%              
%              
% end


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

function makeLinePlot(eccDegs1, ratios1, eccDegs2, ratios2, ...
                 color1,color2, dataNames, ...
                 maxEcc, eccTicks, yRange, yTicks, xLabelName, yLabelName)
    
    plot(eccDegs1, ratios1, ...
        'ko-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'Color', color1(1,:), ...
        'MarkerFaceColor', color1(2,:), ...
        'MarkerEdgeColor', color1(1,:));
    hold on;
    
    plot(eccDegs2, ratios2, ...
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


function [eccXYposDegs, quadrants] = sampleRetinalSpace(WatsonRGCCalc, samplesPerMeridian, maxEcc, spacing)
    
    if (strcmp(spacing, 'log'))
        eccDegs = 0;
        minimalConeSpacingDegs = WatsonRGCCalc.coneRFSpacingAndDensity(eccDegs, 'superior meridian', 'Cones per deg2');
        xPos = logspace(log10(minimalConeSpacingDegs), log10(maxEcc), samplesPerMeridian);
    else
        xPos = linspace(0, maxEcc, samplesPerMeridian);
    end
    
    xPos = [-fliplr(xPos) xPos]; yPos = xPos;
    [X,Y] = meshgrid(xPos, yPos);
    eccXYposDegs = [X(:) Y(:)];
   
    quadrants = containers.Map();
    
    % Points along the nasal meridian
    quadrants('nasal') = struct(...
        'meridianIndices', find((eccXYposDegs(:,1) >= 0) & (eccXYposDegs(:,2) == min(abs(yPos)))), ...
        'colors', [0.0 0.7 0.0; 0.0 1.00 0.0]);
    
    % Points along the temporal meridian
    quadrants('temporal') = struct(...
        'meridianIndices', find((eccXYposDegs(:,1) <= 0) & (eccXYposDegs(:,2) == min(abs(yPos)))), ...
        'colors', [1.0 0.0 0.0; 1.0 0.50 0.50]);

    % Points along the superior meridian
    quadrants('superior') = struct(...
        'meridianIndices',  find((eccXYposDegs(:,2) >= 0) & (eccXYposDegs(:,1) == min(abs(xPos)))), ...
        'colors', [0.0 0.0 1.0; 0.5 0.5 1.0]);
    
    % Points along the inferior meridian
    quadrants('inferior') = struct(...
        'meridianIndices', find((eccXYposDegs(:,2) <= 0) & (eccXYposDegs(:,1) == min(abs(xPos)))), ...
        'colors', [0.0 0.0 0.0; 0.5 0.5 0.5]);
    
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

