function computeConnectivityMatrix

    % Instantiate a WatsonRGCModel object.
    WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);
    
    % Sample the retinal space
    samplesPerMeridian = 50; maxEcc = 30; spacing = 'log';
    [eccXYposDegs, quadrants] = sampleRetinalSpace(WatsonRGCCalc,samplesPerMeridian, maxEcc, spacing);
    
    % Compute midgetRGCRF to cone ratios at the above eccentricities
    whichEye = 'left';
    [midgetRGCRFToConeRatios, coneRFDensitiesDeg2, mRGCRFDensitiesDeg2] ...
        = WatsonRGCCalc.ratioOfMidgetRGCsToCones(eccXYposDegs, whichEye);
    
    % Compute cones per ON or OFF midget RGC
    conesPerMidgetRGCRF = 1 ./ (0.5*midgetRGCRFToConeRatios);
    
    % Compute midget RGC RF spacing from the their density using equation (9) in the Watson (2014) paper.
    singlePolarityMidgetRGCRFDensitiesDeg2 = 0.5*mRGCRFDensitiesDeg2; % only ON or OFF
    singlePolarityMidgetRGCRFSpacingsDeg = sqrt(2.0./(sqrt(3.0)*singlePolarityMidgetRGCRFDensitiesDeg2));
     
    % Compute cone spacing from their density using equation (9) in the Watson (2014) paper.
    coneSpacingsDeg = sqrt(2.0./(sqrt(3.0)*coneRFDensitiesDeg2));
    
    
    % Meridian axes names
    xEccentricityAxisName = sprintf('space (deg)\n   <- nasal  |  temporal ->');
    yEccentricityAxisName = sprintf('space (deg)\n<- inferior  |  superior ->');
    
    % Plot cones per mRGC along meridians
    figureNo = 1;
    extraData = [];
    renderMeridianStats(figureNo, conesPerMidgetRGCRF,eccXYposDegs, quadrants, ...
        extraData, maxEcc*[-1 1], -30:5:30, [1 20], 1:1:20, ...
        xEccentricityAxisName, yEccentricityAxisName, 'cones / mRGC RF center');
    
    
    % Load Bradley-Geisler data
    load('BradleyGeislerRGCspacing.mat','RGCspacingBradleyGeisler');
    
    % Plot mRGC RF spacing along meridians
    figureNo = 2;
    extraData = RGCspacingBradleyGeisler;
    renderMeridianStats(figureNo, singlePolarityMidgetRGCRFSpacingsDeg*60, eccXYposDegs, quadrants, ...
        extraData, 11*[-1 1], -10:1:10, [0 12], 0:1:12, ...
        xEccentricityAxisName, yEccentricityAxisName, ...
        'ON or OFF mRGC RF spacing (arc min)');
    
    % Plot cone spacing along meridians
    figureNo = 3;
    extraData = [];
    renderMeridianStats(figureNo, coneSpacingsDeg*60, eccXYposDegs, quadrants, ...
        extraData, 11*[-1 1], -30:1:30, [0 12], 0:1:12, ...
        xEccentricityAxisName, yEccentricityAxisName, ...
        'cone spacing (arc min)');
    
    % Plot contour map of midgetRGCRF to cones ratio
    figureNo = 4;
    hFig = figure(figureNo);
    clf;
     
    eccRange = [-30 30];
    eccTicks = eccRange(1):10:eccRange(2);

    subplot(2,2,1);
    spacingArcMinLevels(1) = min(coneSpacingsDeg(:))*60;
    spacingArcMinLevels = cat(2, spacingArcMinLevels, (1:0.5:10));
    spacingArcMinLevelsLabeled = [spacingArcMinLevels(1) spacingArcMinLevels(2:2:end)];
    spacingCMap = brewermap(numel(spacingArcMinLevels), 'YlGnBu');
    renderContourPlot(coneSpacingsDeg*60, eccXYposDegs, samplesPerMeridian, ...
         eccRange, eccTicks, spacingArcMinLevels, spacingArcMinLevelsLabeled, spacingCMap, ...
         'cone spacing (arc min)', ...
         xEccentricityAxisName, yEccentricityAxisName);
     
    subplot(2,2,3);
    spacingArcMinLevels = [];
    spacingArcMinLevels(1) = min(singlePolarityMidgetRGCRFSpacingsDeg(:))*60;
    spacingArcMinLevels = cat(2, spacingArcMinLevels, (1:0.5:10));
    spacingArcMinLevelsLabeled = [spacingArcMinLevels(1) spacingArcMinLevels(2:2:end)];
    renderContourPlot(singlePolarityMidgetRGCRFSpacingsDeg*60, eccXYposDegs, samplesPerMeridian, ...
        eccRange, eccTicks, spacingArcMinLevels, spacingArcMinLevelsLabeled, spacingCMap, ...
         'mRGC RF spacing (arc min)', ...
         xEccentricityAxisName, yEccentricityAxisName);
     

    subplot(2,2,2);
    conesPerMidgetRGCLevels = 1:2:51;
    conesPerMidgetRGCLevelsLabeled = [1 5 10 15 20 25 30 35 40 45 50];
    cMap = brewermap(numel(conesPerMidgetRGCLevels), '*YlGnBu');
    renderContourPlot(conesPerMidgetRGCRF, eccXYposDegs, samplesPerMeridian, ...
        eccRange, eccTicks, conesPerMidgetRGCLevels, conesPerMidgetRGCLevelsLabeled, cMap, ...
        'cones per mRGC', ...
        xEccentricityAxisName, yEccentricityAxisName);
    
    subplot(2,2,4);
    eccRange = [-10 10];
    eccTicks = -10:1:10; 
    conesPerMidgetRGCLevels = 1:0.5:6;
    conesPerMidgetRGCLevelsLabeled = conesPerMidgetRGCLevels(1:2:end);
    cMap = brewermap(numel(conesPerMidgetRGCLevels), '*YlGnBu');
    renderContourPlot(conesPerMidgetRGCRF, eccXYposDegs, samplesPerMeridian, ...
        eccRange, eccTicks, conesPerMidgetRGCLevels, conesPerMidgetRGCLevelsLabeled, cMap, ...
        'cones per mRGC', ...
        xEccentricityAxisName, yEccentricityAxisName);
     
end

function renderMeridianStats(figureNo, stats, eccXYposDegs, quadrants, ...
    extraData, eccRange, eccTicks, statsRange, statsTicks, ...
    xEccentricityAxisName, yEccentricityAxisName, statsName)

    nasalStats = stats(quadrants('nasal').meridianIndices);
    nasalEccDegs = eccXYposDegs(quadrants('nasal').meridianIndices,1);
    
    temporalStats = stats(quadrants('temporal').meridianIndices);
    temporalEccDegs = eccXYposDegs(quadrants('temporal').meridianIndices,1);
    
    superiorStats = stats(quadrants('superior').meridianIndices);
    superiorEccDegs = eccXYposDegs(quadrants('superior').meridianIndices,2);
    
    inferiorStats = stats(quadrants('inferior').meridianIndices);
    inferiorEccDegs = eccXYposDegs(quadrants('inferior').meridianIndices,2);
    
    if (~isempty(extraData))
        % visual field is complementary to retinal field
        % superior retina == inferir visual field
        extraNasalEccDegs = extraData('temporalField').ecc;
        extraNasalStats = extraData('temporalField').spacingDegs*60;
        extraTemporalEccDegs = -extraData('nasalField').ecc;
        extraTemporalStats = extraData('nasalField').spacingDegs*60;
        extraSuperiorEccDegs = extraData('inferiorField').ecc;
        extraSuperiorStats = extraData('inferiorField').spacingDegs*60;
        extraInferiorEccDegs = -extraData('superiorField').ecc;
        extraInferiorStats = extraData('superiorField').spacingDegs*60;
    else
        extraNasalEccDegs = [];
        extraNasalStats = [];
        extraTemporalEccDegs = [];
        extraTemporalStats = [];
        extraSuperiorEccDegs = [];
        extraSuperiorStats = [];
        extraInferiorEccDegs = [];
        extraInferiorStats = [];
    end
    
    hFig = figure(figureNo); clf;
    set(hFig, 'Position', [rand(1,1)*100 rand(1,1)*100 850 450]);
    subplot(1,2,1);
    renderLinesPlot(nasalEccDegs, nasalStats, ...
                 temporalEccDegs, temporalStats, ...
                 extraNasalEccDegs, extraNasalStats, ...
                 extraTemporalEccDegs, extraTemporalStats, ...
                 quadrants('nasal').colors, quadrants('temporal').colors, ...
                 {'nasal', 'temporal'}, ...
                 eccRange, eccTicks, statsRange, statsTicks, xEccentricityAxisName, statsName);
             
    subplot(1,2,2);
    renderLinesPlot(superiorEccDegs, superiorStats, ...
                 inferiorEccDegs, inferiorStats, ...
                 extraSuperiorEccDegs, extraSuperiorStats, ...
                 extraInferiorEccDegs, extraInferiorStats, ...
                 quadrants('superior').colors, quadrants('inferior').colors, ...
                 {'superior' 'inferior'}, ...
                 eccRange, eccTicks, statsRange, statsTicks, yEccentricityAxisName, statsName);
             
             
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


function renderContourPlot(stats, eccXYposDegs, samplesPerMeridian, eccRange, eccTicks, ...
    zLevels, zLevelsOnCbar, cMap, titleName, xLabelName, yLabelName)
    
    X = eccXYposDegs(:,1);
    Y = eccXYposDegs(:,2);
    N = 2*samplesPerMeridian;
    contourf(reshape(X, [N N]), reshape(Y, [N N]), reshape(stats, [N N]), zLevels);
    hold on;
    plot([0 0], eccRange, 'k-', 'LineWidth', 1.0);
    plot(eccRange, [0 0], 'k-', 'LineWidth', 1.0);
    axis 'image';
    set(gca, 'XLim', eccRange, 'YLim', eccRange, 'CLim', [min(zLevels) max(zLevels)], ...
        'XTick', eccTicks, 'YTick', eccTicks, 'FontSize', 14);
    
    colormap(gca,cMap)
    hL = colorbar;
    set(hL, 'Ticks', zLevelsOnCbar, 'Limits', [min(zLevels) max(zLevels)]);
    title(titleName);
    xlabel(xLabelName, 'FontName', 'Menlo', 'FontSize', 12, 'FontAngle', 'italic');
    ylabel(yLabelName, 'FontName', 'Menlo', 'FontSize', 12, 'FontAngle', 'italic');
end

function renderLinesPlot(x1, y1, x2, y2, extraX1, extraY1, extraX2, extraY2, color1, color2, dataNames, ...
                 xRange, xTicks, yRange, yTicks, xLabelName, yLabelName)
    
    plot(x1, y1, ...
        'ko-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'Color', color1(1,:), ...
        'MarkerFaceColor', color1(2,:), ...
        'MarkerEdgeColor', color1(1,:));
    hold on;
    
    plot(x2, y2, ...
        'ko-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'Color', color2(1,:), ...
        'MarkerFaceColor', color2(2,:), ...
        'MarkerEdgeColor', color2(1,:));
    
    if (~isempty(extraX1))
        plot(extraX1, extraY1, ...
            'ko--', 'LineWidth', 1.5, 'MarkerSize', 6, ...
            'Color', color1(1,:), ...
            'MarkerFaceColor', color1(2,:), ...
            'MarkerEdgeColor', color1(1,:));
    end
    if (~isempty(extraX2))
        plot(extraX2, extraY2, ...
            'ko-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
            'Color', color2(1,:), ...
            'MarkerFaceColor', color2(2,:), ...
            'MarkerEdgeColor', color2(1,:));
    end
    
    xlabel(xLabelName, 'FontName', 'Menlo', 'FontAngle', 'italic');
    ylabel(yLabelName, 'FontName', 'Menlo', 'FontAngle', 'italic');
    legend(dataNames, 'Location', 'North');
    
    set(gca, 'XLim', xRange, 'YLim', yRange, ...
         'XScale', 'linear', 'YScale', 'linear', ...
         'XTick', xTicks, 'YTick', yTicks, ...
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

