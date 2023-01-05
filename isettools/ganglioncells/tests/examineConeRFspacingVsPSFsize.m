function examineConeRFspacingVsPSFsize()

    retinaQuadrant = 'nasal meridian';
    if (strcmp(retinaQuadrant, 'nasal meridian'))
        radialEccExamined = [0 1 2 3 4 6 8 12 19 24 30];
    else
        radialEccExamined = [0 1 2 3 4 6 8 12 16 20 24 30];
    end

    % Optics subject
    opticsDataBase = 'Artal2012';

    includeConeApertureInPSF = false;

    for iEcc = 1:numel(radialEccExamined)
        if (strcmp(retinaQuadrant, 'nasal meridian')) && (radialEccExamined(iEcc) == 0)
            [psfData{iEcc}, coneData{iEcc}] = retrievePSFandConeData(radialEccExamined(iEcc), opticsDataBase, 'temporal meridian',  includeConeApertureInPSF);
        else
            [psfData{iEcc}, coneData{iEcc}] = retrievePSFandConeData(radialEccExamined(iEcc), opticsDataBase, retinaQuadrant,  includeConeApertureInPSF);
        end

    end

    [radialEccExaminedNasal, psfRcMinNasal, psfRcMaxNasal, coneSpacingNasal] = plotSummaryData(radialEccExamined, retinaQuadrant, psfData, coneData);
    
    hFig = figure(3); clf;
    
    plot(radialEccExaminedNasal, psfRcMinNasal, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [1 0.5 0.5]);
    hold on
    plot(radialEccExaminedNasal, coneSpacingNasal, 'ro--', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [1 0.5 0.5]);
    set(gca, 'FontSize', 14);
    

    hFig = figure(4); clf;
    set(hFig, 'Position', [10 10 330 420], 'Color', [1 1 1]);
    plot([0 10], [0 10], 'k-', 'LineWidth', 1.0);
    hold on
    p1 = plot(coneSpacingNasal, psfRcMinNasal, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [1 0.75 0.75]);

    p2 = plot(coneSpacingNasal, psfRcMaxNasal, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [1 0.5 0.5]);
    set(gca, 'XLim', [0 7], 'YLim', [0 7]);
    set(gca, 'XTick', 0:7, 'YTick', 0:7, 'FontSize', 14);
    grid on
    xlabel('cone spacing (arc min)')
    ylabel('psf characteristic radius (arc min)');


    retinaQuadrant = 'temporal meridian';
    if (strcmp(retinaQuadrant, 'nasal meridian'))
        radialEccExamined = [0 1 2 3 4 6 8 12 19 24 30];
    else
        radialEccExamined = [0 1 2 3 4 6 8 12 16 20 24 30];
    end
    for iEcc = 1:numel(radialEccExamined)
        if (strcmp(retinaQuadrant, 'nasal meridian')) && (radialEccExamined(iEcc) == 0)
            [psfData{iEcc}, coneData{iEcc}] = retrievePSFandConeData(radialEccExamined(iEcc), opticsDataBase, 'temporal meridian',  includeConeApertureInPSF);
        else
            [psfData{iEcc}, coneData{iEcc}] = retrievePSFandConeData(radialEccExamined(iEcc), opticsDataBase, retinaQuadrant,  includeConeApertureInPSF);
        end
    end
    [radialEccExaminedTemporal, psfRcMinTemporal, psfRcMaxTemporal, coneSpacingTemporal] = plotSummaryData(radialEccExamined, retinaQuadrant, psfData, coneData);
   
    figure(4); 
    p3 = plot(coneSpacingTemporal, psfRcMinTemporal, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [0.75 0.75 1]);
    hold on
    p4 = plot(coneSpacingTemporal, psfRcMaxTemporal, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0.5 1]);
    axis 'square';

    legend([p1,p2,p3,p4], {'Rc min (nasal)', 'Rc max (nasal)', 'Rc min (temporal)', 'Rc max (temporal)'}, ...
        'Box', 'off', 'Location', 'NorthOutside', 'NumColumns', 2);

    NicePlot.exportFigToPDF('coneSpacingVsPSFRc.pdf', hFig, 300);

    plotHorizontalMeridianData(radialEccExaminedNasal, psfRcMinNasal, psfRcMaxNasal, coneSpacingNasal, ...
        radialEccExaminedTemporal, psfRcMinTemporal, psfRcMaxTemporal, coneSpacingTemporal);

end

function plotHorizontalMeridianData(radialEccExaminedNasal, psfRcMinNasal, psfRcMaxNasal, coneSpacingNasal, ...
        radialEccExaminedTemporal, psfRcMinTemporal, psfRcMaxTemporal, coneSpacingTemporal)

    radialEcc         = cat(2, fliplr(radialEccExaminedNasal), radialEccExaminedTemporal);
    psfRcMinArcMin    = cat(2,fliplr(psfRcMinNasal), psfRcMinTemporal);
    psfRcMaxArcMin    = cat(2,fliplr(psfRcMaxNasal), psfRcMaxTemporal);
    coneSpacingArcMin = cat(1,flipud(coneSpacingNasal), coneSpacingTemporal);

    hFig = figure(1); clf
    set(hFig, 'Position', [100 100 660 450], 'Color', [1 1 1]);
    ax = subplot('Position', [0.08 0.13 0.91 0.86]);

    includeCones = true;
    includePSF = false;

    if (includePSF)
        p1 = shadedAreaBetweenTwoLines(ax,radialEcc, ...
            psfRcMinArcMin, ...
            psfRcMaxArcMin, ...
            [0.2 1 0.8], [0.1 0.6 0.4], 0.5, 1.0, '-');
    end

    hold(ax, 'on');

    
    if (includeCones)
        p2 = plot(ax, radialEcc, coneSpacingArcMin, 'ro-', ...
        'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', [1 0.5 0.5]);
    end


    characteristicRadiusRange = [0 7];
    eccRange = [-30 30];
    set(ax, 'XLim', eccRange, 'YLim', characteristicRadiusRange, ...
        'XTick', -30:2:30, ...
        'XTickLabel', {'', '-28', '', '-24', '', '-20', '', '-16', '', '-12', '', '-8', '', '-4', '', '0', '', '4', '', '8', '', '12', '', '16', '', '20', '', '24', '', '28', ''}, ...
        'YTick', 0:0.5:7, 'YTickLabel', {'0', '', '1', '', '2', '', '3', '', '4', '', '5', '', '6', '', '7'}, 'FontSize', 18);
    grid(ax, 'on'); box(ax, 'off');
    set(ax, 'LineWidth', 1.0, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2], 'TickDir', 'both');
    xlabel(ax,'\leftarrow {\it nasal retina}         eccentricity (degs)          {\it temporal retina} \rightarrow');
    xtickangle(0);
    ylabel(ax,'arcmin');
    if (includeCones && includePSF)
        legend([p1, p2], {'PSF characteristic radius (55 subjects)', 'cone spacing'}, 'Location', 'NorthOutside', 'NumColumns', 2, 'box', 'off');
    elseif (~includeCones && includePSF)
        legend([p1], {'PSF characteristic radius (55 subjects)'}, 'Location', 'NorthOutside', 'NumColumns', 2, 'box', 'off');
    else
        legend([p2], {'cone spacing'}, 'Location', 'NorthOutside', 'NumColumns', 2, 'box', 'off');
    end
    
    NicePlot.exportFigToPDF('NasalTemporal.pdf', hFig, 300);
    

end

function [radialEccExamined, psfRcMin, psfRcMax, coneSpacing] = plotSummaryData(radialEccExamined, retinaQuadrant, psfData, coneData)

    if (strcmp(retinaQuadrant, 'nasal meridian'))
        radialEccExamined = -radialEccExamined;
    end

    exclusiveCenterConesNum = zeros(numel(radialEccExamined),1);
    midgetRGCRFcenterRadiiDegs = zeros(numel(radialEccExamined), 20, 2);
    opticalPSFRadiiDegs = zeros(numel(radialEccExamined),2);
    allSubjectsOpticalPSFRadiiDegs = zeros(numel(radialEccExamined),55,2);

    for iEcc = 1:numel(radialEccExamined)
        % PSF data from all subjects
        allSubjectsOpticalPSFRadiiDegs(iEcc,:,:) = psfData{iEcc}.allSubjectPSFcharacteristicRadiiDegs;
        localConeSpacings(iEcc,:) = coneData{iEcc}.IQRs;
    end

    eccRange = [0 max(abs(radialEccExamined))+0.5];
    eccTicks = 0:2:30;
    if (strcmp(retinaQuadrant, 'nasal meridian'))
        eccRange = [-eccRange(2) eccRange(1)];
        eccTicks = (-30):2:0;
    end

    characteristicRadiusRange = [0 8];
    markerSize = 10;

    allSubjectsMajorCharacteristicRadiiDegs = squeeze(allSubjectsOpticalPSFRadiiDegs(:,:,1));
    allSubjectsMinorCharacteristicRadiiDegs = squeeze(allSubjectsOpticalPSFRadiiDegs(:,:,2));

    if (strcmp(retinaQuadrant, 'nasal meridian'))
         hFig = figure(1); clf;
         set(hFig, 'Position', [10 10 330 370], 'Color', [1 1 1]);
    else
        hFig = figure(2); clf;
        set(hFig, 'Position', [1000 10 330 370], 'Color', [1 1 1]);
    end
    

    psfRcMin = median(allSubjectsMinorCharacteristicRadiiDegs,2)'*60;
    psfRcMax = median(allSubjectsMajorCharacteristicRadiiDegs,2)'*60;
    coneSpacing = localConeSpacings(:,1)*60;

    ax = subplot('Position', [0.1 0.1 0.89 0.89]);
    p1 = shadedAreaBetweenTwoLines(ax,radialEccExamined, ...
        psfRcMin, ...
        psfRcMax, ...
        [0.2 1 0.8], [0.1 0.6 0.4], 0.5, 1.0, '-');
    hold(ax, 'on');

    %plot(ax, radialEccExamined, localConeSpacings(:,1)*60, 'r--');
    p2 = plot(ax, radialEccExamined, coneSpacing, 'ro-', ...
        'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', [1 0.5 0.5]);
    %plot(ax, radialEccExamined, localConeSpacings(:,3)*60, 'r:');

    set(ax, 'XLim', eccRange, 'YLim', characteristicRadiusRange, ...
        'XTick', eccTicks, 'YTick', 0:1:8, 'FontSize', 14);
    axis(ax, 'square'); grid(ax, 'on'); box(ax, 'off');
    set(ax, 'LineWidth', 1.0, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2], 'FontSize', 14);
    if (strcmp(retinaQuadrant, 'nasal meridian'))
        xlabel(ax,'eccentricity, nasal retina (degs)');
    else
        xlabel(ax,'eccentricity, temporal retina (degs)');
    end
    if (strcmp(retinaQuadrant, 'nasal meridian'))
        ylabel(ax,'characteristic radius (arcmin)');
        legend([p1, p2], {'PSF Rc', 'cone spacing'}, 'Location', 'NorthWest');
    else
        set(ax, 'YTickLabel', {});
    end
    if (strcmp(retinaQuadrant, 'nasal meridian'))
        NicePlot.exportFigToPDF('nasal.pdf', hFig, 300);
    else
        NicePlot.exportFigToPDF('temporal.pdf', hFig, 300);
    end

end


function [psfData, coneData] = retrievePSFandConeData(radialEcc,opticsDataBase, retinaQuadrant, includeConeApertureInPSF)

    whichEye = 'right eye';
    [horizontalEcc, verticalEcc] = cMosaic.eccentricitiesForRetinaMeridianInEye(...
            radialEcc, retinaQuadrant, whichEye);
    eccDegs = [horizontalEcc verticalEcc];

    sizeDegs = max([0.25 0.25+max(abs(eccDegs))*0.1]);
    c = cMosaic(...
        'sourceLatticeSizeDegs', 60, ...
        'eccentricityDegs', eccDegs, ...
        'sizeDegs', sizeDegs*[1 1]);
    coneData.coneSpacingsDegs = c.coneRFspacingsDegs;
    coneData.IQRs = prctile(coneData.coneSpacingsDegs, [1 50 99]);

   
    if (includeConeApertureInPSF)
        fName = sprintf('allSubjectsPSFConeApertureAnalysis_%2.2fdegs_%2.2fdegs_%s', horizontalEcc, verticalEcc, retinaQuadrant);
    else
        fName = sprintf('allSubjectsPSFanalysis_%2.2fdegs_%2.2fdegs_%s', horizontalEcc, verticalEcc, retinaQuadrant);
    end

    %Load computed PSF radii 
    load(sprintf('%s.mat', fName), 'allSubjectPSFcharacteristicRadiiDegs', 'opticsDataBase', 'opticsSubjectRankOrders');


    psfData = struct();
    psfData.allSubjectPSFcharacteristicRadiiDegs = allSubjectPSFcharacteristicRadiiDegs;
    psfData.opticsDataBase = opticsDataBase;
    psfData.subjectRankOrders = opticsSubjectRankOrders;
 
end

function plotHandle = shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    plotHandle = patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end