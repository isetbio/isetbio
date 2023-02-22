function [hFigRcDegs, hFigRsRcRatios, hFigSCintSensRatios] = renderSTFfitPlots(hFigRcDegs, hFigRsRcRatios, hFigSCintSensRatios, ...
    theMidgetRGCMosaic, theMeridianAngles, iMeridianAngle, theMeridianFits)

   
    if (iMeridianAngle == 1)
        hFigRcDegs = figure(1); clf;
        set(hFigRcDegs, 'Position', [10 10 3000 1250], 'Color', [1 1 1]);
        drawnow;

        hFigRsRcRatios = figure(2); clf;
        set(hFigRsRcRatios, 'Position', [100 100 3000 1250], 'Color', [1 1 1]);
        drawnow;

        hFigSCintSensRatios = figure(3); clf;
        set(hFigSCintSensRatios, 'Position', [200 200 3000 1250], 'Color', [1 1 1]);
        drawnow;
    end

    % Retrieve the C&K data
    [CronerKaplanRcDegs.temporalEccDegs, CronerKaplanRcDegs.val] = ...
            RGCmodels.CronerKaplan.digitizedData.parvoCenterRadiusAgainstEccentricity();

    [CronerKaplanRcRsRatios.temporalEccDegs, CronerKaplanRcRsRatios.val] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();

    [CronerKaplanSCintSensRatios.temporalEccDegs, CronerKaplanSCintSensRatios.val] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();

    [CronerKaplanKsKcRatios.temporalEccDegs, CronerKaplanKsKcRatios.val] = ...
            RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterPeakSensisitivityRatioAgainstEccentricity();


    RcDegsLims = [0.5 1]; RcDegsTicks = 0.5:0.1:1.0;
    RsRcLims = [0 16]; RsRcTicks = 0:2:40;
    intSensitivitySCLims = [0 1]; intSensitivitySCTicks = 0:0.1:1.0;

    eccentricityLims = 3*[-1 1];
    eccentricityTicks = -5:1:5;

    % Retrieve RGC positions for RGCs along this meridian
    rgcIndicesAlongThisMeridian = theMeridianFits.rgcIndicesAlongThisMeridian;
    rgcRFpositionsDegs = theMidgetRGCMosaic.rgcRFpositionsDegs(rgcIndicesAlongThisMeridian,:);

    majorityCenterConeType = zeros(1, numel(rgcIndicesAlongThisMeridian));
    parfor i=1:numel(rgcIndicesAlongThisMeridian)
         [~, ~,  majorityCenterConeType(i)]= theMidgetRGCMosaic.centerConeTypeWeights(rgcIndicesAlongThisMeridian(i));
    end

    % Compute signed distances
    radialDistances = sqrt(sum(rgcRFpositionsDegs.^2,2));
    theMeridianAngle = theMeridianAngles(iMeridianAngle);
    if (theMeridianAngle == 90) || (theMeridianAngle == 270)
       signedDistances = radialDistances  .* sign(rgcRFpositionsDegs(:,2));
    else
       signedDistances = radialDistances  .* sign(rgcRFpositionsDegs(:,1));
    end


    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 3, ...
       'colsNum', numel(theMeridianAngles), ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);


     fittedParams = theMeridianFits.fittedParams;
     theFittedSTFs = theMeridianFits.fittedSTFs;


    % The Rc degs values
    figure(hFigRcDegs);       
    ax1 = subplot('Position', subplotPosVectors(1,iMeridianAngle).v);
    ax2 = subplot('Position', subplotPosVectors(2,iMeridianAngle).v);
    ax3 = subplot('Position', subplotPosVectors(3,iMeridianAngle).v);

    plotRGCpositionsAlongMeridian(ax1, ...
            rgcRFpositionsDegs, ...
            eccentricityLims, eccentricityTicks, ...
            sprintf('%d RGCs along the\n%2.1f deg meridian', size(rgcRFpositionsDegs,1), theMeridianAngle));


    plotDataAlongMeridian(ax2, ax3, ...
            signedDistances, fittedParams.achievedRcDegs*60, ...
            eccentricityLims, eccentricityTicks, ...
            RcDegsLims, RcDegsTicks, ...
            [RcDegsLims(1) 8], 0.5:0.5:10, ...
            'RcDegs (arc min)', sprintf('%d RGCs along %2.1f degs', size(rgcRFpositionsDegs,1), theMeridianAngle), ...
            majorityCenterConeType, ...
            fittedParams.conesNumPooledByTheRFcenter, ...
            fittedParams.temporalEquivalentEccDegs, ...
            CronerKaplanRcDegs.temporalEccDegs, CronerKaplanRcDegs.val*60 ...
    );


    % The Rs/Rc ratio values
    figure(hFigRsRcRatios);
    ax1 = subplot('Position', subplotPosVectors(1,iMeridianAngle).v);
    ax2 = subplot('Position', subplotPosVectors(2,iMeridianAngle).v);
    ax3 = subplot('Position', subplotPosVectors(3,iMeridianAngle).v);

    plotRGCpositionsAlongMeridian(ax1, ...
            rgcRFpositionsDegs, ...
            eccentricityLims, eccentricityTicks, ...
            sprintf('%d RGCs along the\n%2.1f deg meridian', size(rgcRFpositionsDegs,1), theMeridianAngle));

        
    plotDataAlongMeridian(ax2, ax3, ...
            signedDistances, fittedParams.achievedRsToRcRatios, ...
            eccentricityLims, eccentricityTicks, ...
            RsRcLims, RsRcTicks, ...
            [RsRcLims(1) 32], RsRcTicks, ...
            'Rs/Rc ratio', sprintf('%d RGCs along %2.1f degs', size(rgcRFpositionsDegs,1), theMeridianAngle), ...
            majorityCenterConeType, ...
            fittedParams.conesNumPooledByTheRFcenter, ...
            fittedParams.temporalEquivalentEccDegs, ...
            CronerKaplanRcRsRatios.temporalEccDegs, 1./CronerKaplanRcRsRatios.val...
            );


    % The S/C integrated sensitivity ratios
    figure(hFigSCintSensRatios);
    ax1 = subplot('Position', subplotPosVectors(1,iMeridianAngle).v);
    ax2 = subplot('Position', subplotPosVectors(2,iMeridianAngle).v);
    ax3 = subplot('Position', subplotPosVectors(3,iMeridianAngle).v);

    plotRGCpositionsAlongMeridian(ax1, ...
            rgcRFpositionsDegs, ...
            eccentricityLims, eccentricityTicks, ...
            sprintf('%d RGCs along the\n%2.1f deg meridian', size(rgcRFpositionsDegs,1), theMeridianAngle));

    plotDataAlongMeridian(ax2, ax3, ...
            signedDistances, fittedParams.achievedSCintSensRatios, ...
            eccentricityLims, eccentricityTicks, ...
            intSensitivitySCLims , intSensitivitySCTicks , ...
            intSensitivitySCLims , intSensitivitySCTicks , ...
            'S/C int. sens. ratio', sprintf('%d RGCs along\nthe %2.1f degs meridian', size(rgcRFpositionsDegs,1), theMeridianAngle), ...
            majorityCenterConeType, ...
            fittedParams.conesNumPooledByTheRFcenter, ...
            fittedParams.temporalEquivalentEccDegs, ...
            CronerKaplanSCintSensRatios.temporalEccDegs, CronerKaplanSCintSensRatios.val ...
            );

end


function plotRGCpositionsAlongMeridian(ax, rgcRFpositionsDegs, ...
    eccentricityLims, eccentricityTicks, plotTitle)

    plot(ax,rgcRFpositionsDegs(:,1),  rgcRFpositionsDegs(:,2), 'k.', ...
         'MarkerSize', 15)
    axis(ax, 'square');
    
    grid(ax, 'on');
    set(ax, 'XLim', eccentricityLims, 'YLim', eccentricityLims, ...
            'XTick', eccentricityTicks, 'YTick', eccentricityTicks, ...
            'FontSize', 16);
    title(ax, plotTitle);
    drawnow;
end


function plotDataAlongMeridian(ax1, ax2, ...
            signedDistances, fittedParams, ...
            XLims, XTicks, ...
            YLims, YTicks, ...
            YLimsCK, YTicksCK, ...
            yLabelString, titleString, ...
            majorityCenterConeTypes, ...
            conesNumPooledByTheRFcenters, ...
            temporalEquivalentEccDegs, ...
            CronerKaplanTemporalEccDegs, CronerKaplanValues)
 

        % Plot on ax1
        majorityCenterConeType = [cMosaic.LCONE_ID cMosaic.MCONE_ID nan];

        for i = 1:numel(majorityCenterConeType)
            if (isnan(majorityCenterConeType(i)))
                rgcIndicesToPlot = find(isnan(majorityCenterConeTypes));
                markerColor = [1 1 0];
            else
                switch (majorityCenterConeType(i))
                    case cMosaic.LCONE_ID
                        markerColor = [1 0 0];
                    case cMosaic.MCONE_ID
                        markerColor = [0 0.7 0];
                end
                rgcIndicesToPlot = find(majorityCenterConeTypes == majorityCenterConeType(i));
            end
           
            
            plot(ax1, signedDistances(rgcIndicesToPlot), fittedParams(rgcIndicesToPlot), 'o', 'MarkerSize', 7, ...
                'MarkerFaceColor', markerColor, 'MarkerEdgeColor', [0 0 0]);
            hold(ax1, 'on');
        end
        
        axis(ax1, 'square');
        set(ax1, 'XLim', XLims, 'XTick', XTicks, 'YLim', YLims, 'YTick', YTicks, 'FontSize', 16);
        grid(ax1, 'on');
        xtickangle(ax1, 0);
        ylabel(ax1, yLabelString);
        xlabel(ax1, 'eccentricity (degs)');
        title(ax1, titleString);

        % Now replot on ax2, along with the Croner & Kaplan data
        for i = 1:numel(majorityCenterConeType)
            if (isnan(majorityCenterConeType(i)))
                rgcIndicesToPlot = find(isnan(majorityCenterConeTypes));
                markerColor = [1 1 0];
            else
                switch (majorityCenterConeType(i))
                    case cMosaic.LCONE_ID
                        markerColor = [1 0 0];
                    case cMosaic.MCONE_ID
                        markerColor = [0 0.7 0];
                end
                rgcIndicesToPlot = find(majorityCenterConeTypes == majorityCenterConeType(i));
            end
            
            plot(ax2, temporalEquivalentEccDegs(rgcIndicesToPlot), fittedParams(rgcIndicesToPlot), 'o', 'MarkerSize', 7, ...
                'MarkerFaceColor', markerColor, 'MarkerEdgeColor', [0 0 0]);
            
            hold(ax2, 'on');
        end

        scatter(ax2, CronerKaplanTemporalEccDegs, CronerKaplanValues, 121, 'o', ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0.3 0.3 0.3], 'MarkerFaceColor', [0.6 0.6 0.6], 'LineWidth', 1.0);
        axis(ax2, 'square');
        set(ax2, 'XLim', [0.01 30], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30], 'XScale', 'log');
        xtickangle(ax2, 0);
        set(ax2, 'YLim', YLimsCK, 'YTick', YTicksCK, 'FontSize', 16);
        grid(ax2, 'on');
        ylabel(ax2, yLabelString);
        xlabel(ax2, 'temporal equiv. eccentricity (degs)');
        drawnow

end