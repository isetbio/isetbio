function visualizeFittedVisualSTFsOfComputeReadyMidgetRGCMosaic(...
    computeReadyMosaicFilename, mRGCMosaicSTFresponsesFilename)

    % Load the compute-ready MRGC mosaic
    load(computeReadyMosaicFilename, 'theComputeReadyMRGCmosaic');

    % Load the computed mRGC  STF responses 
    load(mRGCMosaicSTFresponsesFilename, ...
        'theMRGCMosaicOptimalSTFs', ...
        'visualRcDegsEstimates', ...
        'theMRGCMosaicSTFresponses', ...
        'theMRGCresponseTemporalSupportSeconds', ...
        'orientationsTested', 'spatialFrequenciesTested', ...
        'spatialPhasesDegs', 'coneContrasts');

    
    theSTFdata = theMRGCMosaicOptimalSTFs{1};
    idxRcDegs = find(strcmp(theSTFdata.DoGfitParams.names, 'RcDegs'));
    idxRsRcRatio = find(strcmp(theSTFdata.DoGfitParams.names, 'RsToRc'));
    idxKsKcRatio = find(strcmp(theSTFdata.DoGfitParams.names, 'kS/kC'));
    
    cellEccDegs = zeros(1,theComputeReadyMRGCmosaic.rgcsNum);
    cellRcDegs = cellEccDegs;
    cellRsRcRatios = cellEccDegs;
    cellKsKcRatios = cellEccDegs;
    cellSCintSensRatios = cellEccDegs;

    for iRGC = 1:theComputeReadyMRGCmosaic.rgcsNum

        theSTFdata = theMRGCMosaicOptimalSTFs{iRGC};

%         theSTFdata.measured
%         theSTFdata.DoGfit
%         theSTFdata.DoGfitParams

        cellEccDegs(iRGC) = sqrt(sum((theComputeReadyMRGCmosaic.rgcRFpositionsDegs(iRGC,:)).^2,2));
        cellRcDegs(iRGC) = theSTFdata.DoGfitParams.finalValues(idxRcDegs);
        cellRsRcRatios(iRGC) = theSTFdata.DoGfitParams.finalValues(idxRsRcRatio);
        cellKsKcRatios(iRGC) = theSTFdata.DoGfitParams.finalValues(idxKsKcRatio);
        
    end
    cellSCintSensRatios = cellKsKcRatios .* cellRsRcRatios.^2;


    [CK95RcDegs.eccentricityDegs, CK95RcDegs.data] = ...
        RGCmodels.CronerKaplan.digitizedData.parvoCenterRadiusAgainstEccentricity();

    [CK95RsRcRatios.eccentricityDegs, CK95RsRcRatios.data] = ...
        RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();

    [CK95KsKcRatios.eccentricityDegs, CK95KsKcRatios.data] = ...
        RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterPeakSensisitivityRatioAgainstEccentricity();

    [CK95SCintSensRatios.eccentricityDegs, CK95SCintSensRatios.data] = ...
       RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();


    figure(1); clf;
    subplot(2,4,1);
    plot(CK95RcDegs.eccentricityDegs, CK95RcDegs.data, 's'); hold on
    plot(cellEccDegs, cellRcDegs, 'ro');
    set(gca, 'XLim', [0.03 30], 'XScale', 'log', 'XTick', [0.03 0.1 0.3 1 3 10 30])
    ylabel('Rc (degs)');
    xlabel('ecc (degs)');
    axis('square');

    subplot(2,4,2);
    plot(CK95RsRcRatios.eccentricityDegs, 1./CK95RsRcRatios.data, 's'); hold on
    plot(cellEccDegs, cellRsRcRatios, 'ro');
    set(gca, 'YLim', [0 20]);
    set(gca, 'XLim', [0.03 30], 'XScale', 'log', 'XTick', [0.03 0.1 0.3 1 3 10 30])
    ylabel('Rs/Rc');
    xlabel('ecc (degs)');
    axis('square');
    
    subplot(2,4,6);
    [counts,edges] = histcounts(1./CK95RsRcRatios.data,0:0.5:20);
    countsPercentage = counts / numel(CK95RsRcRatios.data)*100;
    bar(edges(1:end-1), countsPercentage, 1, 'FaceColor','none','EdgeColor',[0 0 0], 'LineWidth', 1.0);
    hold 'on'
    [counts,edges] = histcounts(cellRsRcRatios,0:0.5:20);
    countsPercentage = counts / numel(cellRsRcRatios)*100;
    bar(edges(1:end-1), countsPercentage, 1, 'FaceAlpha', 0.5', 'FaceColor',[1 0.5 0.50],'EdgeColor',[1 0 0], 'LineWidth', 1.0);
    plot(median(1./CK95RsRcRatios.data)*[1 1], [0 20], 'k--', 'LineWidth', 1.5);

    subplot(2,4,3);
    plot(CK95KsKcRatios.eccentricityDegs, CK95KsKcRatios.data, 's'); hold on
    plot(cellEccDegs, cellKsKcRatios, 'ro');
    set(gca, 'YLim', [1e-4 1], 'YScale', 'log', 'YTick', [1e-4 1e-3 1e-2 1e-1 1]);
    set(gca, 'XLim', [0.03 30], 'XScale', 'log', 'XTick', [0.03 0.1 0.3 1 3 10 30])
    ylabel('Ks/Kc');
    xlabel('ecc (degs)');
    axis('square');

    subplot(2,4,4);
    plot(CK95SCintSensRatios.eccentricityDegs, CK95SCintSensRatios.data, 's'); hold on
    plot(cellEccDegs, cellSCintSensRatios, 'ro');
    set(gca, 'YLim', [0 1]);
    set(gca, 'XLim', [0.03 30], 'XScale', 'log', 'XTick', [0.03 0.1 0.3 1 3 10 30])

    ylabel('S/C int sens');
    xlabel('ecc (degs)');
    axis('square');

    subplot(2,4,8);
    [counts,edges] = histcounts(CK95SCintSensRatios.data,0:0.05:1.0);
    countsPercentage = counts / numel(CK95SCintSensRatios.data)*100;
    bar(edges(1:end-1), countsPercentage, 1, 'FaceColor','none','EdgeColor',[0 0 0], 'LineWidth', 1.0);
    hold 'on'
    [counts,edges] = histcounts(cellSCintSensRatios,0:0.05:1.0);
    countsPercentage = counts / numel(cellSCintSensRatios)*100;
    bar(edges(1:end-1), countsPercentage, 1, 'FaceAlpha', 0.5', 'FaceColor',[1 0.5 0.50],'EdgeColor',[1 0 0], 'LineWidth', 1.0);
    plot(median(CK95SCintSensRatios.data)*[1 1], [0 1], 'k--', 'LineWidth', 1.5);
    set(gca, 'XLim', [0 1], 'XTick', 0:0.1:1);

end

