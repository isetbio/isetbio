function compareConeExcitationsToPcurrentCronerKaplan()

    stfResponsesRootDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris/IBIO_rgcMosaicResources/ONcenterMidgetRGCmosaics/intermediateFiles/SLIM/STFresponses';
    analysesRootDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris/IBIO_rgcMosaicResources/ONcenterMidgetRGCmosaics/intermediateFiles/SLIM/demos/CronerKaplanAnalyses';

    stimulusContrast = 25;
    backgroundLuminanceMultiplier = 1;
    stimulusTemporalFrequency = 4;

    %stimulusContrast = 50;
    %backgroundLuminanceMultiplier = 2;
    %stimulusTemporalFrequency = 4;

    stimStringForPhotocurrents = sprintf('%2.0f_x%1.0f_%02.0fHz', stimulusContrast, backgroundLuminanceMultiplier, stimulusTemporalFrequency);

    eccDegs = -0.1; mosaicSpecifier = 'Ecc0.0_0.0_Size2.0x2.0';
    eccDegs = -6.0;   mosaicSpecifier = 'Ecc-7.0_0.0_Size4.0x4.0';

    eccDegs = -12.5;   mosaicSpecifier = 'Ecc-14.0_0.0_Size5.0x5.0';
    eccDegs = -20.0;   mosaicSpecifier = 'Ecc-19.0_0.0_Size6.0x6.0';


    doIt(1, stfResponsesRootDir, analysesRootDir, mosaicSpecifier, stimStringForPhotocurrents, eccDegs)

    if (1==2)
        eccDegs = -1.1; mosaicSpecifier = 'Ecc-2.0_0.0_Size2.0x2.0';
        doIt(2, stfResponsesRootDir, analysesRootDir, mosaicSpecifier, stimStringForPhotocurrents, eccDegs)

        eccDegs = -4.5; mosaicSpecifier = 'Ecc-4.0_0.0_Size3.0x3.0';
        doIt(3, stfResponsesRootDir, analysesRootDir, mosaicSpecifier, stimStringForPhotocurrents, eccDegs)

        eccDegs = -8;   mosaicSpecifier = 'Ecc-7.0_0.0_Size4.0x4.0';
        doIt(4, stfResponsesRootDir, analysesRootDir, mosaicSpecifier, stimStringForPhotocurrents, eccDegs)
    end

end

function doIt(figNo, stfResponsesRootDir, analysesRootDir, mosaicSpecifier, stimStringForPhotocurrents, eccDegs)

    if (abs(eccDegs) < 15)
        basicMosaicString = 'Phi_1.00_Optics_Polans2015-3_maxStrehlRatio_srndModel_PackerDacey2002H1freeLowH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0';
    else
        basicMosaicString = 'Phi_1.00_Optics_Polans2015-3_maxStrehlRatio_srndModel_PackerDacey2002H1freeMidH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0';
    end

    coneExcitationString = sprintf('MRGCMosaic_RE_%s_%s_mRGCMosaic_nativeOptics_Achromatic_@%2.1f_0.0.mat', ...
        mosaicSpecifier, basicMosaicString,  eccDegs);

    fullfile(stfResponsesRootDir, coneExcitationString)
    load(fullfile(stfResponsesRootDir, coneExcitationString), 'theMRGCMosaic');

    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = 'validation';
    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, sprintf('MRGCmosaicAt%2.1fdegs_Stim_%s.pdf', eccDegs, stimStringForPhotocurrents));

    visualizedFOVdegs = 1.8;
    domainVisualizationLimits(1:2) = theMRGCMosaic.eccentricityDegs(1) + [-1 1]*0.5*visualizedFOVdegs;
    domainVisualizationLimits(3:4) = theMRGCMosaic.eccentricityDegs(2) + [-1 1]*0.5*visualizedFOVdegs;
    domainVisualizationTicks.x = theMRGCMosaic.eccentricityDegs(1) + (-1:0.5:1)*0.5*visualizedFOVdegs;
    domainVisualizationTicks.y = theMRGCMosaic.eccentricityDegs(2) + (-1:0.5:1)*0.5*visualizedFOVdegs;

    hFig = figure(100+figNo);
    set(hFig, 'Color', [1 1 1], 'Position', [10+figNo*50, 10+figNo*50, 450 450]);
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);
    theMRGCMosaic.visualize(...
        'figureHandle', hFig, 'axesHandle', ax, ...
        'centerSubregionContourSamples', 40, ...
        'spatialSupportSamples', 1024, ...
        'identifyInputCones', true, ...
        'identifiedConeAperture','lightCollectingArea4sigma', ...
        'identifiedConeApertureThetaSamples', 32, ...
        'plottedRFoutlineLineWidth', 1.0, ...
        'plottedRFoutlineFaceColor', [0.1 0.6 0.1], ...
        'plottedRFoutlineFaceAlpha', 0.5, ...
        'contourGenerationMethod', 'ellipseFitToPooledConePositions', ... %'ellipseFitBasedOnLocalSpacing');  % 'contourOfPooledConeApertureImage'
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks);

    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300);

    [allRGCCenterConesNum, allRGCCenterDominantConeTypes, allRGCCenterRelativeConeWeights] = theMRGCMosaic.allRFcenterConnectivityStats();
    LcenterRGCindices = find(allRGCCenterDominantConeTypes == cMosaic.LCONE_ID);
    McenterRGCindices = find(allRGCCenterDominantConeTypes == cMosaic.MCONE_ID);


    load(fullfile(analysesRootDir, coneExcitationString));
    RcDegsExcitations = RcDegsMosaic;
    KcExcitations = KcMosaic;
    RsToRcExcitations = RsToRcMosaic;
    intStoCsensExcitations = intStoCsensMosaic;


    % Photocurrent data
    coneExcitationString = sprintf('MRGCMosaic_RE_%s_%s_mRGCMosaic_nativeOptics_Achromatic%s_@%2.1f_0.0.mat', ...
        mosaicSpecifier, basicMosaicString, stimStringForPhotocurrents, eccDegs);
    pCurrentString = strrep(coneExcitationString, 'mRGCMosaic', 'mRGCMosaicPhotocurrent');

    fullfile(analysesRootDir, pCurrentString)
    load(fullfile(analysesRootDir, pCurrentString));
    RcDegsPcurrent = RcDegsMosaic;
    KcPcurrent = KcMosaic;
    RsToRcPcurrent = RsToRcMosaic;
    intStoCsensPcurrent = intStoCsensMosaic;

    comparisonPlots(figNo, eccDegs, stimStringForPhotocurrents, ...
        LcenterRGCindices, McenterRGCindices, ...
        intStoCsensExcitations, intStoCsensPcurrent, ...
        RsToRcExcitations, RsToRcPcurrent, ...
        KcExcitations, KcPcurrent, ...
        RcDegsExcitations, RcDegsPcurrent);

end

function comparisonPlots(figNo, eccDegs, stimStringForPhotocurrents, ...
    LcenterRGCindices, McenterRGCindices, ...
    intStoCsensExcitations, intStoCsensPcurrent, ...
    RsToRcExcitations, RsToRcPcurrent, ...
    KcExcitations, KcPcurrent, ...
    RcDegsExcitations, RcDegsPcurrent)


    % Do not plot outliers
    KcPercentile = 100;
    KcPcurrent = KcPcurrent/(prctile(KcPcurrent(:),KcPercentile)/max(KcExcitations(:)));
    maxKc = max(KcExcitations(:));
    idx = find(KcPcurrent<=maxKc);

    LcenterRGCindices = intersect(LcenterRGCindices, idx);
    McenterRGCindices = intersect(McenterRGCindices, idx);

    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10+figNo*50, 10+figNo*50, 1150 950]);
    ax = subplot(2,2,1);
    generateCorrespondencePlot(ax, intStoCsensExcitations, intStoCsensPcurrent, ...
        LcenterRGCindices, McenterRGCindices, 'intS / intC', [0 1.4], 0.2, eccDegs, stimStringForPhotocurrents);

    ax = subplot(2,2,2);
    generateCorrespondencePlot(ax, RsToRcExcitations, RsToRcPcurrent, ...
        LcenterRGCindices, McenterRGCindices, 'Rs/Rc', [3 12], 1, eccDegs, stimStringForPhotocurrents);

    ax = subplot(2,2,3);
    generateCorrespondencePlot(ax, KcExcitations/maxKc, KcPcurrent/maxKc, ...
        LcenterRGCindices, McenterRGCindices, 'Kc', [0 1], 0.2, eccDegs, stimStringForPhotocurrents);

    ax = subplot(2,2,4);
    generateCorrespondencePlot(ax, RcDegsExcitations*60 , RcDegsPcurrent*60, ...
        LcenterRGCindices, McenterRGCindices, 'Rc (arc min)', 3+[0 3], 0.5, eccDegs, stimStringForPhotocurrents);


    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = 'validation';
    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, sprintf('PcurrentVsConeExcitationsAt%2.1fdegsStim_%s.pdf', eccDegs, stimStringForPhotocurrents));

    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300);


end

function generateCorrespondencePlot(ax, excitations, pcurrents, ...
        LcenterRGCindices, McenterRGCindices, characteristic, xyLims, xyTickInterval, eccDegs, stimStringForPhotocurrents)
    scatter(ax, excitations(LcenterRGCindices), pcurrents(LcenterRGCindices), 121, ...
        'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0],  ...
        'MarkerFaceAlpha', 0.5, 'LineWidth', 1.0);
    hold(ax, 'on');
    scatter(ax, excitations(McenterRGCindices), pcurrents(McenterRGCindices), 121, ...
        'MarkerFaceColor', [0.5 1 0.5], 'MarkerEdgeColor', [0 0.8 0],  ...
        'MarkerFaceAlpha', 0.5, 'LineWidth', 1.0);
    plot(ax, xyLims, xyLims, 'k--', 'LineWidth', 1.5);
    set(ax, 'XLim', xyLims, 'YLim', xyLims);
    set(ax, 'XTick', xyLims(1): xyTickInterval: xyLims(2));
    set(ax, 'YTick', xyLims(1): xyTickInterval: xyLims(2));
    hold(ax, 'off');
    xtickangle(ax, 0);
    axis(ax, 'square');
    grid(ax, 'on');
    xlabel(ax, sprintf('photon absorptions - based'));
    ylabel(ax, sprintf('photo-current - based'));
    set(ax, 'FontSize', 20)
    title(ax, characteristic);
end
