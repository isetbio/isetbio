function coneWeightedPSFonTopOfConeMosaic(figNo, theInputConeMosaic, theConeWeightedPSF, thePDFfileName, theTitle)

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    ax = theAxes{1,1};
    visualizedWavelengthIndex = 1;
    maxPSF = 1;
    infoString = '';
    xo = theInputConeMosaic.eccentricityDegs(1);
    yo = theInputConeMosaic.eccentricityDegs(2);
    RGCMosaicConstructor.helper.optics.visualizePSFAtWavelengthOnTopOfConeMosaic(...
        hFig, ax, theInputConeMosaic, theConeWeightedPSF, ...
        visualizedWavelengthIndex, theConeWeightedPSF, ...
        theTitle, ...
        'XLimsArcMin', [-12.5 12.5], ...
        'YLimsArcMin', [-12.5 12.5], ...
        'XTicksArcMin', xo*60 + (-12:6:12), ...
        'YTicksArcMin', yo*60 + (-12:6:12), ...
        'withInfoString', infoString);
    
    ff.box = 'on';
    ff.tickDir = 'in';
    xlabel(ax, 'eccentricity, x (degs)');
    ylabel(ax, 'eccentricity, y (degs)');
    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    %PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

    % OLD WAY
    %theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();

    %if (isdir(theRawFiguresDir))
    %    pdfExportSubDir = 'validation';
    
    %    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, thePDFfileName);
    %    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    %else
    %    fprintf('Raw figures directory not found: %s.\nWill not export PDF\n', theRawFiguresDir);
    %end


    p = getpref('isetbio');
    pdfExportSubDir = fullfile(p.rgcResources.figurePDFsDir);
    theVisualizationPDFfilename = fullfile('modelComponentVisualizations', thePDFfileName);
    
    % Generate the path if we need to
    RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
        pdfExportSubDir, theVisualizationPDFfilename, ...
        'generateMissingSubDirs', true);

    thePDFfileName = fullfile(pdfExportSubDir, theVisualizationPDFfilename);
    NicePlot.exportFigToPDF(thePDFfileName, hFig, 300);

end