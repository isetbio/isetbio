function performVisualizeConePSFsWithinRGCMosaicOp(mosaicParams, tickSeparationArcMin)

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    [~, ~, pdfsDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('pdfsDirectory', ...
            'mosaicParams', mosaicParams);

    % Load the generated center-only connected mRGCmosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');


    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    fprintf('\n---> Select the optics to use for computing the cone PSFs\n');
    [opticsParamsForComputingConePSFs, opticsToEmploy] = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);


    % Generate and set the optics
    theMidgetRGCMosaic.setTheOptics(opticsParamsForComputingConePSFs);

    targetPositionDegs = [];
    while (numel(targetPositionDegs) ~= 2)
        targetPositionDegs = input('Enter position within the RGC mosaic (e.g. [0.7 0.7]) :');
    end

    targetSizeDegs = 1;
    [theLconePSF, theMconePSF, theSconePSF, spatialSupportDegs] = ...
        MosaicPoolingOptimizer.computeConePSFs(theMidgetRGCMosaic, opticsToEmploy, targetSizeDegs, targetPositionDegs);

    % ============== Export to PLOS directory ==========================
    rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';


    
    pdfFileName = sprintf('ConePSF_mosaicEccDegs_%2.1f_%2.1f_positionWithinMosaic_%2.2f_%2.2f.pdf', ...
        mosaicParams.eccDegs(1), mosaicParams.eccDegs(1), targetPositionDegs(1), targetPositionDegs(2));
    psfRangeDegs = 0.5*(tickSeparationArcMin*4)/60;

    
    hFig = figure(555); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    MSreadyPlot.render2DPSF(theAxes{1,1}, spatialSupportDegs, spatialSupportDegs, theLconePSF, psfRangeDegs, 'L-cone PSF', ff, ...
        'tickSeparationArcMin', tickSeparationArcMin);

    pdfFileNameForPLOS = fullfile(rawFiguresRoot, strrep(pdfFileName, 'ConePSF', 'LconePSF'));
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);

    hFig = figure(556); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    MSreadyPlot.render2DPSF(theAxes{1,1}, spatialSupportDegs, spatialSupportDegs, theMconePSF, psfRangeDegs, 'M-cone PSF', ff, ...
        'tickSeparationArcMin', tickSeparationArcMin);

    pdfFileNameForPLOS = fullfile(rawFiguresRoot, strrep(pdfFileName, 'ConePSF', 'MconePSF'));
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);

    hFig = figure(557); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    MSreadyPlot.render2DPSF(theAxes{1,1}, spatialSupportDegs, spatialSupportDegs, theSconePSF, psfRangeDegs, 'S-cone PSF', ff, ...
        'tickSeparationArcMin', tickSeparationArcMin);

    pdfFileNameForPLOS = fullfile(rawFiguresRoot, strrep(pdfFileName, 'ConePSF', 'SconePSF'));
    NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig, 300);

end