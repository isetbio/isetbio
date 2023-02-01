function visualizeMosaicSTFfits(mosaicCenterParams, rfModelParams, opticsParams)

    % Generate the frozen mosaic filename
    frozenMosaicFileName = midgetRGCMosaicInspector.frozenMosaicFileName(...
        mosaicCenterParams, rfModelParams.H1cellIndex, opticsParams);
    
    % Load the frozen midget RGC mosaic
    load(frozenMosaicFileName, 'theMidgetRGCmosaic');


    % Ask the user to specify the optics position for which responses were saved
    opticsPositionDegs = [];
    while (numel(opticsPositionDegs) ~= 2)
        opticsPositionDegs = input('\nEnter the optics position that was used to compute the responses ([x y]): ');
    end

    % Generate the responses filename
    responsesFileName = midgetRGCMosaicInspector.responsesFileName(...
        frozenMosaicFileName, opticsPositionDegs);
       
    load(responsesFileName, 'theMeridianAngles', 'theMeridianFits');


    hFigRcDegs = [];
    hFigRsRcRatios = [];
    hFigSCintSensRatios = [];
    for iMeridianAngle = 1:numel(theMeridianAngles)

        [hFigRcDegs, hFigRsRcRatios, hFigSCintSensRatios] = ...
            midgetRGCMosaicInspector.renderSTFfitPlots(...
                hFigRcDegs, hFigRsRcRatios, hFigSCintSensRatios, ...
                theMidgetRGCmosaic.rgcRFpositionsDegs, ...
                theMeridianAngles, iMeridianAngle, ...
                theMeridianFits{iMeridianAngle});

    end % iMeridianAngle


    % Export to PDF
    hFig =  figure(hFigRcDegs);
    pdfFileName = strrep(responsesFileName, '.mat', '_FittedSTF_RcDegs.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);


    hFig =  figure(hFigRsRcRatios);
    pdfFileName = strrep(responsesFileName, '.mat', '_FittedSTF_RcRsRatios.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);


    hFig =  figure(hFigSCintSensRatios);
    pdfFileName = strrep(responsesFileName, '.mat', '_FittedSTF_SCintSensRatios.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

