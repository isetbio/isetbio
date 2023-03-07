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


    [CK.eccentricityDegs, cK.RcDegs] = RGCmodels.CronerKaplan.digitizedData.parvoCenterRadiusAgainstEccentricity();

    figure(1); clf;
    subplot(2,2,1);
    plot(CK.eccentricityDegs, cK.RcDegs, 's'); hold on
    plot(cellEccDegs, cellRcDegs, 'ro');
    ylabel('Rc (degs)');
    xlabel('ecc (degs)');
    axis('square');

    subplot(2,2,2);
    plot(cellEccDegs, cellRsRcRatios, 'ro');
    ylabel('Rs/Rc');
    xlabel('ecc (degs)');
    axis('square');

    subplot(2,2,3);
    plot(cellEccDegs, cellKsKcRatios, 'ro');
    ylabel('Ks/Kc');
    xlabel('ecc (degs)');
    axis('square');

    subplot(2,2,4);
    plot(cellEccDegs, cellSCintSensRatios, 'ro');
    ylabel('S/C int sens');
    xlabel('ecc (degs)');
    axis('square');

end

