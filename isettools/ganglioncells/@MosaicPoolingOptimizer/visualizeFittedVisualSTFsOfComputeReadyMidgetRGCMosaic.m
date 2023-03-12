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

    % Load the Croner&Kaplan '95 data for Parvocellular neurons
    [CK95RcDegs, CK95RsRcRatios, CK95KsKcRatios, CK95SCintSensRatios] = loadCronerKaplan95ParvoCellularData();
    
    % Extract the corresponding params as measured for the compute-ready mRGCMosaic
    [cellTemporalEquivalentEccDegs, cellEccDegs , cellCenterConeTypeWeights, ...
     cellRcDegs, cellRsRcRatios, cellKsKcRatios, cellSCintSensRatios] = ...
        extractComputeReadyMosaicData(theComputeReadyMRGCmosaic, theMRGCMosaicOptimalSTFs);

    LcenterRGCindices = find(squeeze(cellCenterConeTypeWeights(2,:)) < 0.85*squeeze(cellCenterConeTypeWeights(1,:)));
    McenterRGCindices = find(squeeze(cellCenterConeTypeWeights(1,:)) < 0.85*squeeze(cellCenterConeTypeWeights(2,:)));
    singleCellRGCindices = cat(2, LcenterRGCindices, McenterRGCindices);
    LMcenterRGCindices = setdiff(1:theComputeReadyMRGCmosaic.rgcsNum, singleCellRGCindices);

    pdfFileNamePostfix{1} = 'L-cone center';
    pdfFileNamePostfix{2} = 'M-cone center';
    pdfFileNamePostfix{3} = 'mixed LM-cone center';

    for visualizedISETBioDataSets = 1:3
    pngFileName = strrep(computeReadyMosaicFilename, '.mat', sprintf('_summaryStats_%s.png',pdfFileNamePostfix{visualizedISETBioDataSets}))
    

    hFig = figure(visualizedISETBioDataSets); clf;
    ff = MSreadyPlot.figureFormat('2x4');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);

    eccRange = [0.3 30];
    eccTicks = [0.3 1 3 10 30];
    ISETBioSetsColors = [1 0 0; 0 1 0; 1 0.7 0];
    scatterFaceAlpha = 0.05;
    histogramFaceAlpha = 0.3;

    % The RcDegs as a function of eccentricity plot
    ISETBioData{1}.eccentricityDegs = cellEccDegs(LcenterRGCindices); %cellTemporalEquivalentEccDegs(LcenterRGCindices);
    ISETBioData{1}.data = cellRcDegs(LcenterRGCindices);
    ISETBioData{2}.eccentricityDegs = cellEccDegs(McenterRGCindices); %cellTemporalEquivalentEccDegs(McenterRGCindices);
    ISETBioData{2}.data = cellRcDegs(McenterRGCindices);
    ISETBioData{3}.eccentricityDegs = cellEccDegs(LMcenterRGCindices); % cellTemporalEquivalentEccDegs(LMcenterRGCindices);
    ISETBioData{3}.data = cellRcDegs(LMcenterRGCindices);

    MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParams(theAxes{1,1}, ...
          CK95RcDegs, ISETBioData, visualizedISETBioDataSets, ...
          ISETBioSetsColors, scatterFaceAlpha, ...
          eccRange, eccTicks, ...
          [0 0.1], (0:0.01:0.1), ...
          'linear', 'Rc (degs)', '', ...
          ff);
    

    % The Rs/Rc ratio as a function of eccentricity plot
    ISETBioData{1}.data = cellRsRcRatios(LcenterRGCindices);
    ISETBioData{2}.data = cellRsRcRatios(McenterRGCindices);
    ISETBioData{3}.data = cellRsRcRatios(LMcenterRGCindices);

    MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParams(theAxes{1,2}, ...
          CK95RsRcRatios, ISETBioData, visualizedISETBioDataSets, ...
          ISETBioSetsColors, scatterFaceAlpha, ...
          eccRange, eccTicks, ...
          [0 20], (0:2:20), ...
          'linear', 'Rs/Rc', '', ...
          ff);

    % The Rs/Rc histograms
    histogramEdges = 0:0.5:20;
    MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsHistograms(theAxes{2,2}, ...
          CK95RsRcRatios, ISETBioData, visualizedISETBioDataSets, histogramEdges, ...
          ISETBioSetsColors, histogramFaceAlpha, ...
          [0 20], (0:2:20), ...
          'Rs/Rc', '', ...
          ff);


    % The Ks/Kc int ratio as a function of eccentricity plot
    ISETBioData{1}.data = cellKsKcRatios(LcenterRGCindices);
    ISETBioData{2}.data = cellKsKcRatios(McenterRGCindices);
    ISETBioData{3}.data = cellKsKcRatios(LMcenterRGCindices);
    MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParams(theAxes{1,3}, ...
          CK95KsKcRatios, ISETBioData, visualizedISETBioDataSets,  ...
          ISETBioSetsColors, scatterFaceAlpha, ...
          eccRange, eccTicks, ...
          [1e-4 1e-1], [1e-4 1e-3 1e-2 1e-1], ...
          'log', 'Ks/Kc', '', ...
          ff);


    % The Ks/Kc histograms
    % Transform Ks to log10(k)
    CK95KsKcRatiosLog = CK95KsKcRatios;
    ISETBioDataLog = ISETBioData;
    CK95KsKcRatiosLog.data = log10(CK95KsKcRatios.data);
    ISETBioDataLog{1}.data = log10(ISETBioData{1}.data);
    ISETBioDataLog{2}.data = log10(ISETBioData{2}.data);
    ISETBioDataLog{3}.data = log10(ISETBioData{3}.data);

    histogramEdges = -4:0.1:0;
    MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsHistograms(theAxes{2,3}, ...
          CK95KsKcRatiosLog, ISETBioDataLog, visualizedISETBioDataSets, histogramEdges, ...
          ISETBioSetsColors, histogramFaceAlpha, ...
          [-4 -1], (-4:0.5:0), ...
          'log10(Ks/Kc)', '', ...
          ff);


    % The S/C int ratio as a function of eccentricity plot
    ISETBioData{1}.data = cellSCintSensRatios(LcenterRGCindices);
    ISETBioData{2}.data = cellSCintSensRatios(McenterRGCindices);
    ISETBioData{3}.data = cellSCintSensRatios(LMcenterRGCindices);
    MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParams(theAxes{1,4}, ...
          CK95SCintSensRatios, ISETBioData,  visualizedISETBioDataSets, ...
          ISETBioSetsColors, scatterFaceAlpha, ...
          eccRange, eccTicks, ...
          [0.1 0.9], (0:0.1:1), ...
          'linear', 'S/C int. sens.', '', ...
          ff);

    % The S/c int ratio histograms
    histogramEdges = 0:0.02:1.0;
    MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsHistograms(theAxes{2,4}, ...
          CK95SCintSensRatios, ISETBioData, visualizedISETBioDataSets, histogramEdges, ...
          ISETBioSetsColors, histogramFaceAlpha, ...
          [0.1 0.9], (0:0.1:1), ...
          'S/C int. sens.', '', ...
          ff);

    NicePlot.exportFigToPNG(pngFileName, hFig, 300);
    end
end

function [cellTemporalEquivalentEccDegs, cellEccDegs, cellCenterConeTypeWeights, ...
          cellRcDegs, cellRsRcRatios, cellKsKcRatios, cellSCintSensRatios] = ...
                  extractComputeReadyMosaicData(theComputeReadyMRGCmosaic, theMRGCMosaicOptimalSTFs)

    rgcIndex = 1;
    theSTFdata = theMRGCMosaicOptimalSTFs{rgcIndex};
    idxRcDegs = find(strcmp(theSTFdata.DoGfitParams.names, 'RcDegs'));
    idxRsRcRatio = find(strcmp(theSTFdata.DoGfitParams.names, 'RsToRc'));
    idxKsKcRatio = find(strcmp(theSTFdata.DoGfitParams.names, 'kS/kC'));
    
    cellEccDegs = zeros(1,theComputeReadyMRGCmosaic.rgcsNum);
    cellTemporalEquivalentEccDegs = zeros(1,theComputeReadyMRGCmosaic.rgcsNum);
    cellCenterConeTypeWeights = zeros(3,theComputeReadyMRGCmosaic.rgcsNum);

    cellRcDegs = cellTemporalEquivalentEccDegs;
    cellRsRcRatios = cellTemporalEquivalentEccDegs;
    cellKsKcRatios = cellTemporalEquivalentEccDegs;

    for iRGC = 1:theComputeReadyMRGCmosaic.rgcsNum

        theSTFdata = theMRGCMosaicOptimalSTFs{iRGC};

%         theSTFdata.measured
%         theSTFdata.DoGfit
%         theSTFdata.DoGfitParams

        % Compute equivalent temporal eccentricity
        theRGCRFpositionDegs = theComputeReadyMRGCmosaic.rgcRFpositionsDegs(iRGC,:);
        temporalEquivEccDegs = theComputeReadyMRGCmosaic.temporalEquivalentEccentricityForEccentricity(theRGCRFpositionDegs);
        cellTemporalEquivalentEccDegs(iRGC) = sqrt(sum(temporalEquivEccDegs.^2,2));
        cellEccDegs(iRGC) = sqrt(sum(theRGCRFpositionDegs.^2,2));

        % Cell center cone types
        [centerConeTypeWeightsTmp, ~, ~, cellCenterConeTypes] = theComputeReadyMRGCmosaic.centerConeTypeWeights(iRGC);
        
        % CenterConeTypeWeights are sorted according to decreasing weight
        % Re-sort them so that L is first and M second
        centerConeTypeWeights(1) = centerConeTypeWeightsTmp(find(cellCenterConeTypes == cMosaic.LCONE_ID));
        centerConeTypeWeights(2) = centerConeTypeWeightsTmp(find(cellCenterConeTypes == cMosaic.MCONE_ID));
        centerConeTypeWeights(3) = centerConeTypeWeightsTmp(find(cellCenterConeTypes == cMosaic.SCONE_ID));

        cellCenterConeTypeWeights(:,iRGC) = reshape(centerConeTypeWeights, [3 1]);

        cellRcDegs(iRGC) = theSTFdata.DoGfitParams.finalValues(idxRcDegs);
        cellRsRcRatios(iRGC) = theSTFdata.DoGfitParams.finalValues(idxRsRcRatio);
        cellKsKcRatios(iRGC) = theSTFdata.DoGfitParams.finalValues(idxKsKcRatio);
        
    end

    cellSCintSensRatios = cellKsKcRatios .* cellRsRcRatios.^2;
end


function [CK95RcDegs, CK95RsRcRatios, CK95KsKcRatios, CK95SCintSensRatios] = loadCronerKaplan95ParvoCellularData()
    
    [CK95RcDegs.eccentricityDegs, CK95RcDegs.data] = ...
        RGCmodels.CronerKaplan.digitizedData.parvoCenterRadiusAgainstEccentricity();

    [CK95RsRcRatios.eccentricityDegs, CK95RsRcRatios.data] = ...
        RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();
    CK95RsRcRatios.data = 1./CK95RsRcRatios.data;

    [CK95KsKcRatios.eccentricityDegs, CK95KsKcRatios.data] = ...
        RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterPeakSensisitivityRatioAgainstEccentricity();

    [CK95SCintSensRatios.eccentricityDegs, CK95SCintSensRatios.data] = ...
       RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();
end
