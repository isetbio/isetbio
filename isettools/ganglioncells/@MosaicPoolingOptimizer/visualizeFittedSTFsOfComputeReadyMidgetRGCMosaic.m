function visualizeFittedSTFsOfComputeReadyMidgetRGCMosaic(...
    computeReadyMosaicFilenames, mRGCMosaicSTFresponsesFilenames, ...
    employTemporalEquivalentEccentricity)

    % Load the Croner&Kaplan '95 data for Parvocellular neurons
    [CK95RcDegs, CK95RsRcRatios, CK95KsKcRatios, CK95SCintSensRatios] = loadCronerKaplan95ParvoCellularData();
   
    % Initialize all data sets to empty
    for iDataSet = 1:3
        ISETBioRcDegs{iDataSet}.eccentricityDegs = [];
        ISETBioRsRcRatios{iDataSet}.eccentricityDegs = [];
        ISETBioKsKcRatios{iDataSet}.eccentricityDegs = [];
        ISETBioSCintSensRatios{iDataSet}.eccentricityDegs = [];

        ISETBioRcDegs{iDataSet}.data = [];
        ISETBioRsRcRatios{iDataSet}.data = [];
        ISETBioKsKcRatios{iDataSet}.data = [];
        ISETBioSCintSensRatios{iDataSet}.data = [];
    end
    

    % Concatenate data across mosaics
    for iMosaic = 1:numel(computeReadyMosaicFilenames)
        % Load the compute-ready MRGC mosaic
        load(computeReadyMosaicFilenames{iMosaic}, 'theComputeReadyMRGCmosaic');
    
        % Load the computed mRGC  STF responses 
        load(mRGCMosaicSTFresponsesFilenames{iMosaic}, ...
            'theMRGCMosaicOptimalSTFs', ...
            'visualRcDegsEstimates', ...
            'theMRGCMosaicSTFresponses', ...
            'theMRGCresponseTemporalSupportSeconds', ...
            'orientationsTested', 'spatialFrequenciesTested', ...
            'spatialPhasesDegs', 'coneContrasts');
    
     
        % Extract the corresponding params as measured for the compute-ready mRGCMosaic
        [cellEccDegs , cellCenterConeTypeWeights, ...
         cellRcDegs, cellRsRcRatios, cellKsKcRatios, cellSCintSensRatios] = ...
            extractComputeReadyMosaicData(theComputeReadyMRGCmosaic, theMRGCMosaicOptimalSTFs, employTemporalEquivalentEccentricity);
    

        LcenterRGCindices = find(squeeze(cellCenterConeTypeWeights(2,:)) < 0.85*squeeze(cellCenterConeTypeWeights(1,:)));
        McenterRGCindices = find(squeeze(cellCenterConeTypeWeights(1,:)) < 0.85*squeeze(cellCenterConeTypeWeights(2,:)));
        singleCellRGCindices = cat(2, LcenterRGCindices, McenterRGCindices);
        LMcenterRGCindices = setdiff(1:theComputeReadyMRGCmosaic.rgcsNum, singleCellRGCindices);
    
        % L-center RGCs data
        iDataSet = 1;
        % The RcDegs as a function of eccentricity for L-center RGCs
        ISETBioRcDegs{iDataSet}.eccentricityDegs = cat(2, ISETBioRcDegs{iDataSet}.eccentricityDegs, cellEccDegs(LcenterRGCindices)); %cellTemporalEquivalentEccDegs(LcenterRGCindices);
        ISETBioRcDegs{iDataSet}.data = cat(2, ISETBioRcDegs{iDataSet}.data, cellRcDegs(LcenterRGCindices));
    
        % The RsRcRatios as a function of eccentricity for L-center RGCs
        ISETBioRsRcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioRsRcRatios{iDataSet}.data = cat(2, ISETBioRsRcRatios{iDataSet}.data, cellRsRcRatios(LcenterRGCindices));
    
        % The KsRcRatios as a function of eccentricity for L-center RGCs
        ISETBioKsKcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioKsKcRatios{iDataSet}.data = cat(2, ISETBioKsKcRatios{iDataSet}.data, cellKsKcRatios(LcenterRGCindices));
    
        % The S/C int ratio as a function of eccentricity for L-center RGCs
        ISETBioSCintSensRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioSCintSensRatios{iDataSet}.data = cat(2,ISETBioSCintSensRatios{iDataSet}.data, cellSCintSensRatios(LcenterRGCindices));
    
        
        % M-center RGCs data
        iDataSet = 2;
        % The RcDegs as a function of eccentricity plot for M-center RGCs
        ISETBioRcDegs{iDataSet}.eccentricityDegs = cat(2, ISETBioRcDegs{iDataSet}.eccentricityDegs, cellEccDegs(McenterRGCindices)); %cellTemporalEquivalentEccDegs(McenterRGCindices);
        ISETBioRcDegs{iDataSet}.data = cat(2, ISETBioRcDegs{iDataSet}.data, cellRcDegs(McenterRGCindices));
    
        % The RsRcRatios as a function of eccentricity plot for M-center RGCs
        ISETBioRsRcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioRsRcRatios{iDataSet}.data = cat(2, ISETBioRsRcRatios{iDataSet}.data, cellRsRcRatios(McenterRGCindices));
    
        % The KsRcRatios as a function of eccentricity for M-center RGCs
        ISETBioKsKcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioKsKcRatios{iDataSet}.data = cat(2, ISETBioKsKcRatios{iDataSet}.data, cellKsKcRatios(McenterRGCindices));
    
        % The S/C int ratio as a function of eccentricity for M-center RGCs
        ISETBioSCintSensRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioSCintSensRatios{iDataSet}.data = cat(2,ISETBioSCintSensRatios{iDataSet}.data, cellSCintSensRatios(McenterRGCindices));
    
       
        % LM-center RGCs data
        iDataSet = 3;
        % The RcDegs as a function of eccentricity plot for L/M-center RGCs
        ISETBioRcDegs{iDataSet}.eccentricityDegs = cat(2, ISETBioRcDegs{iDataSet}.eccentricityDegs, cellEccDegs(LMcenterRGCindices)); % cellTemporalEquivalentEccDegs(LMcenterRGCindices);
        ISETBioRcDegs{iDataSet}.data = cat(2, ISETBioRcDegs{iDataSet}.data, cellRcDegs(LMcenterRGCindices));
    
        % The RsRcRatios as a function of eccentricity plot for LM-center RGCs
        ISETBioRsRcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioRsRcRatios{iDataSet}.data = cat(2, ISETBioRsRcRatios{iDataSet}.data, cellRsRcRatios(LMcenterRGCindices));
    
        % The KsRcRatios as a function of eccentricity for LM-center RGCs
        ISETBioKsKcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioKsKcRatios{iDataSet}.data = cat(2, ISETBioKsKcRatios{iDataSet}.data, cellKsKcRatios(LMcenterRGCindices));
    
        % The S/C int ratio as a function of eccentricity for LM-center RGCs
        ISETBioSCintSensRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioSCintSensRatios{iDataSet}.data = cat(2,ISETBioSCintSensRatios{iDataSet}.data, cellSCintSensRatios(LMcenterRGCindices));
    
    end % iMosaic


    % Transform Ks to log10(k)
    CK95KsKcRatiosLog = CK95KsKcRatios;
    ISETBioKsKcRatiosLog = ISETBioKsKcRatios;
    CK95KsKcRatiosLog.data = log10(CK95KsKcRatios.data);
    ISETBioKsKcRatiosLog{1}.data = log10(ISETBioKsKcRatios{1}.data);
    ISETBioKsKcRatiosLog{2}.data = log10(ISETBioKsKcRatios{2}.data);
    ISETBioKsKcRatiosLog{3}.data = log10(ISETBioKsKcRatios{3}.data);

    pdfFileNamePostfix{1} = 'L-cone center';
    pdfFileNamePostfix{2} = 'M-cone center';
    pdfFileNamePostfix{3} = 'mixed LM-cone center';

    eccRange = [0.3 30];
    eccTicks = [0.3 1 3 10 30];
    ISETBioSetsColors = [235 119 129; 82 180 140; 170 160 100]/255;
    scatterFaceAlpha = 0.05;
    histogramFaceAlpha = 0.3;


    for iDataSet = 1:3
        pngFileName{iDataSet} = sprintf('SummaryStats_%s.png',pdfFileNamePostfix{iDataSet});
        
        hFig = figure(iDataSet); clf;
        ff = MSreadyPlot.figureFormat('2x4');
        theAxes = MSreadyPlot.generateAxes(hFig,ff);
        set(hFig, 'Color', [1 1 1]);

        if (iDataSet == 1)
            theVisualizedDataSets = [1 2];
        elseif (iDataSet == 2)
            theVisualizedDataSets = [2,1];
        else
            theVisualizedDataSets = iDataSet;
        end

        % The RcDegs as a function of eccentricity plot
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParams(theAxes{1,1}, ...
                  CK95RcDegs, ISETBioRcDegs, theVisualizedDataSets, ...
                  ISETBioSetsColors, scatterFaceAlpha, ...
                  eccRange, eccTicks, ...
                  [0 0.1], (0:0.01:0.1), ...
                  'linear', 'Rc (degs)', pdfFileNamePostfix{iDataSet}, ...
                  employTemporalEquivalentEccentricity, ...
                  ff);
    

        % The Rs/Rc ratio as a function of eccentricity plot
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParams(theAxes{1,2}, ...
                  CK95RsRcRatios, ISETBioRsRcRatios, theVisualizedDataSets, ...
                  ISETBioSetsColors, scatterFaceAlpha, ...
                  eccRange, eccTicks, ...
                  [0 20], (0:2:20), ...
                  'linear', 'Rs/Rc', '', ...
                  employTemporalEquivalentEccentricity, ...
                  ff);

        % The Rs/Rc ratio histograms
        histogramEdges = 0:0.5:20;
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsHistograms(theAxes{2,2}, ...
              CK95RsRcRatios, ISETBioRsRcRatios, theVisualizedDataSets, histogramEdges, ...
              ISETBioSetsColors, histogramFaceAlpha, ...
              [0 20], (0:2:20), ...
              'Rs/Rc', '', ...
              ff);


        % The Ks/Kc int ratio as a function of eccentricity plot
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParams(theAxes{1,3}, ...
              CK95KsKcRatios, ISETBioKsKcRatios, theVisualizedDataSets,  ...
              ISETBioSetsColors, scatterFaceAlpha, ...
              eccRange, eccTicks, ...
              [1e-4 1e-1], [1e-4 1e-3 1e-2 1e-1], ...
              'log', 'Ks/Kc', '', ...
              employTemporalEquivalentEccentricity, ...
              ff);


        % The Ks/Kc histograms
        histogramEdges = -4:0.1:0;
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsHistograms(theAxes{2,3}, ...
              CK95KsKcRatiosLog, ISETBioKsKcRatiosLog, theVisualizedDataSets, histogramEdges, ...
              ISETBioSetsColors, histogramFaceAlpha, ...
              [-4 -1], (-4:0.5:0), ...
              'log10(Ks/Kc)', '', ...
              ff);


        % The S/C int sens ratio as a function of eccentricity
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParams(theAxes{1,4}, ...
              CK95SCintSensRatios, ISETBioSCintSensRatios, theVisualizedDataSets, ...
              ISETBioSetsColors, scatterFaceAlpha, ...
              eccRange, eccTicks, ...
              [0.1 0.9], (0:0.1:1), ...
              'linear', 'S/C int. sens.', '', ...
              employTemporalEquivalentEccentricity, ...
              ff);

        % The S/c int ratio histograms
        histogramEdges = 0:0.02:1.0;
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsHistograms(theAxes{2,4}, ...
              CK95SCintSensRatios, ISETBioSCintSensRatios, theVisualizedDataSets, histogramEdges, ...
              ISETBioSetsColors, histogramFaceAlpha, ...
              [0.1 0.9], (0:0.1:1), ...
              'S/C int. sens.', '', ...
              ff);

        NicePlot.exportFigToPNG(pngFileName{iDataSet}, hFig, 300);
    end % iDataSet
end

function [cellEccDegs, cellCenterConeTypeWeights, ...
          cellRcDegs, cellRsRcRatios, cellKsKcRatios, ...
          cellSCintSensRatios] = extractComputeReadyMosaicData(...
                                    theComputeReadyMRGCmosaic, ...
                                    theMRGCMosaicOptimalSTFs, ...
                                    employTemporalEquivalentEccentricity)

    rgcIndex = 1;
    theSTFdata = theMRGCMosaicOptimalSTFs{rgcIndex};
    idxRcDegs = find(strcmp(theSTFdata.DoGfitParams.names, 'RcDegs'));
    idxRsRcRatio = find(strcmp(theSTFdata.DoGfitParams.names, 'RsToRc'));
    idxKsKcRatio = find(strcmp(theSTFdata.DoGfitParams.names, 'kS/kC'));
    
    cellEccDegs = zeros(1,theComputeReadyMRGCmosaic.rgcsNum);
    cellCenterConeTypeWeights = zeros(3,theComputeReadyMRGCmosaic.rgcsNum);

    cellRcDegs = cellEccDegs;
    cellRsRcRatios = cellEccDegs;
    cellKsKcRatios = cellEccDegs;

    for iRGC = 1:theComputeReadyMRGCmosaic.rgcsNum

        theSTFdata = theMRGCMosaicOptimalSTFs{iRGC};

%         theSTFdata.measured
%         theSTFdata.DoGfit
%         theSTFdata.DoGfitParams

        % Compute equivalent temporal eccentricity
        theRGCRFpositionDegs = theComputeReadyMRGCmosaic.rgcRFpositionsDegs(iRGC,:);
        if (employTemporalEquivalentEccentricity)
            temporalEquivEccDegs = theComputeReadyMRGCmosaic.temporalEquivalentEccentricityForEccentricity(theRGCRFpositionDegs);
            cellTemporalEquivalentEccDegs = sqrt(sum(temporalEquivEccDegs.^2,2));
            cellEccDegs(iRGC) = cellTemporalEquivalentEccDegs;
        else
            cellEccDegs(iRGC) = sqrt(sum(theRGCRFpositionDegs.^2,2));
        end

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
