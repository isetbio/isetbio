function visualizeFittedSTFsOfComputeReadyMidgetRGCMosaic(...
    computeReadyMosaicFilenames, mRGCMosaicSTFresponsesFilenames, ...
    pdfDirectory, ...
    showZscoresInsteadOfData, onlyShowCronerKaplan95Data, ...
    employTemporalEquivalentEccentricity)

    % Load the Croner&Kaplan '95 data for Parvocellular neurons
    [CK95RcDegs, CK95KcVsRcDegs, CK95RsRcRatios, CK95KsKcRatios, CK95SCintSensRatios] = ...
        loadCronerKaplan95ParvoCellularData();
   
    % Initialize all data sets to empty
    ISETBioRcDegs = cell(1,4);
    ISETBioKc = cell(1,4);
    ISETBioRsRcRatios = cell(1,4);
    ISETBioKsKcRatios = cell(1,4);
    ISETBioSCintSensRatios = cell(1,4);

    for iDataSet = 1:4
        % RcDegs
        ISETBioRcDegs{iDataSet}.eccentricityDegs = [];
        ISETBioRcDegs{iDataSet}.data = [];

        % Kc
        ISETBioKc{iDataSet}.eccentricityDegs = [];
        ISETBioKc{iDataSet}.data = [];

        % RsRc ratios
        ISETBioRsRcRatios{iDataSet}.eccentricityDegs = [];
        ISETBioRsRcRatios{iDataSet}.data = [];

        % KsKc ratios
        ISETBioKsKcRatios{iDataSet}.eccentricityDegs = [];
        ISETBioKsKcRatios{iDataSet}.data = [];

        % Integrated sensitivity ratios
        ISETBioSCintSensRatios{iDataSet}.eccentricityDegs = [];
        ISETBioSCintSensRatios{iDataSet}.data = [];
    end
    

    % Concatenate data across mosaics
    mosaicOutlines = cell(1, numel(computeReadyMosaicFilenames));
    for iMosaic = 1:numel(computeReadyMosaicFilenames)
        % Load the compute-ready MRGC mosaic
        load(computeReadyMosaicFilenames{iMosaic}, 'theComputeReadyMRGCmosaic');
    
        outlineStruct.eccDegs = theComputeReadyMRGCmosaic.eccentricityDegs;
        outlineStruct.sizeDegs = theComputeReadyMRGCmosaic.sizeDegs;
        outlineStruct.sizeDegs(2) = outlineStruct.sizeDegs(2)+(iMosaic-1)*0.2;
        outlineStruct.temporalEquivalentEccDegs = theComputeReadyMRGCmosaic.temporalEquivalentEccentricityDegs;
        mosaicOutlines{iMosaic}.outlineStruct = outlineStruct;
        

        % Load the computed mRGC  STF responses 
        load(mRGCMosaicSTFresponsesFilenames{iMosaic}, ...
            'theMRGCMosaicOptimalSTFs', ...
            'visualRcDegsEstimates', ...
            'theMRGCMosaicSTFresponses', ...
            'theMRGCresponseTemporalSupportSeconds', ...
            'orientationsTested', 'spatialFrequenciesTested', ...
            'spatialPhasesDegs', 'coneContrasts');
    
     
        % Extract the corresponding params as measured for the compute-ready mRGCMosaic
        [cellEccDegs, cellCenterConeTypeWeights, ...
         cellRcDegs, cellKc, cellRsRcRatios, cellKsKcRatios, cellSCintSensRatios] = ...
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
        % All-center RGCs RcDegs data
        ISETBioRcDegs{4}.eccentricityDegs = cat(2, ISETBioRcDegs{4}.eccentricityDegs, ISETBioRcDegs{iDataSet}.eccentricityDegs); %cellTemporalEquivalentEccDegs(LcenterRGCindices);
        ISETBioRcDegs{4}.data = cat(2, ISETBioRcDegs{4}.data, ISETBioRcDegs{iDataSet}.data);
    
        % The Kc as a function of eccentricity for L-center RGCs
        ISETBioKc{iDataSet}.eccentricityDegs = cat(2, ISETBioKc{iDataSet}.eccentricityDegs, cellKc(LcenterRGCindices)); %cellTemporalEquivalentEccDegs(LcenterRGCindices);
        ISETBioKc{iDataSet}.data = cat(2, ISETBioKc{iDataSet}.data, cellKc(LcenterRGCindices));
        % All-center RGCs Kc data
        ISETBioKc{4}.eccentricityDegs = cat(2, ISETBioKc{4}.eccentricityDegs, ISETBioKc{iDataSet}.eccentricityDegs); %cellTemporalEquivalentEccDegs(LcenterRGCindices);
        ISETBioKc{4}.data = cat(2, ISETBioKc{4}.data, ISETBioKc{iDataSet}.data);
    

        % The RsRcRatios as a function of eccentricity for L-center RGCs
        ISETBioRsRcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioRsRcRatios{iDataSet}.data = cat(2, ISETBioRsRcRatios{iDataSet}.data, cellRsRcRatios(LcenterRGCindices));
        % All-center RGCs RsRcRatios data
        ISETBioRsRcRatios{4}.eccentricityDegs = cat(2, ISETBioRsRcRatios{4}.eccentricityDegs, ISETBioRsRcRatios{iDataSet}.eccentricityDegs);
        ISETBioRsRcRatios{4}.data = cat(2, ISETBioRsRcRatios{4}.data, ISETBioRsRcRatios{iDataSet}.data);
    

        % The KsRcRatios as a function of eccentricity for L-center RGCs
        ISETBioKsKcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioKsKcRatios{iDataSet}.data = cat(2, ISETBioKsKcRatios{iDataSet}.data, cellKsKcRatios(LcenterRGCindices));
        % All-center RGCs KsRcRatios data
        ISETBioKsKcRatios{4}.eccentricityDegs = cat(2, ISETBioKsKcRatios{4}.eccentricityDegs, ISETBioKsKcRatios{iDataSet}.eccentricityDegs);
        ISETBioKsKcRatios{4}.data = cat(2, ISETBioKsKcRatios{4}.data, ISETBioKsKcRatios{iDataSet}.data);
    

        % The S/C int ratio as a function of eccentricity for L-center RGCs
        ISETBioSCintSensRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioSCintSensRatios{iDataSet}.data = cat(2,ISETBioSCintSensRatios{iDataSet}.data, cellSCintSensRatios(LcenterRGCindices));
        % All-center RGCs SCintSensRatios data
        ISETBioSCintSensRatios{4}.eccentricityDegs = cat(2, ISETBioSCintSensRatios{4}.eccentricityDegs, ISETBioSCintSensRatios{iDataSet}.eccentricityDegs);
        ISETBioSCintSensRatios{4}.data = cat(2,ISETBioSCintSensRatios{4}.data, ISETBioSCintSensRatios{iDataSet}.data);
        

        % M-center RGCs data
        iDataSet = 2;

        % The RcDegs as a function of eccentricity plot for M-center RGCs
        ISETBioRcDegs{iDataSet}.eccentricityDegs = cat(2, ISETBioRcDegs{iDataSet}.eccentricityDegs, cellEccDegs(McenterRGCindices)); %cellTemporalEquivalentEccDegs(McenterRGCindices);
        ISETBioRcDegs{iDataSet}.data = cat(2, ISETBioRcDegs{iDataSet}.data, cellRcDegs(McenterRGCindices));
        % All-center RGCs RcDegs data
        ISETBioRcDegs{4}.eccentricityDegs = cat(2, ISETBioRcDegs{4}.eccentricityDegs, ISETBioRcDegs{iDataSet}.eccentricityDegs); %cellTemporalEquivalentEccDegs(LcenterRGCindices);
        ISETBioRcDegs{4}.data = cat(2, ISETBioRcDegs{4}.data, ISETBioRcDegs{iDataSet}.data);
    

        % The Kc as a function of eccentricity for M-center RGCs
        ISETBioKc{iDataSet}.eccentricityDegs = cat(2, ISETBioKc{iDataSet}.eccentricityDegs, cellKc(McenterRGCindices)); %cellTemporalEquivalentEccDegs(LcenterRGCindices);
        ISETBioKc{iDataSet}.data = cat(2, ISETBioKc{iDataSet}.data, cellKc(McenterRGCindices));
        % All-center RGCs Kc data
        ISETBioKc{4}.eccentricityDegs = cat(2, ISETBioKc{4}.eccentricityDegs, ISETBioKc{iDataSet}.eccentricityDegs); %cellTemporalEquivalentEccDegs(LcenterRGCindices);
        ISETBioKc{4}.data = cat(2, ISETBioKc{4}.data, ISETBioKc{iDataSet}.data);
    

        % The RsRcRatios as a function of eccentricity plot for M-center RGCs
        ISETBioRsRcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioRsRcRatios{iDataSet}.data = cat(2, ISETBioRsRcRatios{iDataSet}.data, cellRsRcRatios(McenterRGCindices));
        % All-center RGCs RsRcRatios data
        ISETBioRsRcRatios{4}.eccentricityDegs = cat(2, ISETBioRsRcRatios{4}.eccentricityDegs, ISETBioRsRcRatios{iDataSet}.eccentricityDegs);
        ISETBioRsRcRatios{4}.data = cat(2, ISETBioRsRcRatios{4}.data, ISETBioRsRcRatios{iDataSet}.data);
   

        % The KsRcRatios as a function of eccentricity for M-center RGCs
        ISETBioKsKcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioKsKcRatios{iDataSet}.data = cat(2, ISETBioKsKcRatios{iDataSet}.data, cellKsKcRatios(McenterRGCindices));
        % All-center RGCs KsRcRatios data
        ISETBioKsKcRatios{4}.eccentricityDegs = cat(2, ISETBioKsKcRatios{4}.eccentricityDegs, ISETBioKsKcRatios{iDataSet}.eccentricityDegs);
        ISETBioKsKcRatios{4}.data = cat(2, ISETBioKsKcRatios{4}.data, ISETBioKsKcRatios{iDataSet}.data);
    
        % The S/C int ratio as a function of eccentricity for M-center RGCs
        ISETBioSCintSensRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioSCintSensRatios{iDataSet}.data = cat(2,ISETBioSCintSensRatios{iDataSet}.data, cellSCintSensRatios(McenterRGCindices));
        % All-center RGCs SCintSensRatios data
        ISETBioSCintSensRatios{4}.eccentricityDegs = cat(2, ISETBioSCintSensRatios{4}.eccentricityDegs, ISETBioSCintSensRatios{iDataSet}.eccentricityDegs);
        ISETBioSCintSensRatios{4}.data = cat(2,ISETBioSCintSensRatios{4}.data, ISETBioSCintSensRatios{iDataSet}.data);
        
       
        % LM-center RGCs data
        iDataSet = 3;

        % The RcDegs as a function of eccentricity plot for L/M-center RGCs
        ISETBioRcDegs{iDataSet}.eccentricityDegs = cat(2, ISETBioRcDegs{iDataSet}.eccentricityDegs, cellEccDegs(LMcenterRGCindices)); % cellTemporalEquivalentEccDegs(LMcenterRGCindices);
        ISETBioRcDegs{iDataSet}.data = cat(2, ISETBioRcDegs{iDataSet}.data, cellRcDegs(LMcenterRGCindices));
        % All-center RGCs RcDegs data
        ISETBioRcDegs{4}.eccentricityDegs = cat(2, ISETBioRcDegs{4}.eccentricityDegs, ISETBioRcDegs{iDataSet}.eccentricityDegs); %cellTemporalEquivalentEccDegs(LcenterRGCindices);
        ISETBioRcDegs{4}.data = cat(2, ISETBioRcDegs{4}.data, ISETBioRcDegs{iDataSet}.data);
    
        % The Kc as a function of eccentricity for L/M-center RGCs
        ISETBioKc{iDataSet}.eccentricityDegs = cat(2, ISETBioKc{iDataSet}.eccentricityDegs, cellKc(LMcenterRGCindices)); %cellTemporalEquivalentEccDegs(LcenterRGCindices);
        ISETBioKc{iDataSet}.data = cat(2, ISETBioKc{iDataSet}.data, cellKc(LMcenterRGCindices));
        % All-center RGCs Kc data
        ISETBioKc{4}.eccentricityDegs = cat(2, ISETBioKc{4}.eccentricityDegs, ISETBioKc{iDataSet}.eccentricityDegs); %cellTemporalEquivalentEccDegs(LcenterRGCindices);
        ISETBioKc{4}.data = cat(2, ISETBioKc{4}.data, ISETBioKc{iDataSet}.data);
    

        % The RsRcRatios as a function of eccentricity plot for LM-center RGCs
        ISETBioRsRcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioRsRcRatios{iDataSet}.data = cat(2, ISETBioRsRcRatios{iDataSet}.data, cellRsRcRatios(LMcenterRGCindices));
        % All-center RGCs RsRcRatios data
        ISETBioRsRcRatios{4}.eccentricityDegs = cat(2, ISETBioRsRcRatios{4}.eccentricityDegs, ISETBioRsRcRatios{iDataSet}.eccentricityDegs);
        ISETBioRsRcRatios{4}.data = cat(2, ISETBioRsRcRatios{4}.data, ISETBioRsRcRatios{iDataSet}.data);
   
        % The KsRcRatios as a function of eccentricity for LM-center RGCs
        ISETBioKsKcRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioKsKcRatios{iDataSet}.data = cat(2, ISETBioKsKcRatios{iDataSet}.data, cellKsKcRatios(LMcenterRGCindices));
        % All-center RGCs KsRcRatios data
        ISETBioKsKcRatios{4}.eccentricityDegs = cat(2, ISETBioKsKcRatios{4}.eccentricityDegs, ISETBioKsKcRatios{iDataSet}.eccentricityDegs);
        ISETBioKsKcRatios{4}.data = cat(2, ISETBioKsKcRatios{4}.data, ISETBioKsKcRatios{iDataSet}.data);
    
        % The S/C int ratio as a function of eccentricity for LM-center RGCs
        ISETBioSCintSensRatios{iDataSet}.eccentricityDegs = ISETBioRcDegs{iDataSet}.eccentricityDegs;
        ISETBioSCintSensRatios{iDataSet}.data = cat(2,ISETBioSCintSensRatios{iDataSet}.data, cellSCintSensRatios(LMcenterRGCindices));
        % All-center RGCs SCintSensRatios data
        ISETBioSCintSensRatios{4}.eccentricityDegs = cat(2, ISETBioSCintSensRatios{4}.eccentricityDegs, ISETBioSCintSensRatios{iDataSet}.eccentricityDegs);
        ISETBioSCintSensRatios{4}.data = cat(2,ISETBioSCintSensRatios{4}.data, ISETBioSCintSensRatios{iDataSet}.data);
       
    end % iMosaic

    generateAllSynthesizedMosaicOutlinesFigure(mosaicOutlines, false);
    generateAllSynthesizedMosaicOutlinesFigure(mosaicOutlines, true);

    % Transform KsKc ratios to log10(ks/Kc)
    CK95KsKcRatiosLog = CK95KsKcRatios;
    ISETBioKsKcRatiosLog = ISETBioKsKcRatios;
    CK95KsKcRatiosLog.data = log10(CK95KsKcRatios.data);
    ISETBioKsKcRatiosLog{1}.data = log10(ISETBioKsKcRatios{1}.data);
    ISETBioKsKcRatiosLog{2}.data = log10(ISETBioKsKcRatios{2}.data);
    ISETBioKsKcRatiosLog{3}.data = log10(ISETBioKsKcRatios{3}.data);
    ISETBioKsKcRatiosLog{4}.data = log10(ISETBioKsKcRatios{4}.data);

    pdfFileNamePostfix{1} = 'LconeCenter';
    pdfFileNamePostfix{2} = 'MconeCenter';
    pdfFileNamePostfix{3} = 'mixedLMconeCenter';
    pdfFileNamePostfix{4} = 'midget'; %all center types';

    eccRange = [0.1 30];
    eccTicks = [0.1 0.3 1 3 10 30];
    ISETBioSetsColors = [235 119 129; 82 180 140; 170 160 100; 50 130 230]/255;
    scatterFaceAlpha = 0.02;
    histogramFaceAlpha = 0.3;


    for iDataSet = 1:4
        pngFileName = fullfile(pdfDirectory,sprintf('SummaryStats_%s.png',pdfFileNamePostfix{iDataSet}));
        
        hFig = figure(iDataSet); clf;
        ff = MSreadyPlot.figureFormat('2x4-tall');
        theAxes = MSreadyPlot.generateAxes(hFig,ff);
        set(hFig, 'Color', [1 1 1]);

        if (iDataSet == 1)
            theVisualizedDataSetsScatter = [1 2];
            theVisualizedDataSetsMeans = theVisualizedDataSetsScatter;
        elseif (iDataSet == 2)
            theVisualizedDataSetsScatter = [2,1];
            theVisualizedDataSetsMeans = theVisualizedDataSetsScatter;
        elseif (iDataSet == 3) || (iDataSet == 4)
            theVisualizedDataSetsScatter = iDataSet;
            theVisualizedDataSetsMeans = [1 2];
        end

        % The RcDegs as a function of eccentricity plot
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsAsAFunctionOfEcc(theAxes{1,1}, ...
                  CK95RcDegs, ISETBioRcDegs, theVisualizedDataSetsScatter, theVisualizedDataSetsMeans, ...
                  ISETBioSetsColors, scatterFaceAlpha, ...
                  eccRange, eccTicks, ...
                  [0.005 0.1], (0:0.01:0.1), ...
                  'linear', 'Rc (degs)', pdfFileNamePostfix{iDataSet}, ...
                  showZscoresInsteadOfData, ...
                  employTemporalEquivalentEccentricity, onlyShowCronerKaplan95Data, ...
                  ff);

        % The Kc vs RcDegs plot
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParams(theAxes{2,1}, ...
                  CK95KcVsRcDegs.RcDegs, CK95KcVsRcDegs.Kc, ISETBioRcDegs, ISETBioKc, ...
                  theVisualizedDataSetsScatter, ...
                  ISETBioSetsColors, scatterFaceAlpha, ...
                  [0.006 0.105], [0.006 0.01 0.03 0.06 0.1], ...
                  [10 1e4], [10 100 1e3 1e4], ...
                  'log', 'Rc (degs)', ...
                  'log', 'Kc', ...
                  pdfFileNamePostfix{iDataSet}, ...
                  ff);

        % The Rs/Rc ratio as a function of eccentricity plot
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsAsAFunctionOfEcc(theAxes{1,2}, ...
                  CK95RsRcRatios, ISETBioRsRcRatios, theVisualizedDataSetsScatter, theVisualizedDataSetsMeans, ...
                  ISETBioSetsColors, scatterFaceAlpha, ...
                  eccRange, eccTicks, ...
                  [0 22], (0:2:26), ...
                  'linear', 'Rs/Rc', '', ...
                  showZscoresInsteadOfData, onlyShowCronerKaplan95Data, ...
                  employTemporalEquivalentEccentricity, ...
                  ff);

        % The Rs/Rc ratio histograms
        histogramEdges = 0:0.5:20;
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsHistograms(theAxes{2,2}, ...
              CK95RsRcRatios, ISETBioRsRcRatios, theVisualizedDataSetsScatter, histogramEdges, ...
              ISETBioSetsColors, histogramFaceAlpha, ...
              [0 22], (0:2:30), ...
              'Rs/Rc', '', ...
              ff);


        % The Ks/Kc int ratio as a function of eccentricity plot
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsAsAFunctionOfEcc(theAxes{1,4}, ...
              CK95KsKcRatios, ISETBioKsKcRatios, theVisualizedDataSetsScatter, theVisualizedDataSetsMeans, ...
              ISETBioSetsColors, scatterFaceAlpha, ...
              eccRange, eccTicks, ...
              [1e-4 1e-1], [1e-4 1e-3 1e-2 1e-1], ...
              'log', 'Ks/Kc', '', ...
              showZscoresInsteadOfData, ...
              employTemporalEquivalentEccentricity, onlyShowCronerKaplan95Data, ...
              ff);


        % The Ks/Kc histograms
        histogramEdges = -4:0.1:0;
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsHistograms(theAxes{2,4}, ...
              CK95KsKcRatiosLog, ISETBioKsKcRatiosLog, theVisualizedDataSetsScatter, histogramEdges, ...
              ISETBioSetsColors, histogramFaceAlpha, ...
              [-4 -1], (-4:0.5:0), ...
              'log10(Ks/Kc)', '', ...
              ff);


        % The S/C int sens ratio as a function of eccentricity
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsAsAFunctionOfEcc(theAxes{1,3}, ...
              CK95SCintSensRatios, ISETBioSCintSensRatios, theVisualizedDataSetsScatter, theVisualizedDataSetsMeans, ...
              ISETBioSetsColors, scatterFaceAlpha, ...
              eccRange, eccTicks, ...
              [0 1.0], (0:0.1:1), ...
              'linear', 'intS_s/intS_c', '', ...
              showZscoresInsteadOfData, onlyShowCronerKaplan95Data, ...
              employTemporalEquivalentEccentricity, ...
              ff);

        % The S/c int ratio histograms
        histogramEdges = 0:0.02:1.0;
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsHistograms(theAxes{2,3}, ...
              CK95SCintSensRatios, ISETBioSCintSensRatios, theVisualizedDataSetsScatter, histogramEdges, ...
              ISETBioSetsColors, histogramFaceAlpha, ...
              [0.1 0.9], (0:0.1:1), ...
              'intS_s/intS_c', '', ...
              ff);

        if (showZscoresInsteadOfData)
            NicePlot.exportFigToPNG(pngFileName, hFig, 300);
        else
            pdfFileName = fullfile(pdfDirectory,sprintf('SummaryStats_%s.pdf',pdfFileNamePostfix{iDataSet}));
            NicePlot.exportFigToPNG(pdfFileName, hFig, 300); 
        end


        % ============ Now generate separate figs for each parameter (PLOS format)
        rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';

        % The RcDegs as a function of eccentricity plot
        hFig = MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsAsAFunctionOfEcc([], ...
                  CK95RcDegs, ISETBioRcDegs, theVisualizedDataSetsScatter, theVisualizedDataSetsMeans, ...
                  ISETBioSetsColors, scatterFaceAlpha, ...
                  eccRange, eccTicks, ...
                  [0.005 0.1], (0:0.01:0.1), ...
                  'linear', 'Rc (degs)', pdfFileNamePostfix{iDataSet}, ...
                  showZscoresInsteadOfData, onlyShowCronerKaplan95Data, ...
                  employTemporalEquivalentEccentricity, ...
                  []);


        % Export fig 
        PNGFullFileName = strrep(strrep(pngFileName, '.png', 'RcDegs.png'), pdfDirectory, rawFiguresRoot);
        
        if (showZscoresInsteadOfData)
            PDFFullFileName = strrep(PNGFullFileName, '.png', '.pdf');
            NicePlot.exportFigToPDF(PDFFullFileName, hFig, 300);
        else
            NicePlot.exportFigToPNG(PNGFullFileName, hFig, 300);
        end
        

        % The S/C int sens ratio as a function of eccentricity
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsAsAFunctionOfEcc([], ...
              CK95SCintSensRatios, ISETBioSCintSensRatios, theVisualizedDataSetsScatter, theVisualizedDataSetsMeans, ...
              ISETBioSetsColors, scatterFaceAlpha, ...
              eccRange, eccTicks, ...
              [0 1.0], (0:0.1:1), ...
              'linear', 'intS_s/intS_c', '', ...
              showZscoresInsteadOfData, onlyShowCronerKaplan95Data, ...
              employTemporalEquivalentEccentricity, ...
              []);

        % Export fig 
        PNGFullFileName = strrep(strrep(pngFileName, '.png', 'intSCRatio.png'), pdfDirectory, rawFiguresRoot);

        if (showZscoresInsteadOfData)
            PDFFullFileName = strrep(PNGFullFileName, '.png', '.pdf');
            NicePlot.exportFigToPDF(PDFFullFileName, hFig, 300);
        else
            NicePlot.exportFigToPNG(PNGFullFileName, hFig, 300);
        end
        

        % The Ks/Kc int ratio as a function of eccentricity plot
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsAsAFunctionOfEcc([], ...
              CK95KsKcRatios, ISETBioKsKcRatios, theVisualizedDataSetsScatter, theVisualizedDataSetsMeans, ...
              ISETBioSetsColors, scatterFaceAlpha, ...
              eccRange, eccTicks, ...
              [1e-4 1e-1], [1e-4 1e-3 1e-2 1e-1], ...
              'log', 'Ks/Kc', '', ...
              showZscoresInsteadOfData, onlyShowCronerKaplan95Data, ...
              employTemporalEquivalentEccentricity, ...
              []);

        % Export fig 
        PNGFullFileName = strrep(strrep(pngFileName, '.png', 'KsKcRatio.png'), pdfDirectory, rawFiguresRoot);
        

        if (showZscoresInsteadOfData)
            PDFFullFileName = strrep(PNGFullFileName, '.png', '.pdf');
            NicePlot.exportFigToPDF(PDFFullFileName, hFig, 300);
        else
            NicePlot.exportFigToPNG(PNGFullFileName, hFig, 300);
        end


        % The Rs/Rc ratio as a function of eccentricity plot
        MSreadyPlot.renderAchievedVsCronerKaplanDoGmodelParamsAsAFunctionOfEcc([], ...
                  CK95RsRcRatios, ISETBioRsRcRatios, theVisualizedDataSetsScatter, theVisualizedDataSetsMeans, ...
                  ISETBioSetsColors, scatterFaceAlpha, ...
                  eccRange, eccTicks, ...
                  [0 22], (0:2:26), ...
                  'linear', 'Rs/Rc', '', ...
                  showZscoresInsteadOfData, onlyShowCronerKaplan95Data, ...
                  employTemporalEquivalentEccentricity, ...
                  []);

        % Export fig 
        PNGFullFileName = strrep(strrep(pngFileName, '.png', 'RsRcRatio.png'), pdfDirectory, rawFiguresRoot);
        
        if (showZscoresInsteadOfData)
            PDFFullFileName = strrep(PNGFullFileName, '.png', '.pdf');
            NicePlot.exportFigToPDF(PDFFullFileName, hFig, 300);
        else
            NicePlot.exportFigToPNG(PNGFullFileName, hFig, 300);
        end
    end % iDataSet

    % Generate paper-ready figures (scaled versions of the figures i
    % nrawFiguresRoot directory) which are stored in the PaperReady folder
    PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
    commandString = sprintf('%s/cpdf -args %s/generatePLOSOnePaperReadyFigures.txt', PLOSdirectory, PLOSdirectory);
    system(commandString);

end


function generateAllSynthesizedMosaicOutlinesFigure(mosaicOutlines, employTEE)

    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x1 small');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    hold(theAxes{1,1}, 'on');

    mosaicColors = brewermap(numel(mosaicOutlines), 'Pastel1');

    for iMosaic = 1:numel(mosaicOutlines)
        outlineStruct = mosaicOutlines{iMosaic}.outlineStruct;
        if (employTEE)
            xo = abs(outlineStruct.temporalEquivalentEccDegs(1));
            yo = abs(outlineStruct.temporalEquivalentEccDegs(2));
        else
            xo = outlineStruct.eccDegs(1);
            yo = outlineStruct.eccDegs(2);
        end

        xOutline = xo + outlineStruct.sizeDegs(1)*0.5*[-1 1 1 -1 -1];
        if (employTEE)
            xOutline = max(xOutline, 0.05);
        end

        yOutline = yo + outlineStruct.sizeDegs(2)*0.5*[-1 -1 1 1 -1];
        patch(theAxes{1,1}, xOutline, yOutline, mosaicColors(iMosaic,:), 'FaceAlpha', 0.5);
        plot(theAxes{1,1}, xOutline, yOutline, '-', 'LineWidth', ff.lineWidth, 'Color', 0.3*mosaicColors(iMosaic,:));
        
    end

    set(theAxes{1,1}, 'XTick', -20:5:20, 'YTick', -10:2:10);
    axis(theAxes{1,1}, 'equal');

    % Limits and ticks
    if (employTEE)
        set(theAxes{1,1}, 'XLim', [0.05 30], 'XTick', [0.1 0.3 1 3 10 30], 'XScale', 'log', 'YLim', [-2.5 2.5]);
    else
        set(theAxes{1,1}, 'XLim', [-20 20], 'YLim', [-2.5 2.5]);
    end

    

    % Grids
    grid(theAxes{1,1}, 'on');
    box(theAxes{1,1}, 'off');

    if (employTEE)
        xlabel(theAxes{1,1},'temporal equivalent eccentricity (degs)', 'FontAngle', ff.axisFontAngle);
    else
        xlabel(theAxes{1,1},'eccentricity (degs)', 'FontAngle', ff.axisFontAngle);
    end

    ylabel(theAxes{1,1},'eccentricity (degs)', 'FontAngle', ff.axisFontAngle);

    % axis color and width
    set(theAxes{1,1}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    % Font size
    set(theAxes{1,1}, 'FontSize', ff.fontSize);


    rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';
    if (employTEE)
        pdfFileName = sprintf('%s/synthesizedMosaicOutlinesTEE.pdf', rawFiguresRoot);
    else
        pdfFileName = sprintf('%s/synthesizedMosaicOutlines.pdf', rawFiguresRoot);
    end

    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
   

    % Generate paper-ready figures (scaled versions of the figures i
    % nrawFiguresRoot directory) which are stored in the PaperReady folder
    PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
    commandString = sprintf('%s/cpdf -args %s/generatePLOSOnePaperReadyFigures.txt', PLOSdirectory, PLOSdirectory);
    system(commandString);

end

function [cellEccDegs, cellCenterConeTypeWeights, ...
          cellRcDegs, cellKc, cellRsRcRatios, cellKsKcRatios, ...
          cellSCintSensRatios] = extractComputeReadyMosaicData(...
                                    theComputeReadyMRGCmosaic, ...
                                    theMRGCMosaicOptimalSTFs, ...
                                    employTemporalEquivalentEccentricity)

    rgcIndex = 1;
    theSTFdata = theMRGCMosaicOptimalSTFs{rgcIndex};

    idxKc = find(strcmp(theSTFdata.DoGfitParams.names, 'Kc'));
    idxRcDegs = find(strcmp(theSTFdata.DoGfitParams.names, 'RcDegs'));
    idxRsRcRatio = find(strcmp(theSTFdata.DoGfitParams.names, 'RsToRc'));
    idxKsKcRatio = find(strcmp(theSTFdata.DoGfitParams.names, 'kS/kC'));
    
    cellEccDegs = zeros(1,theComputeReadyMRGCmosaic.rgcsNum);
    cellCenterConeTypeWeights = zeros(3,theComputeReadyMRGCmosaic.rgcsNum);
    cellKc = cellEccDegs;
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

        cellKc(iRGC) = theSTFdata.DoGfitParams.finalValues(idxKc);
        cellRcDegs(iRGC) = theSTFdata.DoGfitParams.finalValues(idxRcDegs);
        cellRsRcRatios(iRGC) = theSTFdata.DoGfitParams.finalValues(idxRsRcRatio);
        cellKsKcRatios(iRGC) = theSTFdata.DoGfitParams.finalValues(idxKsKcRatio);
        
    end

    cellSCintSensRatios = cellKsKcRatios .* cellRsRcRatios.^2;
end


function [CK95RcDegs, CK95KcVsRcDegs, CK95RsRcRatios, CK95KsKcRatios, CK95SCintSensRatios] = ...
    loadCronerKaplan95ParvoCellularData()
    
    [CK95RcDegs.eccentricityDegs, CK95RcDegs.data] = ...
        RGCmodels.CronerKaplan.digitizedData.parvoCenterRadiusAgainstEccentricity();

    CK95KcVsRcDegs = ...
        RGCmodels.CronerKaplan.digitizedData.parvoCenterKcVsRadiusDegs();

    [CK95RsRcRatios.eccentricityDegs, CK95RsRcRatios.data] = ...
        RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();
    CK95RsRcRatios.data = 1./CK95RsRcRatios.data;

    [CK95KsKcRatios.eccentricityDegs, CK95KsKcRatios.data] = ...
        RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterPeakSensisitivityRatioAgainstEccentricity();

    [CK95SCintSensRatios.eccentricityDegs, CK95SCintSensRatios.data] = ...
       RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();
end
