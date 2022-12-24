function analyzeEccentricityComputations()

    mosaicEccDegs = [...
        0 0; ...
        -1 0; ...
        -2 0; ...
        -4 0; ...
        -6 0; ...
        -8 0]; 

    mosaicEccDegs = [...
        0 0; ...
        -1 0];

    inspectTheSpatialRFs = ~true;
    inspectTheSTFs = ~true;
    contrastModelToCronerAndKaplan = true;

    % Αctions
    centerMostRGCsNumToAnalyze = input('Enter # of center-most cells to include in the analysis. [] for all: ');

    for H1index = 1:4
        doIt(mosaicEccDegs, H1index, centerMostRGCsNumToAnalyze, ...
        inspectTheSpatialRFs, inspectTheSTFs, contrastModelToCronerAndKaplan);
    end

end

function doIt(mosaicEccDegs, H1cellIndex, centerMostRGCsNumToAnalyze, ...
    inspectTheSpatialRFs, inspectTheSTFs, contrastModelToCronerAndKaplan)

   % Get dropboxDir & intermediate data files location
    computerInfo = GetComputerInfo();
    switch (computerInfo.localHostName)
        case 'Ithaka'
            dropboxDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris/midgetRGCMosaics';

        case 'Crete'
            dropboxDir = '/Volumes/Dropbox/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris/midgetRGCMosaics';

        otherwise
            if (contains(computerInfo.networkName, 'leviathan'))
                dropboxDir = '/media/dropbox_disk/Aguirre-Brainard Lab Dropbox/isetbio isetbio/midgetRGCMosaics';
            else
                error('Could not establish dropbox location')
            end
    end
    mappedRFsDir = sprintf('%s/RGCMosaicsWithFixedParamsH%d', dropboxDir, H1cellIndex);

    ZernikeDataBase = 'Polans2015';
    subjectRankOrder = 6;
    pupilDiameterMM = 3.0;


    % STF stimulus parameters
    % L+M contrast for gratings used to measure the STFs 
    coneContrasts = [1 1 0];
    
    

    if (inspectTheSTFs)
        for iEcc = 1:size(mosaicEccDegs,1)
            inspectSTFs(mosaicEccDegs(iEcc,:), mappedRFsDir, ...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                coneContrasts, centerMostRGCsNumToAnalyze);
        end
    end

   
    if (inspectTheSpatialRFs)
        coVisualizeMeasuredSTFs = ~true;
        coVisualizeFittedSTFs = true;

        for iEcc = 1:size(mosaicEccDegs,1)
            inspectSpatialRFs(mosaicEccDegs(iEcc,:), mappedRFsDir, ...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                coneContrasts, centerMostRGCsNumToAnalyze, ...
                coVisualizeFittedSTFs, coVisualizeMeasuredSTFs);
        end
    end
    

    
    if (contrastModelToCronerAndKaplan)
        for iEcc = 1:size(mosaicEccDegs,1)
            dataOut{iEcc} = extractPatchDoGParams(mosaicEccDegs(iEcc,:), mappedRFsDir, ...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                coneContrasts, centerMostRGCsNumToAnalyze);
        end
        contastDerivedDoGparamsToCronerAndKaplanDoGParams(mosaicEccDegs, dataOut, mappedRFsDir);
    end
end


function contastDerivedDoGparamsToCronerAndKaplanDoGParams(mosaicEccDegs, dataOut, mappedRFsDir)

    mosaicTemporalEccDegs = [];
    modelRGCRcDegs = [];
    modelRGCtemporalEccDegs = [];
    modelRGCRsRcRatios = [];
    modelRGCKsKcRatios = [];
    modelRGCSCintSensRatios = [];
    targetRsRcRatios = [];
    targetSCintSensRatios = [];

    for iEcc = 1:size(mosaicEccDegs,1)
        modelRGCRcDegs = cat(2, modelRGCRcDegs, dataOut{iEcc}.modelRGCRcDegs);
        modelRGCtemporalEccDegs = cat(2, modelRGCtemporalEccDegs, dataOut{iEcc}.modelRGCtemporalEccDegs);
        modelRGCRsRcRatios = cat(2, modelRGCRsRcRatios, dataOut{iEcc}.modelRGCRsRcRatios);
        modelRGCKsKcRatios = cat(2, modelRGCKsKcRatios, dataOut{iEcc}.modelRGCKsKcRatios);
        modelRGCSCintSensRatios = cat(2, modelRGCSCintSensRatios, dataOut{iEcc}.modelRGCSCintSensRatios);
        targetRsRcRatios = cat(2, targetRsRcRatios, dataOut{iEcc}.targetRsRcRatio);
        targetSCintSensRatios = cat(2, targetSCintSensRatios, dataOut{iEcc}.targetSCintSensRatio);
        mosaicTemporalEccDegs = cat(2, mosaicTemporalEccDegs, dataOut{iEcc}.mosaicTemporalEccDegs);
    end

    mosaicTemporalEccDegs(mosaicTemporalEccDegs<0.1) = 0.1;
    modelRGCtemporalEccDegs(modelRGCtemporalEccDegs<0.1) = 0.1;

    % Sort according to ecc
    [~, idx] = sort(mosaicTemporalEccDegs, 'ascend');
    mosaicTemporalEccDegs = mosaicTemporalEccDegs(idx);

    hFig = figure(4); clf;
    set(hFig, 'Position', [10 10 1200 1100], 'Color', [1 1 1]);


    % Top-left: RcDegs 
    panelPositionX = 0.08;
    panelPositionY = 0.67;
    ax = subplot('Position', [panelPositionX panelPositionY 0.37 0.3]); 

    % Retrieve the C&K data
    [CronerKaplanTemporalEccDegs, CronerKaplanRcDegs] = ...
            RGCmodels.CronerKaplan.digitizedData.parvoCenterRadiusAgainstEccentricity();
    
    % Plot model against C&K
    plotModelAgainstCronerKaplan(ax, modelRGCtemporalEccDegs, modelRGCRcDegs*60, ...
        CronerKaplanTemporalEccDegs, CronerKaplanRcDegs*60, [0 8], 0:1:14, sprintf('%d\n', 0:2:14), 'linear', 'Rc (arc min)')

  
    % Compute the Z-scores
    [eccDegsTable, zScores] = computeModelDataZscore(mosaicTemporalEccDegs, ...
        CronerKaplanTemporalEccDegs, CronerKaplanRcDegs, ...
        modelRGCtemporalEccDegs, modelRGCRcDegs);

    % Plot the Z-scores
    ax = subplot('Position', [panelPositionX panelPositionY-0.1 0.37 0.08]); 
    plotZScores(ax, eccDegsTable, zScores, true);


    % Top - right: S/C int sensitivity
    panelPositionX = 0.60;
    panelPositionY = 0.67;
    ax = subplot('Position', [panelPositionX panelPositionY 0.37 0.3]); 

    % Retrieve the C&K data
    [CronerKaplanTemporalEccDegs, CronerKaplanSCintSensRatios] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();

    % Plot model against C&K
    plotModelAgainstCronerKaplan(ax, modelRGCtemporalEccDegs, modelRGCSCintSensRatios, ...
        CronerKaplanTemporalEccDegs, CronerKaplanSCintSensRatios, [0 1], 0:0.2:1.0, sprintf('%2.1f\n', 0:0.2:1.0), 'linear', 'S/C int. sensitivity');


    % Compute the Z-scores
    [eccDegsTable, zScores] = computeModelDataZscore(mosaicTemporalEccDegs, ...
        CronerKaplanTemporalEccDegs, CronerKaplanSCintSensRatios, ...
        modelRGCtemporalEccDegs, modelRGCSCintSensRatios);

    % Plot the Z-scores
    ax = subplot('Position', [panelPositionX panelPositionY-0.1 0.37 0.08]); 
    plotZScores(ax, eccDegsTable, zScores, true);


    % Bottom-left: Rs/Rc ratio
    panelPositionX = 0.08;
    panelPositionY = 0.17;
    ax = subplot('Position', [panelPositionX panelPositionY 0.37 0.3]); 

    % Retrieve the C&K data
    [CronerKaplanTemporalEccDegs, CronerKaplanRcRsRatios] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();

    % Plot model against C&K
    plotModelAgainstCronerKaplan(ax, modelRGCtemporalEccDegs, modelRGCRsRcRatios, ...
        CronerKaplanTemporalEccDegs, 1./CronerKaplanRcRsRatios, [0 50], 0:5:50, sprintf('%2.0f\n', 0:5:50), 'linear', 'Rs/Rc');


    % Compute the Z-scores
    [eccDegsTable, zScores] = computeModelDataZscore(mosaicTemporalEccDegs, ...
        CronerKaplanTemporalEccDegs, 1./CronerKaplanRcRsRatios, ...
        modelRGCtemporalEccDegs, modelRGCRsRcRatios);

    % Plot the Z-scores
    ax = subplot('Position', [panelPositionX panelPositionY-0.1 0.37 0.08]); 
    plotZScores(ax, eccDegsTable, zScores, false);

    
  
    % Bottom-right: Ks/Kc ratio
    panelPositionX = 0.60;
    panelPositionY = 0.17;
    ax = subplot('Position', [panelPositionX panelPositionY 0.37 0.3]); 

    % Retrieve the C&K data
    [CronerKaplanTemporalEccDegs, CronerKaplanKsKcRatios] = ...
            RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterPeakSensisitivityRatioAgainstEccentricity();

    % Plot model against C&K
    plotModelAgainstCronerKaplan(ax, modelRGCtemporalEccDegs, modelRGCKsKcRatios, ...
        CronerKaplanTemporalEccDegs, CronerKaplanKsKcRatios, [1e-4 1], [1e-4 1e-3 1e-2 1e-1 1], {'', '1e-3', '1e-2', '1e-1', '1'}, ...
         'log', 'Ks/Kc');


    % Compute the Z-scores
    [eccDegsTable, zScores] = computeModelDataZscore(mosaicTemporalEccDegs, ...
        CronerKaplanTemporalEccDegs, CronerKaplanKsKcRatios, ...
        modelRGCtemporalEccDegs, modelRGCKsKcRatios);

    % Plot the Z-scores
    ax = subplot('Position', [panelPositionX panelPositionY-0.1 0.37 0.08]); 
    plotZScores(ax, eccDegsTable, zScores, false);
   

    pdfFileName = fullfile(mappedRFsDir, 'ComparisonToCronerKaplan.pdf');
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
    fprintf('PDF saved at %s\n', pdfFileName);
end


function plotModelAgainstCronerKaplan(ax, modelRGCtemporalEccDegs, modelParam, ...
        CronerKaplanTemporalEccDegs, CronerKaplanParam, YLims, YTicks, YTickLabel, YScale, plotTitle)

    p1 = scatter(ax,modelRGCtemporalEccDegs, modelParam, 100, ...
             'LineWidth', 1.0, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor',[1 0.5 0.5]);
    hold(ax, 'on');
    p2 = scatter(ax,CronerKaplanTemporalEccDegs, CronerKaplanParam, 144, ...
             'LineWidth', 1.0, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[0.7 0.7 0.7]);
    grid(ax, 'on')
    set(ax, 'XLim', [0.1 30], 'XScale', 'log', 'XTick', [0.1 0.3 1 3 10 30], 'XTickLabel', {}, ...
             'YLim', YLims, 'YTick', YTicks, 'YTickLabel', YTickLabel, 'YScale', YScale', 'FontSize', 20);
    set(ax, 'TickDir', 'in', 'YColor', [0 0 0]);
    
    ylabel(ax,plotTitle);

    legend(ax, [p1 p2], {'ISETBio', 'macaque (Croner & Kaplan)',}, ...
        'NumColumns', 1, 'Location', 'NorthWest', 'box', 'off', 'FontSize', 15);
end


function plotZScores(ax, eccDegsTable, zScores, noXLabel)
    validZscoreIndices = find((~isnan(zScores))&(~isinf(zScores)));

    hold(ax, 'on');
    xx = [0.1 1:30];
    y1 = xx*0+1;
    y2 = xx*0-1;
    shadedAreaBetweenTwoLines(ax, xx, y1, y2, [1 0.5 0.5], [0.7 0.0 0], 0.3, 1.0, 'none');
    stem(ax, eccDegsTable, zScores, 'o', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1.0 0.2 0.2]*0.95, ...
        'MarkerSize', 14, 'LineWidth', 1.5, 'Color', [1.0 0.2 0.2]*0.95);
    plot(ax, eccDegsTable(validZscoreIndices), zScores(validZscoreIndices), 'o', 'MarkerFaceColor', [1.0 0.5 0.5], 'MarkerEdgeColor', [1.0 0.2 0.2]*0.95, ...
        'MarkerSize', 14, 'LineWidth', 1.5, 'Color', [ 1.0 0.2 0.2]*0.95);

    zScoreText = sprintf('Z($\\mu,\\sigma$) = %2.2f, %2.2f', ...
        mean(zScores(validZscoreIndices)), std(zScores(validZscoreIndices)));
    text(ax, 0.11,-3.5, zScoreText, 'Color', [0.5 0 0], 'Interpreter', 'latex', 'FontSize', 16);
    grid(ax, 'on')
    set(ax, 'XLim', [0.1 30], 'XScale', 'log', 'XTick', [0.1 0.3 1 3 10 30]);
    set(ax, 'YLim', [-4 2], 'YTick', -4:1:4, 'FontSize', 20, 'YTickLabel', {'-4', '', '-2', '', '0', '', '+2', '', '+4'});
   
    set(ax, 'TickDir', 'both', 'YColor', [1 0 0.0]);
    ylabel(ax,sprintf('Zscore\n(signed)'));
    if (noXLabel)
    else
        xlabel(ax,'temporal equivalent eccentericity (degs)');
    end
end

function [eccRanges,  zScores] = computeModelDataZscore(...
         mosaicTemporalEccDegs, ...
         eccDegsData, parameterData, ...
         eccDegsModel, parameterModel)

    eccRanges = mosaicTemporalEccDegs;

    zScores = nan(1,numel(eccRanges));
    dataMean = nan(1,numel(eccRanges));
    modelMean = nan(1,numel(eccRanges)); 
    for iEcc = 1:numel(eccRanges)
        centerEcc = eccRanges(iEcc);
        if (centerEcc < 2)
            eccBinWidthDegs = 0.5;
        elseif (centerEcc < 4)
            eccBinWidthDegs = 0.5;
        elseif (centerEcc < 8)
            eccBinWidthDegs = 1.0;
        elseif (centerEcc < 16)
            eccBinWidthDegs = 2.0;
        else
            eccBinWidthDegs = 2.5;
        end

        eccBinWidthDegs = 0.75*eccBinWidthDegs;

        idx = find(abs(eccDegsData-centerEcc) <= eccBinWidthDegs);
        idx2 = find(abs(eccDegsModel-centerEcc) <= eccBinWidthDegs);
        if (~isempty(idx)) && (~isempty(idx2))
            weights = 1 - abs(centerEcc-eccDegsData(idx))/eccBinWidthDegs;
            dataMean(iEcc) = sum(parameterData(idx).*weights)/sum(weights);
            dataStd = std(parameterData(idx), weights);

            weights = 1 - abs(centerEcc-eccDegsModel(idx2))/eccBinWidthDegs;
            modelMean(iEcc) = sum(parameterModel(idx2).*weights)/sum(weights);
            modelStd = std(parameterModel(idx2), weights);

            signedDistance =  modelMean(iEcc) - dataMean(iEcc);
            meanStd = 0.5*(dataStd+modelStd);
            
            zScores(iEcc) = signedDistance/dataStd;
        end
    end
end


function dataOut = extractPatchDoGParams(mosaicEccDegs, mappedRFsDir, ZernikeDataBase, subjectRankOrder, ...
    pupilDiameterMM, coneContrasts, centerMostRGCsNumToAnalyze)

    
    % Load the midgetRGCMosaic
    fName = fullfile(mappedRFsDir, sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2)));
    load(fName, 'theMidgetRGCmosaic');

    % Assemble the responses filename
    responsesFileNamePostFix = sprintf('_STFresponses_cLMS_%2.2f_%2.2f_%2.2f.mat', coneContrasts(1), coneContrasts(2), coneContrasts(3));
    fNameResponses = strrep(fName, '.mat', responsesFileNamePostFix);
    load(fNameResponses, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'fittedSTFs');

    modelRGCRsRcRatios = nan(1, numel(fittedSTFs));
    modelRGCKsKcRatios = nan(1, numel(fittedSTFs));
    modelRGCSCintSensRatios = nan(1, numel(fittedSTFs));
    modelRGCRcDegs = nan(1, numel(fittedSTFs));
    modelRGCtemporalEccDegs = nan(1, numel(fittedSTFs));

    if (~isempty(centerMostRGCsNumToAnalyze))
        fprintf('There are %d RGCs in the patch. Analyzing the %d center-most.\n', numel(fittedSTFs),centerMostRGCsNumToAnalyze);
        fittedSTFs = fittedSTFs(1:centerMostRGCsNumToAnalyze);
    end


    for iRGC = 1:numel(fittedSTFs)
        f = fittedSTFs{iRGC};
        theRGCindex = f.targetRGC;
        connectivityVector = full(squeeze(theMidgetRGCmosaic.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
        inputConeIndices = find(abs(connectivityVector) > 0.0001);
        inputConeTypes = theMidgetRGCmosaic.inputConeMosaic.coneTypes(inputConeIndices);

        targetVisualParams = f.targetVisualRFDoGparams;    
        theFittedSTFDoGparams = f.theFittedSTFDoGparams;

        % Rs/Rc
        idx = find(strcmp(theFittedSTFDoGparams.names, 'RsToRc'));
        modelRsRcRatio = theFittedSTFDoGparams.finalValues(idx);

        % Ks/Kc
        idx = find(strcmp(theFittedSTFDoGparams.names, 'kS/kC'));
        modelKsKcRatio = theFittedSTFDoGparams.finalValues(idx);
        modelRGCKsKcRatios(iRGC) = modelKsKcRatio;

        % Rc degs
        idx = find(strcmp(theFittedSTFDoGparams.names, 'RcDegs'));
        modelRGCRcDegs(iRGC) = theFittedSTFDoGparams.finalValues(idx);
        
        % RsRc ratio
        modelRGCRsRcRatios(iRGC) = modelRsRcRatio;

        % S/C int sens ratio
        modelRGCSCintSensRatios(iRGC) = (modelRsRcRatio)^2 * modelKsKcRatio;

        % Temporal equivalent eccentricity
        rgcRFposDegs = theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,:);
        temporalEccDegs = theMidgetRGCmosaic.temporalEquivalentEccentricityForEccentricity(rgcRFposDegs);
        modelRGCtemporalEccDegs(iRGC) = sqrt(sum(temporalEccDegs.^2,2));
    end

    % Mosaic temporal equivalent eccentricity
    temporalEccDegs = theMidgetRGCmosaic.temporalEquivalentEccentricityForEccentricity(mosaicEccDegs);
    mosaicTemporalEccDegs = sqrt(sum(temporalEccDegs.^2,2));

    dataOut = struct(...
        'targetRsRcRatio', targetVisualParams.surroundToCenterRcRatio, ...
        'targetSCintSensRatio', targetVisualParams.surroundToCenterIntegratedSensitivityRatio, ...
        'modelRGCRcDegs', modelRGCRcDegs, ...
        'modelRGCRsRcRatios', modelRGCRsRcRatios, ...
        'modelRGCKsKcRatios', modelRGCKsKcRatios, ...
        'modelRGCSCintSensRatios', modelRGCSCintSensRatios , ...
        'modelRGCtemporalEccDegs', modelRGCtemporalEccDegs, ...
        'mosaicTemporalEccDegs', mosaicTemporalEccDegs ...
        );
end


function inspectSpatialRFs(mosaicEccDegs, mappedRFsDir, ZernikeDataBase, subjectRankOrder, ...
    pupilDiameterMM, coneContrasts, centerMostRGCsNumToAnalyze, coVisualizeFittedSTFs, coVisualizeMeasuredSTFs)

    % Load the midgetRGCMosaic
    fName = fullfile(mappedRFsDir, sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2)));
    load(fName, 'theMidgetRGCmosaic');

    thePSFData = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{1}.theVlambdaWeightedPSFData;
    maxVisualizedRFs = centerMostRGCsNumToAnalyze;
    videoFileName = strrep(fName, 'mat', '_SpatialRFs.mp4');

    
    if (coVisualizeMeasuredSTFs == false)
        theMidgetRGCmosaic.visualizeSpatialRFs(...
                'maxVisualizedRFs', maxVisualizedRFs, ...
                'generateVideo', true, ...
                'videoFileName', videoFileName);
    else

        % Assemble the responses filename
        responsesFileNamePostFix = sprintf('_STFresponses_cLMS_%2.2f_%2.2f_%2.2f.mat', coneContrasts(1), coneContrasts(2), coneContrasts(3));
        fNameResponses = strrep(fName, '.mat', responsesFileNamePostFix);
        load(fNameResponses, 'theMidgetRGCMosaicResponses', ...
             'orientationsTested', 'spatialFrequenciesTested', ...
             'spatialPhasesDegs', 'fittedSTFs');

        videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();

        % Only analyze the centerMostRGCsNumToAnalyze
        fittedSTFs = fittedSTFs(1:min([numel(fittedSTFs) centerMostRGCsNumToAnalyze]));

        for iRGC = 1:numel(fittedSTFs)
            if (iRGC > maxVisualizedRFs)
                continue;
            end

            f = fittedSTFs{iRGC};
            theRGCindex = f.targetRGC;
            connectivityVector = full(squeeze(theMidgetRGCmosaic.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
            inputConeIndices = find(abs(connectivityVector) > 0.0001);
            inputConeTypes = theMidgetRGCmosaic.inputConeMosaic.coneTypes(inputConeIndices);
    
            theFittedSTF = f.theFittedSTF;
            targetVisualParams = f.targetVisualRFDoGparams;
            
            theTargetSTFsupport = f.theTargetSTFdata(:,1);
            theTargetSTFDoGModelApproximation = f.theTargetSTFdata(:,2);
            theTargetSTFmeasured = f.theTargetSTFdata(:,3);        
            theFittedSTFDoGparams = f.theFittedSTFDoGparams;
    
    
            [hFig, allAxes] = theMidgetRGCmosaic.visualizeSpatialRFs(...
                'onlyForRGCwithIndex', theRGCindex, ...
                'generateVideo', false, ...
                'withEccentricityCrossHairs', true, ...
                'withPSFData', thePSFData, ...
                'fontSize', 16);
    
            % replace graphic is (1,1) with the STFs
            noXLabel = false; noYLabel = false;
            cla(allAxes{1,1});
            visualizeSTFs(allAxes{1,1}, spatialFrequenciesTested, ...
                        f.allMeasuredSTFs, f.theMeasuredSTFtoFit, ...
                        theFittedSTF.sfHiRes, theFittedSTF.compositeSTFHiRes, ...
                        theTargetSTFsupport, theTargetSTFmeasured, ...
                        noXLabel, noYLabel, inputConeTypes, theRGCindex, iRGC, numel(fittedSTFs));
            set(allAxes{1,1}, 'FontSize', 16);
            drawnow;
            videoOBJ.writeVideo(getframe(hFig));
    
        end

        videoOBJ.close();
        fprintf('spatial RFs video saved at %s\n', videoFileName );
    end

end


function inspectSTFs(mosaicEccDegs, mappedRFsDir, ZernikeDataBase, subjectRankOrder, ...
    pupilDiameterMM, coneContrasts, centerMostRGCsNumToAnalyze)

    % Load the midgetRGCMosaic
    fName = fullfile(mappedRFsDir, sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2)));
    load(fName, 'theMidgetRGCmosaic');

    % Assemble the responses filename
    responsesFileNamePostFix = sprintf('_STFresponses_cLMS_%2.2f_%2.2f_%2.2f.mat', coneContrasts(1), coneContrasts(2), coneContrasts(3));
    fNameResponses = strrep(fName, '.mat', responsesFileNamePostFix);
    load(fNameResponses, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'fittedSTFs');
    
    
    videoFileName = strrep(fName, 'mat', '_STFs.mp4');

    videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1500 450], 'Color', [1 1 1]);

    % Only analyze the centerMostRGCsNumToAnalyze
    fittedSTFs = fittedSTFs(1:min([numel(fittedSTFs) centerMostRGCsNumToAnalyze]));

    for iRGC = 1:numel(fittedSTFs)
        f = fittedSTFs{iRGC};
        theRGCindex = f.targetRGC;
        connectivityVector = full(squeeze(theMidgetRGCmosaic.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
        inputConeIndices = find(abs(connectivityVector) > 0.0001);
        inputConeTypes = theMidgetRGCmosaic.inputConeMosaic.coneTypes(inputConeIndices);

        theFittedSTF = f.theFittedSTF;
        targetVisualParams = f.targetVisualRFDoGparams;
        
        theTargetSTFsupport = f.theTargetSTFdata(:,1);
        theTargetSTFDoGModelApproximation = f.theTargetSTFdata(:,2);
        theTargetSTFmeasured = f.theTargetSTFdata(:,3);        
        theFittedSTFDoGparams = f.theFittedSTFDoGparams;
        
        % CronerKaplan targets
        % Rs/Rc ratio
        targetCronerKaplanSurroundToCenterRcRatio = RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio;

        % Temporal-equivalent eccentricity based SCint sensitivity ratio
        temporalEquivalentEccDegs = theMidgetRGCmosaic.temporalEquivalentEccentricityForEccentricity(theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,:));
        radialTemporalEquivalentEccDegs = sqrt(sum(temporalEquivalentEccDegs.^2,2));
        targetCronerKaplanSCIntSensitivity = RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(radialTemporalEquivalentEccDegs);
        
        ax = subplot('Position', [0.05 0.13 0.25 0.8]);
        noXLabel = false; noYLabel = false;
        visualizeSTFs(ax, spatialFrequenciesTested, ...
                f.allMeasuredSTFs, f.theMeasuredSTFtoFit, ...
                theFittedSTF.sfHiRes, theFittedSTF.compositeSTFHiRes, ...
                theTargetSTFsupport, theTargetSTFmeasured, ...
                noXLabel, noYLabel, inputConeTypes, theRGCindex, iRGC, numel(fittedSTFs));


        ax = subplot('Position', [0.39 0.13 0.25 0.8]);
        hold(ax, 'off')
        paramName = 'RsToRc';
        idxRsRc = find(strcmp(theFittedSTFDoGparams.names, paramName));
        visualizeParameter(ax, targetVisualParams.surroundToCenterRcRatio, ...
            theFittedSTFDoGparams.finalValues(idxRsRc), sprintf('Rs/Rc ratio (%d/%d)',iRGC, numel(fittedSTFs)),  [0 20], 2);
        allRsRcRatios(iRGC) = theFittedSTFDoGparams.finalValues(idxRsRc);
        [N,edges] = histcounts(allRsRcRatios,0:0.2:20);
        y1 = edges(1:end-1);
        shadedAreaBetweenTwoLines(ax, 20-N/max(N)*5, y1, y1*0+20, [0 0 1], [0 0 1], 0.3, 1.0, '-');
       
        % Targets (affected by CronerKaplan multipliers)
        p1 = plot(ax, targetVisualParams.surroundToCenterRcRatio*[1 1], [0 targetVisualParams.surroundToCenterRcRatio], 'r-', 'LineWidth', 1.0);
        plot(ax, [targetVisualParams.surroundToCenterRcRatio 20], targetVisualParams.surroundToCenterRcRatio*[1 1], 'r-', 'LineWidth', 1.0); 

        % Targets (CronerKaplan)
        p2 = plot(ax, targetCronerKaplanSurroundToCenterRcRatio*[1 1], [0 targetCronerKaplanSurroundToCenterRcRatio], 'k--', 'LineWidth', 1.0);
        plot(ax, [targetCronerKaplanSurroundToCenterRcRatio 20], targetCronerKaplanSurroundToCenterRcRatio*[1 1], 'k--', 'LineWidth', 1.0);
        
        legend(ax, [p1 p2], 'target', 'Croner&Kaplan mean', 'Location', 'NorthWest');


        ax = subplot('Position', [0.73 0.13 0.25 0.8]);
        hold(ax, 'off');
        idxKsKc = find(strcmp(theFittedSTFDoGparams.names, 'kS/kC'));
        theFittedSCintSensitivity = theFittedSTFDoGparams.finalValues(idxKsKc) * (theFittedSTFDoGparams.finalValues(idxRsRc))^2;
        visualizeParameter(ax, targetVisualParams.surroundToCenterIntegratedSensitivityRatio, ...
            theFittedSCintSensitivity, sprintf('S/C int. sensitivity ratio (%d/%d)',iRGC, numel(fittedSTFs)), [0 1], 0.1);
        allSCRatios(iRGC) = theFittedSCintSensitivity;
        [N,edges] = histcounts(allSCRatios,0:0.01:1);
        y1 = edges(1:end-1);
        shadedAreaBetweenTwoLines(ax, 1-N/max(N)*0.25, y1, y1*0 +1.0, [0 0 1], [0 0 1], 0.3, 1.0, '-');

        % Targets (affected by CronerKaplan multipliers)
        p1 = plot(ax, targetVisualParams.surroundToCenterIntegratedSensitivityRatio*[1 1], [0 targetVisualParams.surroundToCenterIntegratedSensitivityRatio], 'r-', 'LineWidth', 1.0);
        plot(ax, [targetVisualParams.surroundToCenterIntegratedSensitivityRatio 20], targetVisualParams.surroundToCenterIntegratedSensitivityRatio*[1 1], 'r-', 'LineWidth', 1.0); 

        % Targets (CronerKaplan)
        p2 = plot(ax, targetCronerKaplanSCIntSensitivity*[1 1], [0 targetCronerKaplanSCIntSensitivity], 'k--', 'LineWidth', 1.0);
        plot(ax, [targetCronerKaplanSCIntSensitivity 1], targetCronerKaplanSCIntSensitivity*[1 1], 'k--', 'LineWidth', 1.0);
        
        legend(ax, [p1 p2], 'target', 'Croner&Kaplan mean', 'Location', 'NorthWest');


        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
    end
    videoOBJ.close();

    fprintf('STFs video saved at %s\n', videoFileName );

end

function shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)

    idx = ~(isnan(y1) | (isnan(y2)));
    x = x(idx==true);
    y1 = y1(idx==true);
    y2 = y2(idx==true);

    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end

function visualizeParameter(ax, targetValue, achievedValue, paramName, paramRange, paramTick)
    plot(ax, paramRange, paramRange, 'k-', 'LineWidth', 1.0);
    hold(ax, 'on')
    scatter(ax,targetValue, achievedValue, 14*14, ...
        'MarkerFaceColor', [0. 0 1], 'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.25, 'LineWidth', 1.0);
    set(ax, 'XLim', paramRange, 'XTick', paramRange(1):paramTick:paramRange(2));
    set(ax, 'YLim', paramRange, 'YTick', paramRange(1):paramTick:paramRange(2));
    set(ax, 'XColor', [0.7 0 0], 'YColor', [0 0 0.7], 'LineWidth', 1.0);
    set(ax, 'TickDir', 'both', 'FontSize', 18);
    xtickangle(ax, 0);
    axis(ax, 'square');
    box(ax, 'off')
    grid(ax, 'on');
    title(ax, paramName, 'FontWeight', 'Normal');
    xlabel(ax,'target');
    ylabel(ax, 'achieved');
end


function visualizeSTFs(ax, measuredDpatialFrequencySupport, ...
    allMeasuredSTFs, theMeasuredSTFtoFit,  ...
    theFittedSTFsupport, theFittedSTF, ...
    theTargetSTFsupport, theTargetSTFmeasured, ...
    noXLabel, noYLabel, inputConeTypes, theRGCindex, iRGC, totalRGCs)

    
    plot(ax, measuredDpatialFrequencySupport, allMeasuredSTFs', ...
        '-', 'Color', [0.5 0.5 0.8], 'LineWidth', 1.5);
    hold(ax, 'on');
    
    p1 = plot(ax, theTargetSTFsupport, theTargetSTFmeasured, 'r-', 'LineWidth', 3.0);
    p2 = plot(ax, theFittedSTFsupport, theFittedSTF, '-', 'LineWidth', 3.0, 'Color', [0 0 1]);

    p3 = scatter(ax, measuredDpatialFrequencySupport, theMeasuredSTFtoFit, ...
        200, 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', [0 0 1], ...
        'MarkerFaceAlpha', 0.9, 'LineWidth', 1.5);

    legend(ax, [p1 p2 p3], {'target (DoG model)', 'achieved (DoG model fit)', 'achieved (measured)'}, ...
        'Location', 'SouthWest', 'box', 'off', 'FontSize', 16, 'NumColumns', 1);
    set(ax, 'XLim', [0.1 100], 'XScale', 'log', 'XTick', [0.1 0.3 1 3 10 30 100]);
    set(ax, 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3], 'LineWidth', 1.0);
    set(ax, 'YLim', [0 1.0], 'YTick', 0:0.1:2);
    set(ax, 'TickDir', 'both', 'FontSize', 18);
    xtickangle(ax, 0);
    grid(ax, 'on');
    box(ax, 'off');
    axis(ax, 'square');
    hold(ax, 'off');
    if (noXLabel)
        set(ax, 'XTickLabel', {});
    else
        xlabel(ax, 'spatial frequency (cpd)');
    end
    if (noYLabel)
        set(ax, 'YTickLabel', {});
    else
        ylabel(ax, 'STF (visual space)');
    end

    if (numel(inputConeTypes) <6)
        if (numel(inputConeTypes) == 1)
            coneInfoString = sprintf('input cone: ');
        else
            coneInfoString = sprintf('input cones: ', numel(inputConeTypes));
        end
        for iCone = 1:numel(inputConeTypes)
            switch (inputConeTypes(iCone))
                case cMosaic.LCONE_ID
                    coneInfoString = sprintf('%s L', coneInfoString);
                case cMosaic.MCONE_ID
                    coneInfoString = sprintf('%s M', coneInfoString);
                case cMosaic.SCONE_ID
                    coneInfoString = sprintf('%s S', coneInfoString);
            end
        end
    else
        coneInfoString = sprintf('%d input cones', numel(inputConeTypes));
    end

    coneInfoString = sprintf('RGC %d (%d/%d): %s', theRGCindex, iRGC, totalRGCs, coneInfoString);
    title(ax, sprintf('%s', coneInfoString), 'FontWeight', 'Normal');
end
