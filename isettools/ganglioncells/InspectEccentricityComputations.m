function inspectEccentricityComputations

    % Get dropboxDir & intermediate data files location
    computerInfo = GetComputerInfo();
    switch (computerInfo.localHostName)
        case 'Ithaka'
            dropboxDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris/midgetRGCMosaics';
            mappedRFsDir = '/Volumes/SSDdisk/MATLAB/toolboxes/isetbio/isettools/ganglioncells/JohannesAnalysesDataOLD';
            %mappedRFsDir = '/Volumes/SSDdisk/MATLAB/toolboxes/isetbio/isettools/ganglioncells/JohannesAnalysesData';
            mappedRFsDir = '/Volumes/SSDdisk/MATLAB/toolboxes/isetbio/isettools/ganglioncells/JohannesAnalysesDataBoostedKSParams';
            mappedRFsDir = '/Volumes/SSDdisk/MATLAB/toolboxes/isetbio/isettools/ganglioncells/JohannesAnalysesDataBoost2';
            %mappedRFsDir = '/Volumes/SSDdisk/MATLAB/toolboxes/isetbio/isettools/ganglioncells/JohannesAnalysesDataBoost2NoCompensation';
            mappedRFsDir = '/Volumes/SSDdisk/MATLAB/toolboxes/isetbio/isettools/ganglioncells/JohannesAnalysesDataFixedVnVw';
   
        case 'Crete'
            dropboxDir = '/Volumes/Dropbox/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris/midgetRGCMosaics';
            mappedRFsDir = '/Volumes/MATLAB/toolboxes/isetbio/isettools/ganglioncells/JohannesAnalysesDataBoostedKSPrams';

        otherwise
            if (contains(computerInfo.networkName, 'leviathan'))
                dropboxDir = '/media/dropbox_disk/Aguirre-Brainard Lab Dropbox/isetbio isetbio/midgetRGCMosaics';
            else
                error('Could not establish dropbox location')
            end
    end

    ZernikeDataBase = 'Polans2015';
    subjectRankOrder = 6;
    pupilDiameterMM = 3.0;


    % STF stimulus parameters
    % L+M contrast for gratings used to measure the STFs 
    coneContrasts = [1 1 0];
    
    % 3 degrees in the temporal retina (negative horizontal ecc)
    mosaicEccDegs = [ 6.0 0; ...
                      5.0 0; ...
                      4.0 0; ...
                      3.0 0; ...
                      2.0 0; ...
                      1.0 0; ...
                      0.0 0; ...
                     -0.5 0; ...
                     -1.0 0; ...
                     -2.0 0; ...
                     -3.0 0; ...
                     -4.0 0; ...
                     -5.0 0; ...
                     -6.0 0; ...
                     -8.0 0; ...
                    -10.0 0; ...
                    -12.0 0 ...
                    ];

    if (contains(mappedRFsDir, 'OLD'))
        mosaicEccDegs = [0 0; ...
                     -1 0; ...
                     -2 0; ...
                     -3 0; ...
                     -4.0 0; ...
                     -5 0; ...
                     -6 0; ...
                     -8 0; ...
                     -10 0; ...
                     -12 0 ...
                     ];
        mosaicEccDegs = [6 0];
    else

        mosaicEccDegs = [0 0; ...
                     1 0; ...
                     2 0; ...
                     4.0 0; ...
                     6 0; ...
                     8 0; ...
                     10 0 ...
                     ];

        mosaicEccDegs = [2 0];
    end


%     mosaicEccDegs = [ ...
%         -20 0; ...
%         -16 0; ...
%         -12 0; ...
%         -10 0; ...
%          -8 0; ...
%          -6 0; ...
%          -5 0; ...
%          -4 0; ...
%          -2 0; ...
%          -1 0; ...
%        -0.5 0; ...
%         0.0 0; ...
%         0.5 0; ...
%           1 0; ...
%           2 0; ...
%           3 0; ...
%           4 0; ...
%           5 0; ...
%           6 0; ...
%           8 0; ...
%          10 0; ...
%          12 0; ...
%          20 0];



    % Αctions
    centerMostRGCsNumToAnalyze = 100;
    inspectTheSTFs = true;
    inspectTheSpatialRFs = true;
    contrastModelToCronerAndKaplan = true;
    
    if (inspectTheSTFs)
        for iEcc = 1:size(mosaicEccDegs,1)
            inspectSTFs(mosaicEccDegs(iEcc,:), mappedRFsDir, ...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                coneContrasts, centerMostRGCsNumToAnalyze);
        end
    end

   
    if (inspectTheSpatialRFs)
        for iEcc = 1:size(mosaicEccDegs,1)
            inspectSpatialRFs(mosaicEccDegs(iEcc,:), mappedRFsDir, ...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                coneContrasts, centerMostRGCsNumToAnalyze);
        end
    end
    

    
    if (contrastModelToCronerAndKaplan)
        for iEcc = 1:size(mosaicEccDegs,1)
            dataOut{iEcc} = extractPatchDoGParams(mosaicEccDegs(iEcc,:), mappedRFsDir, ...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                coneContrasts, centerMostRGCsNumToAnalyze);
        end
        contastDerivedDoGparamsToCronerAndKaplanDoGParams(dataOut);
    end
end


function contastDerivedDoGparamsToCronerAndKaplanDoGParams(dataOut)

    mosaicTemporalEccDegs = [];
    modelRGCRcDegs = [];
    modelRGCtemporalEccDegs = [];
    modelRGCRsRcRatios = [];
    modelRGCKsKcRatios = [];
    modelRGCSCintSensRatios = [];
    targetRsRcRatios = [];
    targetSCintSensRatios = [];

    for iEcc = 1:numel(dataOut)
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
    targetRsRcRatios  = targetRsRcRatios(idx);
    targetSCintSensRatios = targetSCintSensRatios(idx);


    hFig = figure(4); clf;
    set(hFig, 'Position', [10 10 1650 600], 'Color', [1 1 1]);

    plotKcKs = true;
    ax = subplot('Position', [0.05 0.13 0.25 0.8]);  
    if (plotKcKs)
        [CronerKaplanTemporalEccDegs, CronerKaplanKsKcDegs] = ...
            RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterPeakSensisitivityRatioAgainstEccentricity();

        scatter(modelRGCtemporalEccDegs, modelRGCKsKcRatios, 100, ...
             'LineWidth', 1.0, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor',[1 0.5 0.5]);
        hold(ax, 'on');

        scatter(CronerKaplanTemporalEccDegs, CronerKaplanKsKcDegs, 144, ...
             'LineWidth', 1.0, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[0.7 0.7 0.7]);
    
        legend({'ISETBio midget RGCs', 'macaque midget RGCs (Croner & Kaplan)'}, 'Location', 'NorthOutside', 'box', 'off', 'FontSize', 15);
       
        grid 'on'
        set(gca, 'XLim', [0.1 30], 'XScale', 'log', 'XTick', [0.1 0.3 1 3 10 30], ...
             'YLim', [1e-5 1], 'YScale', 'log', 'YTick', [1e-5 1e-4 1e-3 1e-2 1e-1 1], 'FontSize', 20);
        set(gca, 'TickDir', 'both');
        xlabel('temporal equivalent eccentericity (degs)');
        ylabel('Ks/Kc');

    else

        [CronerKaplanTemporalEccDegs, CronerKaplanRcDegs] = ...
            RGCmodels.CronerKaplan.digitizedData.parvoCenterRadiusAgainstEccentricity();
    
        scatter(modelRGCtemporalEccDegs, modelRGCRcDegs*60, 100, ...
             'LineWidth', 1.0, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor',[1 0.5 0.5]);
        hold(ax, 'on');
    
        scatter(CronerKaplanTemporalEccDegs, CronerKaplanRcDegs*60, 144, ...
             'LineWidth', 1.0, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[0.7 0.7 0.7]);
    
        legend({'ISETBio midget RGCs', 'macaque midget RGCs (Croner & Kaplan)'}, 'Location', 'NorthOutside', 'box', 'off', 'FontSize', 15);
        
        grid 'on'
        set(gca, 'XLim', [0.1 30], 'XScale', 'log', 'XTick', [0.1 0.3 1 3 10 30], ...
             'YLim', [0 14], 'YTick', 0:2:20, 'FontSize', 20);
        set(gca, 'TickDir', 'both');
        xlabel('temporal equivalent eccentericity (degs)');
        ylabel('Rc (arc min)');
    end



    ax = subplot('Position', [0.39 0.13 0.25 0.8]); 
    % Retrieve Croner&Kaplan Rs/Rc ratios
    [CronerKaplanTemporalEccDegs, CronerKaplanRcRsRatios] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();

    scatter(modelRGCtemporalEccDegs, modelRGCRsRcRatios, 100, ...
         'LineWidth', 1.0, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor',[1 0.5 0.5]);
    hold(ax, 'on');

    scatter(CronerKaplanTemporalEccDegs, 1./CronerKaplanRcRsRatios, 144, ...
         'LineWidth', 1.0, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[0.7 0.7 0.7]);

    %plot(mosaicTemporalEccDegs, targetRsRcRatios, 'b-', 'LineWidth',2);
    legend({'ISETBio midget RGCs', 'macaque midget RGCs (Croner & Kaplan)'}, 'NumColumns', 1, ...
        'Location', 'NorthOutside', 'box', 'off', 'FontSize', 15);
    
    grid 'on'
    set(gca, 'XLim', [0.1 30], 'XScale', 'log', 'XTick', [0.1 0.3 1 3 10 30], ...
         'YLim', [0 24], 'YTick', 0:2:25, 'FontSize', 20);
    set(gca, 'TickDir', 'both');
    xlabel('temporal equivalent eccentericity (degs)');
    ylabel('Rs/Rc ratio');


    
    ax = subplot('Position', [0.73 0.13 0.25 0.8]);
    % Retrieve Croner&Kaplan S/C int. sensititivity data
    [CronerKaplanTemporalEccDegs, CronerKaplanSCintSensRatios] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();

    scatter(modelRGCtemporalEccDegs, modelRGCSCintSensRatios, 100, ...
         'LineWidth', 1.0, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor',[1 0.5 0.5]);
    hold(ax, 'on');

    scatter(CronerKaplanTemporalEccDegs, CronerKaplanSCintSensRatios, 144, ...
         'LineWidth', 1.0, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[0.7 0.7 0.7]);

    %plot(mosaicTemporalEccDegs, targetSCintSensRatios , 'b-', 'LineWidth', 2);
    legend({'ISETBio midget RGCs', 'macaque midget RGCs (Croner & Kaplan)'}, ...
        'NumColumns', 1, 'Location', 'NorthOutside', 'box', 'off', 'FontSize', 15);
    
    grid 'on'
    set(gca, 'XLim', [0.1 30], 'XScale', 'log', 'XTick', [0.1 0.3 1 3 10 30], ...
         'YLim', [0 1.0], 'YTick', 0:0.1:1.2, 'FontSize', 20);
    set(gca, 'TickDir', 'both');
    xlabel('temporal equivalent eccentericity (degs)');
    ylabel('S/C int. sensitivity ratio');

    
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
        numel(fittedSTFs)
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
    pupilDiameterMM, coneContrasts, centerMostRGCsNumToAnalyze)

    % Load the midgetRGCMosaic
    fName = fullfile(mappedRFsDir, sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2)));
    load(fName, 'theMidgetRGCmosaic');

    videoFileName = strrep(fName, 'mat', '_SpatialRFs.mp4');

    maxVisualizedRFs = 100;
    coVisualizeSTFs = true;

    if (coVisualizeSTFs == false)
    
        theMidgetRGCmosaic.visualizeSpatialRFs(...
                'maxVisualizedRFs', maxVisualizedRFs, ...
                'generateVideo', true);
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
                'generateVideo', false);
    
            % replace graphic is (1,1) with the STFs
            noXLabel = false; noYLabel = false;
            cla(allAxes{1,1});
            visualizeSTFs(allAxes{1,1}, spatialFrequenciesTested, ...
                        f.allMeasuredSTFs, f.theMeasuredSTFtoFit, ...
                        theFittedSTF.sfHiRes, theFittedSTF.compositeSTFHiRes, ...
                        theTargetSTFsupport, theTargetSTFmeasured, ...
                        noXLabel, noYLabel, inputConeTypes, theRGCindex);
    
            drawnow;
            videoOBJ.writeVideo(getframe(hFig));
    
        end

        videoOBJ.close();
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
                noXLabel, noYLabel, inputConeTypes, theRGCindex);

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


end

function shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
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
    noXLabel, noYLabel, inputConeTypes, theRGCindex)

    
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
        ylabel(ax, 'STF');
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

    coneInfoString = sprintf('RGC %d: %s', theRGCindex, coneInfoString);
    title(ax, sprintf('%s', coneInfoString), 'FontWeight', 'Normal');
end
