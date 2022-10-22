function generateMidgetRGCmosaicComponents
    
    % Operations
    operations = {...
        'generateMosaic' ...
        'generateRTVFTobjList' ...
        'generateCenterSurroundRFs' ...
        'computeSTF' ...
        'visualizeSTFs' ...
        'fitSTFs' ...
        'summarizeSTFfits'
        };

    % Operation to compute
    %operations = {'fitSTFs'};
    %operations = {'summarizeSTFfits'};
    operations = operations(1:2);

    % L-M gratings
    %coneContrasts = [0.12 -0.12 0];

    % L+M gratings
    coneContrasts = [1 1 0];

    % Temporal retina (negative x-coords)
    eccSizeDegsExamined = [...
         0  0.5; ...
        -1  0.5; ...
        -2  0.5; ...
        -3  0.6; ...
        -4  0.8; ...
        -6  1.0; ...
        -8  1.2; ...
        -10 1.4; ...
        -12 1.6; ...
        -14 1.8; ...
        -16 2.0; ...
        -20 2.2; ...
        -25 2.5 ...
       ];
    
    % Get dropboxDir location
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


    resetSummaryFigure  = true; hFigSummary = [];
   
    for iEcc = 8:8 %size(eccSizeDegsExamined,1)
        
        fprintf('Generating components for mosaic %d of %d\n', iEcc, size(eccSizeDegsExamined,1));
        eccDegs  = eccSizeDegsExamined(iEcc,1) * [1 0];
        sizeDegs = eccSizeDegsExamined(iEcc,2) * [1 1];
        hFigSummary = doIt(operations, eccDegs, sizeDegs, coneContrasts, dropboxDir, resetSummaryFigure, hFigSummary);
        resetSummaryFigure = false;
    end

    if (strcmp(operations{1}, 'summarizeSTFfits'))
        NicePlot.exportFigToPDF('ModelComparisonToCronerKaplan.pdf', hFigSummary, 300);
    end


end

function hFigSummary = doIt(operations, eccDegs, sizeDegs, coneContrasts, dropboxDir, resetSummaryFigure, hFigSummary)

    fName = fullfile(dropboxDir, sprintf('mRGCmosaicComponents_eccDegs_%2.2f.mat', eccDegs(1)));

    for iOp = 1:numel(operations)

        switch (operations{iOp})
            case 'generateMosaic'
                fprintf('Generating midgetRGCmosaic ...\n');
                % Generate the midgetRGCmosaic
                theMidgetRGCmosaic = midgetRGCMosaic(...
                        'sourceLatticeSizeDegs', 60, ...
                        'eccentricityDegs', eccDegs, ...
                        'sizeDegs', sizeDegs ...
                        );
                save(fName, 'theMidgetRGCmosaic', '-v7.3');
                fprintf('Exported the computed midgetRGCMosaic to %s\n', fName);

            case 'generateRTVFTobjList'
                fprintf('Generating RTVFTobjList ... \n');

                % Load the midget RGCmosaic
                load(fName, 'theMidgetRGCmosaic');

                % Find out the range of cones in the RF center
                conesNumPooledByTheRFcenter = unique(full(sum(theMidgetRGCmosaic.rgcRFcenterConeConnectivityMatrix,1)));
                fprintf('Cones/RF center for this mosaic: %d\n', conesNumPooledByTheRFcenter);

                % Generate the grids for different
                % opticsPositions & conesNumInRFcenters
                for iGridPosition = 1:numel(conesNumPooledByTheRFcenter)
                    % Optics position
                    eccDegsGrid(iGridPosition,:) = eccDegs;
                    
                    % Cones num in RF center
                    conesNumPooledByTheRFcenterGrid(iGridPosition) = conesNumPooledByTheRFcenter(iGridPosition);

                    % From Croner & Kaplan '95 (Figure 4c and text)
                    % "P surrounds were on average 6.7 times wider than the centers of
                    % the same cells, or about 45 times larger in area".
                    surroundToCenterRcRatioGrid(iGridPosition) = RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio;
    
                    % Also from Croner & Kaplan '95 (Figure 10b)
                    % "These mean ratios for P and M cells are not significantly different
                    % (Student's t-test: P = 0.482). The overall mean ratio is 0.55.
                    temporalEquivalentEccDegs = theMidgetRGCmosaic.temporalEquivalentEccentricityForEccentricity(eccDegsGrid(iGridPosition,:));
                    radialTemporalEquivalentEccDegs = sqrt(sum(temporalEquivalentEccDegs.^2,2));
                    scIntSensitivity = RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(radialTemporalEquivalentEccDegs);
                    surroundToCenterIntegratedSensitivityRatioGrid(iGridPosition) = scIntSensitivity;
                    
                end


                % Compute a list of RTVFTobj for each of the examined grid positions
                RTVFTobjList = generateRTVFTobjects(theMidgetRGCmosaic, ...
                    eccDegsGrid, conesNumPooledByTheRFcenterGrid, ...
                    surroundToCenterRcRatioGrid, surroundToCenterIntegratedSensitivityRatioGrid);

                % Save the computed list of RTVFTobj for each of the examined grid positions
                save(fName, ...
                    'RTVFTobjList', ...
                    'eccDegsGrid', ...
                    'conesNumPooledByTheRFcenterGrid', ...
                    'surroundToCenterRcRatioGrid', ...
                    'surroundToCenterIntegratedSensitivityRatioGrid', ...
                    '-append');
                fprintf('Appended the computed RTVFTobjList to %s\n', fName);

            case 'generateCenterSurroundRFs'
                load(fName, ...
                    'theMidgetRGCmosaic', ...
                    'RTVFTobjList', ...
                    'eccDegsGrid', ...
                    'conesNumPooledByTheRFcenterGrid', ...
                    'surroundToCenterRcRatioGrid', ...
                    'surroundToCenterIntegratedSensitivityRatioGrid');

                    % Generate C/S spatial RFs for all cells in the
                    % midgetRGCmosaic
                    theMidgetRGCmosaic.generateCenterSurroundSpatialPoolingRF(RTVFTobjList, ...
                        eccDegsGrid, conesNumPooledByTheRFcenterGrid, ...
                        surroundToCenterRcRatioGrid, surroundToCenterIntegratedSensitivityRatioGrid);

                    clear 'RTVFTobjList';

                    % Save the updated midgetRGCmosaic which now includes
                    % the computed RTVFTobjList
                    save(fName, ...
                        'theMidgetRGCmosaic', ...
                        'eccDegsGrid', ...
                        'conesNumPooledByTheRFcenterGrid', ...
                        'surroundToCenterRcRatioGrid', ...
                        'surroundToCenterIntegratedSensitivityRatioGrid', ...
                        '-v7.3');

            case 'computeSTF'
                load(fName, 'theMidgetRGCmosaic');
                
                [theMidgetRGCMosaicResponses, spatialFrequenciesTested, spatialPhasesDegs] = ...
                    computeTheSTF(theMidgetRGCmosaic, coneContrasts);
                
                % Save the responses to a separate file
                responsesPostfix = sprintf('_Responses_%2.2f_%2.2f_%2.2f.mat', ...
                    coneContrasts(1), coneContrasts(2), coneContrasts(3));
                fNameResponses = strrep(fName, '.mat', responsesPostfix);
                save(fNameResponses, ...
                    'theMidgetRGCMosaicResponses', ...
                    'spatialFrequenciesTested', ...
                    'spatialPhasesDegs', ...
                    'coneContrasts', ...
                    '-v7.3');

           
            case {'visualizeSTFs','fitSTFs'} 
                load(fName, 'theMidgetRGCmosaic');
                
                % Load the responses to a separate file
                responsesPostfix = sprintf('_Responses_%2.2f_%2.2f_%2.2f.mat', ...
                    coneContrasts(1), coneContrasts(2), coneContrasts(3));
                fNameResponses = strrep(fName, '.mat', responsesPostfix);
                load(fNameResponses, ...
                    'theMidgetRGCMosaicResponses', ...
                    'spatialFrequenciesTested', ...
                    'spatialPhasesDegs', ...
                    'coneContrasts');

                % Sort RGCs according to their eccentricity
                mRGCmosaicCenterDegs = mean(theMidgetRGCmosaic.rgcRFpositionsDegs,1);
   
                ecc = sum((bsxfun(@minus, theMidgetRGCmosaic.rgcRFpositionsDegs, mRGCmosaicCenterDegs)).^2,2);
                [~,sortedRGCindices] = sort(ecc, 'ascend');

                
                hFig = figure(66); clf;
                set(hFig, 'Position', [90 10 855 990], 'Color', [1 1 1]);

                % Video setup
                if (strcmp(operations{iOp}, 'visualizeSTFs'))
                    fNameVideo = strrep(fName, '.mat', '_Video');
                    videoOBJ = VideoWriter(fNameVideo, 'MPEG-4');
                    videoOBJ.FrameRate = 10;
                    videoOBJ.Quality = 100;
                    videoOBJ.open();
                end

                skips = 1;
                if (strcmp(operations{iOp}, 'fitSTFs'))
                    
                    maxFitsNum = 100;
                    skips = round(numel(sortedRGCindices)/maxFitsNum);
                    if (skips < 1)
                        skips = 1;
                    end
                    DoGparams = cell(1, numel(1:skips:numel(sortedRGCindices)));
                    DoGfitResults = cell(1,numel(1:skips:numel(sortedRGCindices)));

                    
                end

                iFit = 0;
                for iii = 1:skips:numel(sortedRGCindices)
                    iRGC = sortedRGCindices(iii);
                    
                    connectivityVector = full(squeeze(theMidgetRGCmosaic.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
                    indicesOfCenterCones = find(abs(connectivityVector) > 0.0001);
                    
                    % Retrieve the correct RTVFTobj based on this cells position and
                    % #of center cones. For now only checking the centerConesNum
                    iObj = find(...
                        (theMidgetRGCmosaic.theConesNumPooledByTheRFcenterGrid == numel(indicesOfCenterCones)) ...  % match the conesNum in the center
                    );
                    theRTVFTobj = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{iObj};

                    maxResponse = max(max(abs(squeeze(theMidgetRGCMosaicResponses(:, :, iRGC)))));
                    minResponse = -maxResponse;
                    meanResponse = 0;

                    theResponseModulation = zeros(1, numel(spatialFrequenciesTested));
                    for iSF = 1:numel(spatialFrequenciesTested)
                        % Retrieve the mRGC response time-series

                        theResponse = squeeze(theMidgetRGCMosaicResponses(iSF, :, iRGC));
                        % Compute the response modulation for this SF
                        theResponseModulation(iSF) = max(theResponse)-min(theResponse);
                        
                        % Plot the time-series response for this SF
                        ax = subplot(numel(spatialFrequenciesTested),2,(iSF-1)*2+1);
                        plot(ax, 1:numel(spatialPhasesDegs), 0*theResponse, 'k-', 'LineWidth', 1.0);
                        hold(ax, 'on');
                        plot(ax, 1:numel(spatialPhasesDegs), theResponse, 'bo-', 'MarkerSize', 10, 'MarkerFaceColor', [0.3 0.8 0.8], 'LineWidth', 1.0);
                        hold(ax, 'off');
                        set(ax,  'YLim', [minResponse maxResponse], 'XTick', [], ...
                            'YTick', [minResponse meanResponse maxResponse], ...
                            'YTickLabel', sprintf('%2.2f\n',[minResponse meanResponse maxResponse]), ...
                            'XColor', 'none');

                        title(ax, sprintf('%2.1f c/deg', spatialFrequenciesTested(iSF)));
                        box(ax, 'off');
                        if (iSF == numel(spatialFrequenciesTested))
                            xlabel(ax, 'time');
                        end
                        ylabel(ax, 'response');
                    end
                    theMeasuredSTF = theResponseModulation/max(theResponseModulation);

                    % Visualize the examined retinal RF
                    ax = subplot(numel(spatialFrequenciesTested),2, [2 4 6 8]);
                    theMidgetRGCmosaic.visualizeSingleRetinalRF(iRGC, ...
                        'plotTitle', sprintf('RGC: %d of %d', iii, numel(sortedRGCindices)), ...
                        'figureHandle', hFig, ...
                        'axesHandle', ax);

                    % Visualize the computed STF
                    ax = subplot(numel(spatialFrequenciesTested),2, ((6:numel(spatialFrequenciesTested))-1)*2+2);
                    
                    % The target STF
                    normalizedTargetSTF = theRTVFTobj.rfComputeStruct.theSTF.target;
                    normVal = max(normalizedTargetSTF);
                    normalizedTargetSTF1 = normalizedTargetSTF / normVal;

                    % The STF achieved by the RTVT - this is what the cone pooling weights are computed from
                    normalizedTargetSTF = theRTVFTobj.rfComputeStruct.theSTF.fitted;
                    normVal = max(normalizedTargetSTF);
                    normalizedTargetSTF2 = normalizedTargetSTF / normVal;
                    

                    % Plot the target STF
                    p1 = plot(ax, theRTVFTobj.rfComputeStruct.theSTF.support, normalizedTargetSTF1, 'r-', 'Color', [1 0.5 0.5], 'LineWidth', 2.0);
                    hold(ax, 'on')
                    p2 = plot(ax, theRTVFTobj.rfComputeStruct.theSTF.support, normalizedTargetSTF2, 'r--', 'LineWidth', 2.0);
                    
                    % Plot the measured STF
                    p3 = plot(ax,spatialFrequenciesTested, theMeasuredSTF, 'ko', ...
                        'MarkerSize', 14, 'MarkerFaceColor', [0.2 0.9 0.9], 'MarkerEdgeColor', [0 0.4 1], ...
                        'LineWidth', 1.0);

                    % Fit the measured STF to extract the DoGparams
                    if (strcmp(operations{iOp}, 'fitSTFs'))
                        % Estimate the retinal RF center
                        conesNumPooledByTheRFcenter = numel(indicesOfCenterCones);
                        coneRcDegs = mean(theMidgetRGCmosaic.inputConeMosaic.coneApertureDiametersDegs(indicesOfCenterCones)) * ...
                                          theMidgetRGCmosaic.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;
                        retinalRFcenterRcDegs = sqrt(conesNumPooledByTheRFcenter)*coneRcDegs;

                        iFit = iFit + 1;
                        % Fit the DoG model to the measured STF
                        [DoGparams{iFit}, theFittedSTF] = fitDoGmodelToMeasuredSTF(spatialFrequenciesTested, theMeasuredSTF, retinalRFcenterRcDegs);
                        
                        % Save the fit results
                        DoGfitResults{iFit} = struct(...
                             'targetRGC', iRGC, ...
                             'targetRGCeccentricityDegs', theMidgetRGCmosaic.rgcRFpositionsDegs(iRGC,:), ...
                             'targetVisualRFDoGparams', theRTVFTobj.targetVisualRFDoGparams, ...
                             'spatialFrequenciesTested', spatialFrequenciesTested, ...
                             'theMeasuredSTF',theMeasuredSTF, ...
                             'theFittedSTF', theFittedSTF ...
                             );
                        
                        % Plot the fitted STF
                        plot(ax, theFittedSTF.sfHiRes, theFittedSTF.compositeSTFHiRes, 'k-', 'LineWidth', 3.0);
                        p4 = plot(ax, theFittedSTF.sfHiRes, theFittedSTF.compositeSTFHiRes, 'c-', 'LineWidth', 1.5);
                    end

                    hold(ax, 'off');
                    if (strcmp(operations{iOp}, 'fitSTFs'))
                        legend([p1, p2, p3, p4], {'target', 'fitted', 'measured', 'fitToMeasured'}, 'Location', 'SouthWest');
                    else
                        legend([p1, p2, p3], {'target', 'fitted', 'measured'}, 'Location', 'SouthWest');
                    end

                    title(ax, sprintf('stim LMS contrast = < %2.1f, %2.1f, %2.1f >', coneContrasts(1), coneContrasts(2), coneContrasts(3)));
                    set(ax, 'XLim', [0.3 70], 'XTick', [0.1 0.3 1 3 10 30 100], 'YLim', [0 1.02], 'YTick', 0:0.1:1.0);
                    xlabel(ax, 'spatial frequency (c/deg)');
                    ylabel(ax, 'STF');
                    grid(ax, 'on');
                    set(ax, 'XLim', [0.1 100], 'XScale', 'Log', 'FontSize', 16);
                    
                    if (strcmp(operations{iOp}, 'visualizeSTFs'))
                        drawnow;
                        videoOBJ.writeVideo(getframe(hFig));
                    end
                end

                if (strcmp(operations{iOp}, 'visualizeSTFs'))
                    videoOBJ.close();
                end

                if (strcmp(operations{iOp}, 'fitSTFs'))
                    % Save the DoGparams and the fits to a separate file
                    dogParamsPostfix = sprintf('_fittedDoGmodels_%2.2f_%2.2f_%2.2f.mat', ...
                        coneContrasts(1), coneContrasts(2), coneContrasts(3));
                    fNameDoGparams = strrep(fName, '.mat', dogParamsPostfix );
                    save(fNameDoGparams, 'DoGparams', 'DoGfitResults');
                end

            case 'summarizeSTFfits'
                % Load the DoGparams and the fits to a separate file
                dogParamsPostfix = sprintf('_fittedDoGmodels_%2.2f_%2.2f_%2.2f.mat', ...
                        coneContrasts(1), coneContrasts(2), coneContrasts(3));
                fNameDoGparams = strrep(fName, '.mat', dogParamsPostfix );
                load(fNameDoGparams, 'DoGparams', 'DoGfitResults');

                if (resetSummaryFigure)
                    hFigSummary = figure(1);
                    clf;
                    set(hFigSummary, 'Position', [10 10 1100 1100], 'Color', [1 1 1]);
                end

                addModelData = true;
                hFigSummary = visualizeSTFfitsSummary(DoGparams, DoGfitResults, addModelData, resetSummaryFigure, hFigSummary);

            otherwise
                error('Unknown operation: ''%s''.', operations{iOp});
        end % Switch
    end


end


function hFig = visualizeSTFfitsSummary(DoGparams, DoGfitResults, addModelData, reset, hFig)

    ksKcRatio = zeros(numel(DoGparams),1);
    rcDegs = zeros(numel(DoGparams),1);
    rsDegs = zeros(numel(DoGparams),1);
    targetRGCeccentricityDegs = zeros(numel(DoGparams),2);
    for iii = 1:numel(DoGparams)
        rcDegs(iii) = DoGparams{iii}.bestFitValues(4);
        rsDegs(iii) = rcDegs(iii) * DoGparams{iii}.bestFitValues(3);
        ksKcRatio(iii) = DoGparams{iii}.bestFitValues(2);
        targetSurroundToCenterRcRatio = DoGfitResults{iii}.targetVisualRFDoGparams.surroundToCenterRcRatio;
        targetSurroundToCenterIntegratedSensitivityRatio = DoGfitResults{iii}.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;
        targetRGCeccentricityDegs(iii,:) = DoGfitResults{iii}.targetRGCeccentricityDegs;
    end

    targetRGCeccentricityDegs = sqrt(sum(targetRGCeccentricityDegs.^2,2));
    rsRcRatio = rsDegs./rcDegs;
    integratedSCsensitivity = ksKcRatio .* (rsRcRatio.^2);

    targetKsKcRatio = targetSurroundToCenterIntegratedSensitivityRatio / ((targetSurroundToCenterRcRatio)^2);


    eccRange = [0 30];
    eccTicks = 0:5:30;

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 2, ...
       'heightMargin',  0.10, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.04, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.06, ...
       'topMargin',      0.02);

    % Rc, Rs as a function of eccentricity
    ax = subplot('Position', subplotPosVectors(1,1).v);

    % The C&K data
    [eccDegs, RcDegsCronerKaplan] = RGCmodels.CronerKaplan.digitizedData.parvoCenterRadiusAgainstEccentricity();
    p1 = plot(ax,eccDegs,   RcDegsCronerKaplan,   'ro',  ...
        'MarkerFaceColor', [0.65 0.65 0.65], 'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'MarkerSize', 8, 'LineWidth', 1.0);
    hold(ax, 'on');
    [eccDegs, RsDegsCronerKaplan] = RGCmodels.CronerKaplan.digitizedData.parvoSurroundRadiusAgainstEccentricity();
    p2 = plot(ax,eccDegs, RsDegsCronerKaplan, 'rs', ...
        'MarkerFaceColor', [0.85 0.85 0.85], 'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'MarkerSize', 10, 'LineWidth', 1.0);
    
    % The model data
    if (addModelData)
        p3 = scatter(ax,targetRGCeccentricityDegs, rcDegs, 81, ...
            'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.65, ...
            'MarkerFaceColor', [1 0.5 0.7],  'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 0.75);
        p4 = scatter(ax,targetRGCeccentricityDegs, rsDegs, 81, 'Marker', 's', ...
            'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.65, ...
            'MarkerFaceColor', [0.3 0.8 1],  'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 0.75);
    end

    if (addModelData)
        legend(ax,[p1,p2,p3,p4], {'C&K centers', 'C&K surrounds', '@mRGCMosaic centers', '@mRGCMosaic surrounds'}, ...
            'NumColumns', 2, 'Location', 'NorthOutside', 'box', 'off');
    else
        legend(ax,[p1,p2], {'C&K centers', 'C&K surrounds'}, ...
            'NumColumns', 1, 'Location', 'NorthOutside', 'box', 'off');
    end

    % Finalize plot
    axis(ax, 'square')
    grid(ax, 'on');
    box(ax, 'off');
    set(ax, 'XLim', [0.03 eccRange(2)], 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30], ...
            'YLim', [0.003 1], 'YTick', [0.001 0.003 0.01 0.03 0.1 0.3 1 3 10], ...
            'XScale', 'log', 'YScale', 'log', ...
            'FontSize', 20, 'LineWidth', 1.0, 'TickDir', 'both');
    
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'radius (degs)');


    % Integrated S/C ratio as a function of eccentricity
     ax = subplot('Position', subplotPosVectors(1,2).v);
    p1 = plot(ax,targetRGCeccentricityDegs, integratedSCsensitivity, 'r.', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.5]);
    hold(ax, 'on');
    box(ax, 'off');
    % The C&K data
    %eccCont = linspace(eccRange(1), eccRange(2), 100);
    %scIntSensCont = 0.466 + 0.007*eccCont; % C&K figure 11
    %p2 = plot(eccRange, targetSurroundToCenterIntegratedSensitivityRatio*[1 1], 'r-', 'LineWidth', 1.0);
    %p2 = plot(eccCont, scIntSensCont, 'r-', 'LineWidth', 1.0);
    
    [ecc, intSCratioCronerKaplan] = RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();
    p2 = plot(ax,ecc, intSCratioCronerKaplan, 'rs',  'MarkerEdgeColor', [0.5 0. 0], 'MarkerFaceColor', [0.85 0.85 0.85], 'MarkerSize', 10, 'LineWidth', 1.0);
    plot(ax,targetRGCeccentricityDegs, integratedSCsensitivity, 'r.', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.5]);
    legend(ax,[p1,p2], {'@mRGCMosaic', 'C&K'}, 'NumColumns', 2, 'Location', 'NorthOutside', 'box', 'off');
    
    axis(ax, 'square')
    grid(ax, 'on');
    box(ax, 'off');
    set(ax, 'YLim', [0 1], 'XLim', eccRange, 'XTick', eccTicks, 'FontSize', 20, 'LineWidth', 1.0, 'TickDir', 'both');
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'integrated surround/center sensitivity');


    % Center/Surround radius as a function of eccentricity
     ax = subplot('Position', subplotPosVectors(2,1).v);
    p1=plot(ax,targetRGCeccentricityDegs, 1./rsRcRatio, 'r.', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.5]);
    hold(ax, 'on');
    % The C&K data
    %plot(eccRange, 1./targetSurroundToCenterRcRatio*[1 1], 'r-', 'LineWidth', 1.0);
    [ecc, RcRsRatioCronerKaplan] = RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();
    p2 = plot(ax,ecc, RcRsRatioCronerKaplan, 'rs',  'MarkerEdgeColor', [0.5 0. 0], 'MarkerFaceColor', [0.85 0.85 0.85], 'MarkerSize', 10, 'LineWidth', 1.0);
    plot(ax,targetRGCeccentricityDegs, 1./rsRcRatio, 'r.', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.5]);

    legend(ax,[p1,p2], {'@mRGCMosaic', 'C&K'}, 'NumColumns', 2, 'Location', 'NorthOutside', 'box', 'off');
    axis(ax, 'square')
    grid(ax, 'on');
    box(ax, 'off');
    set(ax, 'YLim', [0 0.5], 'YTick', 0:0.1:1.0, 'XLim', eccRange, 'XTick', eccTicks, 'FontSize', 20, 'LineWidth', 1.0, 'TickDir', 'both');
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'center/surround radius');

    % Ks/Kc as a function of eccentricity
     ax = subplot('Position', subplotPosVectors(2,2).v);
    p1 = plot(ax,targetRGCeccentricityDegs, ksKcRatio, 'r.', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.5]);
    hold(ax, 'on');
    % The C&K data
    [ecc, ksKcRatioCronerKaplan] = RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterPeakSensisitivityRatioAgainstEccentricity();
    p2 = plot(ax,ecc, ksKcRatioCronerKaplan, 'ks',  'MarkerEdgeColor', [0.5 0. 0], 'MarkerFaceColor', [0.85 0.85 0.85], 'MarkerSize', 10, 'LineWidth', 1.0);
    plot(ax,targetRGCeccentricityDegs, ksKcRatio, 'r.', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.5]);
    %plot(eccRange, targetKsKcRatio*[1 1], 'r-', 'LineWidth', 1.0);
    legend(ax,[p1,p2], {'@mRGCMosaic', 'C&K'}, 'NumColumns', 2, 'Location', 'NorthOutside', 'box', 'off');
    set(ax, 'YLim', [1e-4 1], 'YScale', 'log', 'YTick', [1e-4 1e-3 1e-2 1e-1 1], 'YTickLabel',{'1e-4', '1e-3', '1e-2', '1e-1', '1'}, ...
        'XLim', eccRange, 'XTick', eccTicks, 'FontSize', 20,  'LineWidth', 1.0, 'TickDir', 'both');
    axis(ax, 'square')
    grid(ax, 'on'); box(ax, 'off');
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'kS/Kc ratio');

end



function [DoGparams, theFittedSTF] = fitDoGmodelToMeasuredSTF(sf, theMeasuredSTF, retinalRFcenterRcDegs)
    % DoG param initial values and limits: center gain, kc
    Kc = struct(...    
        'low', 1e-4, ...
        'high', 1e5, ...
        'initial', 1);

    % DoG param initial values and limits: Ks/Kc ratio
    KsToKc = struct(...
        'low', 1e-6, ...
        'high', 1, ...
        'initial', 0.1);

    % DoG param initial values and limits: RsToRc ratio
    RsToRc = struct(...
        'low', 1.5, ...
        'high', 10, ...
        'initial', 5);

    % DoG param initial values and limits: RcDegs
    RcDegs = struct(...
        'low', retinalRFcenterRcDegs/10, ...
        'high', retinalRFcenterRcDegs*200, ...
        'initial', retinalRFcenterRcDegs*5);
    
     %                          Kc           kS/kC             RsToRc            RcDegs    
     DoGparams.initialValues = [Kc.initial   KsToKc.initial    RsToRc.initial    RcDegs.initial];
     DoGparams.lowerBounds   = [Kc.low       KsToKc.low        RsToRc.low        RcDegs.low];
     DoGparams.upperBounds   = [Kc.high      KsToKc.high       RsToRc.high       RcDegs.high];
     DoGparams.names         = {'Kc',        'kS/kC',         'RsToRc',         'RcDegs'};
     DoGparams.scale         = {'log',       'log',           'linear',         'linear'};
     
     % The DoG model in the frequency domain
     DoGSTF = @(params,sf)(...
                    abs(params(1)       * ( pi * params(4)^2             * exp(-(pi*params(4)*sf).^2) ) - ...
                    params(1)*params(2) * ( pi * (params(4)*params(3))^2 * exp(-(pi*params(4)*params(3)*sf).^2) )));
        
     % The optimization objective
     objective = @(p) sum((DoGSTF(p, sf) - theMeasuredSTF).^2);

     % Ready to fit
     options = optimset(...
            'Display', 'off', ...
            'Algorithm', 'interior-point',... % 'sqp', ... % 'interior-point',...
            'GradObj', 'off', ...
            'DerivativeCheck', 'off', ...
            'MaxFunEvals', 10^5, ...
            'MaxIter', 10^3);
        
     % Multi-start
     problem = createOptimProblem('fmincon',...
          'objective', objective, ...
          'x0', DoGparams.initialValues, ...
          'lb', DoGparams.lowerBounds, ...
          'ub', DoGparams.upperBounds, ...
          'options', options...
          );
      
     ms = MultiStart(...
          'Display', 'off', ...
          'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
          'UseParallel', true);
      
     % Run the multi-start
     multiStartsNum = 16;
     DoGparams.bestFitValues = run(ms, problem, multiStartsNum);

     theFittedSTF.compositeSTF = DoGSTF(DoGparams.bestFitValues, sf);
     theFittedSTF.centerSTF = DoGparams.bestFitValues(1) * ( pi * DoGparams.bestFitValues(4)^2 * exp(-(pi*DoGparams.bestFitValues(4)*sf).^2) );
     theFittedSTF.surroundSTF = DoGparams.bestFitValues(1)*DoGparams.bestFitValues(2) * ( pi * (DoGparams.bestFitValues(4)*DoGparams.bestFitValues(3))^2 * exp(-(pi*DoGparams.bestFitValues(4)*DoGparams.bestFitValues(3)*sf).^2) );
     
     sfHiRes = logspace(log10(0.1), log10(100), 64);
     theFittedSTF.sfHiRes = sfHiRes;
     theFittedSTF.compositeSTFHiRes = DoGSTF(DoGparams.bestFitValues, sfHiRes);
     theFittedSTF.centerSTFHiRes = DoGparams.bestFitValues(1) * ( pi * DoGparams.bestFitValues(4)^2 * exp(-(pi*DoGparams.bestFitValues(4)*sfHiRes).^2) );
     theFittedSTF.surroundSTFHiRes = DoGparams.bestFitValues(1)*DoGparams.bestFitValues(2) * ( pi * (DoGparams.bestFitValues(4)*DoGparams.bestFitValues(3))^2 * exp(-(pi*DoGparams.bestFitValues(4)*DoGparams.bestFitValues(3)*sfHiRes).^2) );
     
end


function [theMidgetRGCMosaicResponses, spatialFrequenciesTested, spatialPhasesDegs] = ...
    computeTheSTF(theMidgetRGCmosaic, coneContrasts)

    sceneFOVdegs = theMidgetRGCmosaic.inputConeMosaic.sizeDegs;

    % Generate a presentation display with a desired resolution
    pixelsNum = 512;
    retinalImageResolutionDegs = max(sceneFOVdegs)/pixelsNum;
    viewingDistanceMeters = 4;
    theDisplay = rfMappingStimulusGenerator.presentationDisplay(...
            theMidgetRGCmosaic.inputConeMosaic.wave, retinalImageResolutionDegs, ...
            viewingDistanceMeters);

    % Stim params for the STF mapping
    stimParams = struct(...
            'backgroundLuminanceCdM2', 50.0, ...
            'backgroundChromaticity', [0.301 0.301], ...
            'coneContrasts', coneContrasts, ...
            'contrast', 0.75, ...
            'spatialFrequencyCPD', [], ...
            'orientationDegs', 0, ...
            'spatialPhaseIncrementDegs', 30, ...
            'pixelSizeDegs', retinalImageResolutionDegs, ...
            'stimSizeDegs', max(sceneFOVdegs), ...
            'wavelengthSupport', displayGet(theDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(theDisplay, 'viewing distance') ...
            );

    
    spatialFrequenciesTested = [0.25 0.5 1 2 4 6 8 10 12 16 24 32 64];

    rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsMicrons,1); 
    
    for iFreq = 1:numel(spatialFrequenciesTested)
   
        stimParams.spatialFrequencyCPD = spatialFrequenciesTested(iFreq);
        fprintf('Generating scenes for the frames of the %2.3f c/deg pattern.\n', stimParams.spatialFrequencyCPD);

        [theDriftingGratingSpatialModulationPatterns, spatialPhasesDegs] = ...
            rfMappingStimulusGenerator.driftingGratingFrames(stimParams);

        if (iFreq == 1)
            theMidgetRGCMosaicResponses = zeros(numel(spatialFrequenciesTested), numel(spatialPhasesDegs), rgcsNum);
        end

        % Generate scenes for the different spatial phasaes
        [theDriftingGratingFrameScenes, theNullStimulusScene, spatialSupportDegs] = ...
            rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                theDisplay, stimParams, theDriftingGratingSpatialModulationPatterns, ...
                'validateScenes', false);
   
        % Compute mRGCmosaic responses
        for iFrame = 1:numel(theDriftingGratingFrameScenes)
            theScene = theDriftingGratingFrameScenes{iFrame};

            r = theMidgetRGCmosaic.compute(...
                theScene, ...
                'nTrials', 1, ...
                'theNullScene', theNullStimulusScene, ...
                'normalizeConeResponsesWithRespectToNullScene', true);
            theMidgetRGCMosaicResponses(iFreq, iFrame,:) = r;
        end
        
    end

end


function RTVFTobjList = generateRTVFTobjects(theMidgetRGCmosaic, ...
    eccDegsGrid, conesNumPooledByTheRFcenterGrid, ...
    surroundToCenterRcRatioGrid, surroundToCenterIntegratedSensitivityRatioGrid)

    gridPositionsNum = size(eccDegsGrid,1);
    RTVFTobjList = cell(1, gridPositionsNum);

    % Visual RF model to match. Choose between: 
    % {'ellipsoidal gaussian center, gaussian surround', ...
    %  'gaussian center, gaussian surround', ...
    %  'arbitrary center, gaussian surround'}
    % When simulating the Croner&Kaplan assessment this must be set to 'gaussian center, gaussian surround';
    visualRFmodel = 'gaussian center, gaussian surround';

    % Retinal cone pooling model to use. Choose between:
    % 'arbitrary center cone weights, variable exponential surround weights';
    % 'arbitrary center cone weights, double exponential surround weights-free'
    % 'arbitrary center cone weights, double exponential surround weights-meanVnVwRatio'
    % 'arbitrary center cone weights, double exponential surround weights-meanRnRwRatio'
    % 'arbitrary center cone weights, double gaussian surround weights'
    % 'arbitrary center cone weights, gaussian surround weights'
    % 'arbitrary center cone weights, gaussian surround weights with adjustments'   % takes a long time - not very beneficial 
    % retinalConePoolingModel = 'arbitrary center cone weights, double exponential surround weights-meanRnRwRatio';
    retinalConePoolingModel = 'arbitrary center cone weights, double exponential surround weights-free';


    targetVisualRFDoGparams = struct(...
        'conesNumPooledByTheRFcenter', [], ...  % this will depend on the connectivity betwen cone/mRGC mosaics
        'surroundToCenterRcRatio', [], ...
        'surroundToCenterIntegratedSensitivityRatio', [], ... 
        'visualRFmodel', visualRFmodel, ...  
        'retinalConePoolingModel', retinalConePoolingModel ...
        );

    ZernikeDataBase = 'Artal2012';
    subjectRankOrder = 3;

    % Struct with the various optics params
    opticsParams = struct(...
        'positionDegs', [], ...           % (x,y) eccentricity for the PSF, in degrees
        'ZernikeDataBase', ZernikeDataBase, ...
        'examinedSubjectRankOrder', subjectRankOrder, ...
        'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
        'analyzedEye', 'right eye', ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', 3.0, ...
        'wavefrontSpatialSamples', [], ...
        'psfUpsampleFactor', [] ...
        );

    % Extract the cone mosaic from the midgetRGCmosaic
    theConeMosaic = theMidgetRGCmosaic.inputConeMosaic;

    for iGridPosition = 1:gridPositionsNum
        % Copy params structs
        theGridOpticsParams = opticsParams;
        theGridTargetVisualRFDoGparams = targetVisualRFDoGparams;
        
        % Update params structs for this grid position

        % Update opticsParams position for this grid position
        theGridOpticsParams.positionDegs = eccDegsGrid(iGridPosition,:);

        % Update psf upsample and spatial extend based on the eccentricity
        if (theGridOpticsParams.positionDegs <= 1)
            % Cones are tiny, so upsample the spatial resolution
            theGridOpticsParams.psfUpsampleFactor = 2;
            theGridOpticsParams.wavefrontSpatialSamples = 301;
        elseif (theGridOpticsParams.positionDegs <= 8)
            theGridOpticsParams.psfUpsampleFactor = 1;
            theGridOpticsParams.wavefrontSpatialSamples = 401;
        elseif (theGridOpticsParams.positionDegs <= 14)
            theGridOpticsParams.psfUpsampleFactor = 1;
            theGridOpticsParams.wavefrontSpatialSamples = 601;   
        else
            theGridOpticsParams.psfUpsampleFactor = 1;
            theGridOpticsParams.wavefrontSpatialSamples = 801;
        end

        % Update targetVisualRFDoGparams conesNum for this grid position
        theGridTargetVisualRFDoGparams.conesNumPooledByTheRFcenter = conesNumPooledByTheRFcenterGrid(iGridPosition);
        
        % Update targetVisualRFDoGparams surroundToCenterRcRatio for this grid position
        theGridTargetVisualRFDoGparams.surroundToCenterRcRatio = surroundToCenterRcRatioGrid(iGridPosition);

        % Update targetVisualRFDoGparams surroundToCenterIntegratedSensitivityRatio
        theGridTargetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio = surroundToCenterIntegratedSensitivityRatioGrid(iGridPosition);

        % Compute the RetinaToVisualFieldTransformer for this grid position
        tic
        
        % Select from:
        % - 1 (Single start run, fastest results), 
        % - some number (Multi-start), or 
        % - inf (Global search)

        multiStartsNum = 4;  
        doDryRunFirst = true;

        RTVFTobjList{iGridPosition} = RetinaToVisualFieldTransformer(...
            theConeMosaic, ...
            theGridOpticsParams, theGridTargetVisualRFDoGparams, ...
            'simulateCronerKaplanEstimation', true, ...
            'multiStartsNum', multiStartsNum, ...
            'doDryRunFirst', doDryRunFirst);
        toc

    end % iGridPosition

end