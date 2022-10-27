function generateMidgetRGCmosaicComponents
    
    % Operations
    operations = {...
        'generateMosaic' ...
        'visualizeMosaic' ...
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
    operations = operations(2:2);

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
   
    eccentricityIndices =  2:8 %[1 2 3 4 5 6 8 9 11 12];
    for ii = 1:numel(eccentricityIndices)
        
        iEcc = eccentricityIndices(ii);
        fprintf('Generating components for mosaic %d of %d\n', iEcc, numel(eccentricityIndices));
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

            case 'visualizeMosaic'
                % Load the midget RGCmosaic
                load(fName, 'theMidgetRGCmosaic');
                hFig = theMidgetRGCmosaic.visualize(...
                    'xRangeDegs', 0.75, ...
                    'yRangeDegs', 0.75, ...
                    'maxVisualizedRFs', 512);

                fNameMosaicPDF = strrep(fName, '.mat', '.pdf');
                NicePlot.exportFigToPDF(fNameMosaicPDF, hFig, 300);
                fprintf('Saved pdf of mosaic to %s\n', fNameMosaicPDF);

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
                    % Here we compute the temporal-equivalent eccentricity
                    % based SCint sensitivity ratio
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

                    % If we re-run this step a second time, the
                    % 'RTVFTobjList' variable does not exist in the file
                    % but we can get it from the midgetRGCMosaic itself
                    if (~exist('RTVFTobjList', 'var')) && (~isempty(theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList))
                        RTVFTobjList = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList;
                    end

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
                
                % Compute responses to a varietry of spatial frequencies and orientations
                [theMidgetRGCMosaicResponses, orientationsTested, spatialFrequenciesTested, spatialPhasesDegs] = ...
                    computeTheSTFs(theMidgetRGCmosaic, coneContrasts);
                
                % Assemble responses savefilename
                responsesPostfix = sprintf('_Responses_%2.2f_%2.2f_%2.2f.mat', ...
                    coneContrasts(1), coneContrasts(2), coneContrasts(3));
                fNameResponses = strrep(fName, '.mat', responsesPostfix);

                % Save the responses
                save(fNameResponses, ...
                    'theMidgetRGCMosaicResponses', ...
                    'orientationsTested', ...
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
                    'orientationsTested', ...
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

                    maxResponse = max(max(max(abs(squeeze(theMidgetRGCMosaicResponses(:, :, :, iRGC))))));
                    minResponse = -maxResponse;
                    meanResponse = 0;

                    % Generate circular LUT
                    oriColorLUT = phasemap(numel(orientationsTested));

                    theResponseModulation = zeros(numel(orientationsTested), numel(spatialFrequenciesTested));
                    for iSF = 1:numel(spatialFrequenciesTested)
                        for iOri = 1:numel(orientationsTested)
                            % Retrieve the mRGC response time-series
                            theResponse = squeeze(theMidgetRGCMosaicResponses(iOri, iSF, :, iRGC));
    
                            % Compute the response modulation for this SF
                            theResponseModulation(iOri, iSF) = max(theResponse)-min(theResponse);
                            
                            % Plot the time-series response for this SF and orientation
                            oriColor = oriColorLUT(iOri,:);
           
                            if (iOri == 1)
                                ax = subplot(numel(spatialFrequenciesTested),2,(iSF-1)*2+1);
                                plot(ax, 1:numel(spatialPhasesDegs), 0*theResponse, 'k-', 'LineWidth', 1.0);
                                hold(ax, 'on');
                            end

                            plot(ax, 1:numel(spatialPhasesDegs), theResponse, 'ko-', 'MarkerSize', 10, ...
                                'MarkerFaceColor', oriColor, 'LineWidth', 1.0);
                                
                            if (iOri == numel(orientationsTested))
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
    
                        end % iOri
                    end % iSF

                    theMeasuredSTF = theResponseModulation/max(theResponseModulation(:));

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
                    
                    % Plot the measured STF (for all orientations)
                    maxSF = zeros(1,numel(orientationsTested));
                    for iOri = 1:numel(orientationsTested)
                        theSTFatThisOri = squeeze(theMeasuredSTF(iOri,:));
                        plot(ax,spatialFrequenciesTested, theSTFatThisOri, '-', ...
                            'Color', oriColorLUT(iOri,:), 'LineWidth', 1.0);

                        % Find spatial frequency at which STF drops to 20%
                        % of max
                        [mag, iSFpeak] = max(theSTFatThisOri);
                        thresholdSTF = mag * 0.2;

                        ii = iSFpeak;
                        keepGoing = true;
                        iStop = [];
                        while (ii<numel(spatialFrequenciesTested))&&(keepGoing)
                            ii = ii + 1;
                            if (theSTFatThisOri(ii)>=thresholdSTF) && (theSTFatThisOri(ii+1)<thresholdSTF)
                                keepGoing = false;
                                iStop = ii;
                            end
                        end
                        maxSF(iOri) = spatialFrequenciesTested(iStop);

                    end

                    % Best orientation
                    [~, iBestOri] = max(maxSF);

                    p3 = plot(ax,spatialFrequenciesTested, theMeasuredSTF(iBestOri,:), 'ko', ...
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
                        [DoGparams{iFit}, theFittedSTF] = fitDoGmodelToMeasuredSTF(spatialFrequenciesTested, squeeze(theMeasuredSTF(iBestOri,:)), retinalRFcenterRcDegs);
                        
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
        'high', 100, ...
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
     multiStartsNum = 24;
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


function [theMidgetRGCMosaicResponses, orientationsTested, spatialFrequenciesTested, spatialPhasesDegs] = ...
    computeTheSTFs(theMidgetRGCmosaic, coneContrasts)

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

    
    rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsMicrons,1); 

    orientationsTested = 0:30:150;
    spatialFrequenciesTested = [0.25 0.5 1 2 4 6 8 10 12 16 24 32 64];

    % Allocate memory
    stimParams.orientationDegs = 0;
    stimParams.spatialFrequencyCPD = spatialFrequenciesTested(1);
    [~, spatialPhasesDegs] = rfMappingStimulusGenerator.driftingGratingFrames(stimParams);
    theMidgetRGCMosaicResponses = ...
        zeros(numel(orientationsTested), numel(spatialFrequenciesTested), numel(spatialPhasesDegs), rgcsNum);
     
    for iOri = 1:numel(orientationsTested)
        stimParams.orientationDegs = orientationsTested(iOri);

        fprintf('Computing STF for the %d degs orientation patterns.\n', ...
                stimParams.orientationDegs);
    
        parfor iFreq = 1:numel(spatialFrequenciesTested)
            theStimParams = stimParams;
            theStimParams.spatialFrequencyCPD = spatialFrequenciesTested(iFreq);
            
            
            % Generate spatial modulation patterns
            theDriftingGratingSpatialModulationPatterns = rfMappingStimulusGenerator.driftingGratingFrames(theStimParams);
    
            % Generate scenes for the different spatial phases
            [theDriftingGratingFrameScenes, theNullStimulusScene] = ...
                rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    theDisplay, theStimParams, theDriftingGratingSpatialModulationPatterns, ...
                    'validateScenes', false);
       
            % Allocate memory
            theFrameResponses = zeros(numel(spatialPhasesDegs), rgcsNum);

            % Compute mRGCmosaic responses
            for iFrame = 1:numel(spatialPhasesDegs)
                % Get scene corresponding to this stimulus frame
                theScene = theDriftingGratingFrameScenes{iFrame};
    
                % Compute the mosaic's response to this stimulus frame
                r = theMidgetRGCmosaic.compute(...
                    theScene, ...
                    'nTrials', 1, ...
                    'theNullScene', theNullStimulusScene, ...
                    'normalizeConeResponsesWithRespectToNullScene', true);
    
                % Store the mosaic's responses
                theFrameResponses(iFrame,:) = r;
            end

            % Save memory
            theDriftingGratingFrameScenes = [];
            theNullStimulusScene = [];

            theMidgetRGCMosaicResponses(iOri, iFreq,:,:) = theFrameResponses;
            
        end % iFreq
    end % iOri

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
            theGridOpticsParams, ...
            theGridTargetVisualRFDoGparams, ...
            'simulateCronerKaplanEstimation', true, ...
            'multiStartsNum', multiStartsNum, ...
            'doDryRunFirst', doDryRunFirst);
        toc

    end % iGridPosition

end


function M = phasemap(varargin)
% phasemap returns or sets a cyclic colormap with a constant lightness profile
% appropriate for plotting phase.  This is the phase colormap from the cmocean package,
% which is described here: 
% 
% Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 
% 2016. True colors of oceanography: Guidelines for effective and accurate 
% colormap selection. Oceanography 29(3):9?13. 
% http://dx.doi.org/10.5670/oceanog.2016.66
% 
%% Syntax 
% 
%  phasemap
%  M = phasemap(N)
%  phasemap(...,'rad') 
%  phasemap(...,'deg') 
% 
%% Description
% 
% phasemap sets the current colormap to a 256-level phase map.  
%
% M = phasemap(N) returns an Nx3 matrix of RGB values for a phase colormap. 
%
% phasemap(...,'rad') sets axis limits of current colormap to caxis([-pi pi]).
%
% phasemap(...,'deg') sets axis limits of current colormap to caxis([-180 180]).
% 
%% Example 1: Radians: 
% % Imagine some phase map ph in radians, which we'll wrap with phasewrap: 
% 
% ph = phasewrap(2*peaks(900)); 
% imagesc(ph) 
% colorbar
% phasemap
% 
%% Example 2: N-levels of degrees: 
% % Imagine some phase map ph in degrees, which we'll wrap with phasewrap: 
% 
% ph = phasewrap(100*peaks(900),'degrees'); 
% imagesc(ph) 
% colorbar
% phasemap(12)
% phasebar('location','se')
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas 
% at Austin's Institute for Geophysics (UTIG), May 2016. 
% http://www.chadagreene.com. 
% 
% RGB values are from Kristen Thyng's cmocean package, which is described 
% here: http://matplotlib.org/cmocean/.  We've also written about cmocean 
% for a peer-reviewed paper which came out in the journal Oceanography. If
% this colormap is useful for you, please do us a kindness and cite us: 
% 
% Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 
% 2016. True colors of oceanography: Guidelines for effective and accurate 
% colormap selection. Oceanography 29(3):9?13. 
% http://dx.doi.org/10.5670/oceanog.2016.66
% 
% See also phasewrap and phasebar. 
%% Set defaults: 
useRadians = true; 
N = 256; 
%% Parse inputs: 
if nargout==0
   
   % Try to automatically determine if current displayed data exist and are in radians or degrees: 
   if max(abs(caxis))>pi
      useRadians = false;
   else
      useRadians = true; 
   end
   
   % If no data are already displayed use radians: 
   if isequal(caxis,[0 1])
      useRadians = true; 
   end
end
tmp = isscalar(varargin); 
if any(tmp) 
   N = varargin{tmp}; 
end
tmp = strncmpi(varargin,'degrees',3); 
if any(tmp) 
   useRadians = false; 
end
%% Define RGB values
% These are from Kristen Thyng's cmocean package. 
M = [5.672433191540513508e-01 4.100182395068523444e-01 7.083093311283020221e-02
5.736953301772044211e-01 4.063138532087779420e-01 6.763640231601228403e-02
5.799480967447653290e-01 4.026276742797504915e-01 6.473393177099490137e-02
5.862172695866233463e-01 3.988330724695699492e-01 6.205458968303614181e-02
5.923144925108665415e-01 3.950439510462298576e-01 5.971205801157861842e-02
5.984111408418428812e-01 3.911540258605458176e-01 5.767411608227714787e-02
6.043424290206336202e-01 3.872685128597095350e-01 5.602885984725362173e-02
6.102730395075426362e-01 3.832795072149322024e-01 5.476296754995452903e-02
6.160285278727513081e-01 3.793045934655587259e-01 5.394141797415289569e-02
6.217946768383846301e-01 3.752152844285729927e-01 5.356325267543649710e-02
6.273981944803832533e-01 3.711339411256204812e-01 5.366088692600886473e-02
6.329634850050215311e-01 3.669706922951990746e-01 5.424568721283973965e-02
6.384108107505974683e-01 3.627845776164700942e-01 5.531730789006169408e-02
6.437627089659085922e-01 3.585593190882822023e-01 5.687403160453335244e-02
6.490458410279109636e-01 3.542732744485371232e-01 5.891870808260257081e-02
6.541890190896787471e-01 3.499849852073002587e-01 6.140802385263100654e-02
6.592788747503696145e-01 3.456223550919270138e-01 6.436753134213843430e-02
6.642461005174192801e-01 3.412442211863969677e-01 6.774174550344227996e-02
6.690819958208311657e-01 3.368604110213779812e-01 7.149835754351707706e-02
6.738554888085569461e-01 3.324082052078400351e-01 7.567487403043354766e-02
6.784905737890469801e-01 3.279583943075117247e-01 8.019026053099251317e-02
6.829905229231653108e-01 3.235105768774509949e-01 8.502570015729454811e-02
6.873973768055010591e-01 3.190238938082787246e-01 9.021504188680756764e-02
6.916667866641007523e-01 3.145441232419453059e-01 9.569860743379610124e-02
6.957922329372184800e-01 3.100809319167464606e-01 1.014568051402635063e-01
6.997776299670888100e-01 3.056329480851149860e-01 1.074875553773054770e-01
7.036357370628869567e-01 3.011877595522590711e-01 1.138098973510021372e-01
7.073365995938747375e-01 2.967824317970563230e-01 1.203755425082677633e-01
7.108777354295376938e-01 2.924241442509455946e-01 1.271782950774464904e-01
7.142560749513957585e-01 2.881211307727623572e-01 1.342125762305202696e-01
7.174739941560477341e-01 2.838747035534406327e-01 1.414868318931639668e-01
7.205242830202237547e-01 2.796988374089715013e-01 1.489903050379454141e-01
7.233972835735811291e-01 2.756122552678655913e-01 1.567027299209689251e-01
7.260902290803170622e-01 2.716252424744705074e-01 1.646182356269083324e-01
7.286005794465572061e-01 2.677482802858139888e-01 1.727304060861980239e-01
7.309261576887562395e-01 2.639918300158962117e-01 1.810321892738678551e-01
7.330653870302871189e-01 2.603659015050753256e-01 1.895163786061728473e-01
7.350186794657903588e-01 2.568773938658681089e-01 1.981820666380681439e-01
7.367827869904646221e-01 2.535402047817919002e-01 2.070088577247218842e-01
7.383581933944027842e-01 2.503622696871710129e-01 2.159876612538133012e-01
7.397459704801194746e-01 2.473506448008883085e-01 2.251089093005075314e-01
7.409477672978027618e-01 2.445114018731244232e-01 2.343626473565274537e-01
7.419657814271449769e-01 2.418495533170682554e-01 2.437386318766578963e-01
7.428027162250009363e-01 2.393690074290438785e-01 2.532264286379571905e-01
7.434617280753010871e-01 2.370725520147589083e-01 2.628155072505718448e-01
7.439463672555698404e-01 2.349618648405441546e-01 2.724953230491475309e-01
7.442604955662253108e-01 2.330377131867389928e-01 2.822544938256838498e-01
7.444083109771580942e-01 2.312996623822202169e-01 2.920823711509678011e-01
7.443942345178126141e-01 2.297462095038382279e-01 3.019690621658213669e-01
7.442228242253035031e-01 2.283749886783954519e-01 3.119048996322280609e-01
7.438987278365055689e-01 2.271828522038071840e-01 3.218804798008250700e-01
7.434266406175144004e-01 2.261659565030307850e-01 3.318866918844731173e-01
7.428112685586756303e-01 2.253198503030483524e-01 3.419147402717135353e-01
7.420572968700177574e-01 2.246395628281664947e-01 3.519561606537970344e-01
7.411693634997024160e-01 2.241196901870717295e-01 3.620028311945611144e-01
7.401524036663914563e-01 2.237545768080423936e-01 3.720435749650954360e-01
7.390107164863894962e-01 2.235380126256438782e-01 3.820735508345692866e-01
7.377485954521836309e-01 2.234637459623912692e-01 3.920865605716371549e-01
7.363702801567822975e-01 2.235254350456643002e-01 4.020760659682286464e-01
7.348798843643120637e-01 2.237166777905215120e-01 4.120358686970734263e-01
7.332813905831826462e-01 2.240310655045488208e-01 4.219601084098942301e-01
7.315786460893067833e-01 2.244622306283280966e-01 4.318432603708243778e-01
7.297753600299499155e-01 2.250038886803640281e-01 4.416801327482167139e-01
7.278757902192694029e-01 2.256496284295224086e-01 4.514624116253083685e-01
7.258826732888272737e-01 2.263936401540491139e-01 4.611893704869291510e-01
7.237992011402827330e-01 2.272301345525115646e-01 4.708571564055735181e-01
7.216284690463036222e-01 2.281535065561967612e-01 4.804619650485102977e-01
7.193734270864984293e-01 2.291583678228772403e-01 4.900002903581174851e-01
7.170368791833915401e-01 2.302395605311605142e-01 4.994689200594377154e-01
7.146214817984472001e-01 2.313921682581480121e-01 5.088649304139889473e-01
7.121300712110125719e-01 2.326113615814442426e-01 5.181844726895714626e-01
7.095645191802300022e-01 2.338929647122393596e-01 5.274270250454293762e-01
7.069269088332329476e-01 2.352328934084947365e-01 5.365908280501098249e-01
7.042192497772316040e-01 2.366272849868884531e-01 5.456740589098282301e-01
7.014433942735369687e-01 2.380725388567397904e-01 5.546751290266969114e-01
6.986010335708933150e-01 2.395653172009881704e-01 5.635926726134320441e-01
6.956936866865383040e-01 2.411025487728192251e-01 5.724255550539750770e-01
6.927222115248060641e-01 2.426816864587864053e-01 5.811742593820377056e-01
6.896879229931356381e-01 2.443001105764274206e-01 5.898374988574298650e-01
6.865918014911074341e-01 2.459555424906482712e-01 5.984146270185066729e-01
6.834346456921173152e-01 2.476459602268146765e-01 6.069051304628041432e-01
6.802170668044646984e-01 2.493695974967287121e-01 6.153086083569707654e-01
6.769394825945098670e-01 2.511249427769637332e-01 6.236247501806896354e-01
6.736020505347796172e-01 2.529107709208097687e-01 6.318534591911313392e-01
6.702032268029198825e-01 2.547269081210946085e-01 6.399981954948240626e-01
6.667440181482768846e-01 2.565719542580951473e-01 6.480556766297400628e-01
6.632240032303863275e-01 2.584454084861053658e-01 6.560256306953086147e-01
6.596425384904228695e-01 2.603470198353451392e-01 6.639077331315591524e-01
6.559987514287347610e-01 2.622767870376229715e-01 6.717015744726071436e-01
6.522915338733621393e-01 2.642349583566807047e-01 6.794066263454738852e-01
6.485195352742859631e-01 2.662220313649069636e-01 6.870222057113775094e-01
6.446806025240062743e-01 2.682390430547165194e-01 6.945485065851293438e-01
6.407712044067549462e-01 2.702878635294215792e-01 7.019874685983423790e-01
6.367911507125239012e-01 2.723687210815862736e-01 7.093338316074968564e-01
6.327380807122843231e-01 2.744830916476992400e-01 7.165858102551939668e-01
6.286093641096192064e-01 2.766326937674883912e-01 7.237412561247542619e-01
6.244020959215076383e-01 2.788194857486135936e-01 7.307976115653572746e-01
6.201130919933868224e-01 2.810456617826825876e-01 7.377518617228360220e-01
6.157388853340158841e-01 2.833136467579401296e-01 7.446004848213658711e-01
6.112757235007779677e-01 2.856260894714388043e-01 7.513394007827567389e-01
6.067195673197314232e-01 2.879858538961362902e-01 7.579639183289200721e-01
6.020660912878182947e-01 2.903960081063858834e-01 7.644686807933391837e-01
5.973106860777335214e-01 2.928598104095548882e-01 7.708476109730414416e-01
5.924478666488723899e-01 2.953810018358738465e-01 7.770946054703944395e-01
5.874718273398599200e-01 2.979635029516978784e-01 7.832026524810360435e-01
5.823786412617427688e-01 3.006102519317957467e-01 7.891612643940841831e-01
5.771627200394546797e-01 3.033250182270930684e-01 7.949608965741956634e-01
5.718182588451700132e-01 3.061116287038831429e-01 8.005910315557441814e-01
5.663392645516818202e-01 3.089739250027142559e-01 8.060401374993888535e-01
5.607195914193257025e-01 3.119157138110425498e-01 8.112956338568207970e-01
5.549529855475957563e-01 3.149407094165454501e-01 8.163438666523304965e-01
5.490331394145471222e-01 3.180524679993645409e-01 8.211700962113627211e-01
5.429537578898111505e-01 3.212543132710714011e-01 8.257585005973195891e-01
5.367081743584852793e-01 3.245494976399965759e-01 8.300925055155035093e-01
5.302876326058616474e-01 3.279420696272150604e-01 8.341557801338073119e-01
5.236887893045435449e-01 3.314328654907214844e-01 8.379276125633778882e-01
5.169063460899447904e-01 3.350233559526123450e-01 8.413881804452462143e-01
5.099355929853703895e-01 3.387142672664584242e-01 8.445169634123261826e-01
5.027725725347598207e-01 3.425054625546620768e-01 8.472929342850022971e-01
4.954142602591064537e-01 3.463958234006226933e-01 8.496947959170987330e-01
4.878587595988513326e-01 3.503831358364098980e-01 8.517012656547067184e-01
4.801055083588288697e-01 3.544639856683631796e-01 8.532914075636400808e-01
4.721554923719734620e-01 3.586336687735280915e-01 8.544450102471590203e-01
4.640044749743072461e-01 3.628897608659185647e-01 8.551433986394967324e-01
4.556622651494331433e-01 3.672221012185919453e-01 8.553678849867077938e-01
4.471383276851411126e-01 3.716204653438348049e-01 8.551029025080878476e-01
4.384422744277990391e-01 3.760744076888676291e-01 8.543357673497274929e-01
4.295862071986853437e-01 3.805721840819369373e-01 8.530568872354411525e-01
4.205817729540208272e-01 3.851023944161593993e-01 8.512595016023745131e-01
4.114469045503038047e-01 3.896507488703780386e-01 8.489411196051931396e-01
4.022060497114997024e-01 3.942002477832268204e-01 8.461048966175331865e-01
3.928817094017691969e-01 3.987355949127092125e-01 8.427576180215047286e-01
3.834985523654112494e-01 4.032413019893585915e-01 8.389104772668918297e-01
3.740836020579017540e-01 4.077017739713900135e-01 8.345792042879973671e-01
3.646672169663981444e-01 4.121010163064194876e-01 8.297846549404398475e-01
3.552786781091737400e-01 4.164248241297358044e-01 8.245505460684040555e-01
3.459478747064330673e-01 4.206600190520577986e-01 8.189037829436274230e-01
3.367081015244690612e-01 4.247932282895754974e-01 8.128760845516961320e-01
3.275908439376793990e-01 4.288132465017488459e-01 8.065006765026215829e-01
3.186221457587909978e-01 4.327124829824056107e-01 7.998088272458536707e-01
3.098298271614847166e-01 4.364836390012560297e-01 7.928343311444766561e-01
3.012402910763991026e-01 4.401209996756487719e-01 7.856111406094337113e-01
2.928823313409779638e-01 4.436186871784305041e-01 7.781766034977661839e-01
2.847793948471887826e-01 4.469737539371314017e-01 7.705648357916961011e-01
2.769467007747948850e-01 4.501870831571345710e-01 7.628028303724243564e-01
2.694019299196109829e-01 4.532587645627866313e-01 7.549201981024484809e-01
2.621602390310989183e-01 4.561900001239682090e-01 7.469447502290214036e-01
2.552341461123000532e-01 4.589829673886487993e-01 7.389023362912160442e-01
2.486334450106276739e-01 4.616406826381265760e-01 7.308167439984587510e-01
2.423662314394859907e-01 4.641664308850603216e-01 7.227110939283054591e-01
2.364388695871720847e-01 4.665636310323913460e-01 7.146082501062388515e-01
2.308485415758319559e-01 4.688387341053135704e-01 7.065202838204042157e-01
2.255938839599332946e-01 4.709968473866177896e-01 6.984625903578935979e-01
2.206708265023217819e-01 4.730433385072548291e-01 6.904486767769802968e-01
2.160726675659729312e-01 4.749837571162059402e-01 6.824902568657742474e-01
2.117901806110671359e-01 4.768237675084572857e-01 6.745973531916734656e-01
2.078117520620332837e-01 4.785690914533750617e-01 6.667784017923048534e-01
2.041235496378702718e-01 4.802254602742371259e-01 6.590403559750543927e-01
2.007097189181233365e-01 4.817985751209007961e-01 6.513887865326223325e-01
1.975526046085452347e-01 4.832940743336505141e-01 6.438279764300510744e-01
1.946329918332628250e-01 4.847175068012891663e-01 6.363610086687176981e-01
1.919303619466051491e-01 4.860743102570944107e-01 6.289898465867946875e-01
1.894231569164997186e-01 4.873697935214787447e-01 6.217154063202443570e-01
1.870890463166731688e-01 4.886091217818528398e-01 6.145376215326621150e-01
1.849051913577256789e-01 4.897973040922220522e-01 6.074555008360708053e-01
1.828485011223850432e-01 4.909391823733177396e-01 6.004671785775045345e-01
1.808958771515883079e-01 4.920394212959531832e-01 5.935699598659075482e-01
1.790244436455807608e-01 4.931024985341831424e-01 5.867603608667124604e-01
1.772117616891990266e-01 4.941326949799630941e-01 5.800341455014857983e-01
1.754360269913468329e-01 4.951340846166410459e-01 5.733863597592738781e-01
1.736762515758797087e-01 4.961105238543395690e-01 5.668113648548980255e-01
1.719124306347093734e-01 4.970656402350968617e-01 5.603028704563300932e-01
1.701256963402829336e-01 4.980028205181701861e-01 5.538539691464678949e-01
1.682986669160662918e-01 4.989250963808533057e-01 5.474578820340219032e-01
1.664149064733197181e-01 4.998354730877956897e-01 5.411056294295618629e-01
1.644597156463551424e-01 5.007365672557706482e-01 5.347884863149882095e-01
1.624201076780373842e-01 5.016306791238075435e-01 5.284974569553574364e-01
1.602848864385326977e-01 5.025198058433570925e-01 5.222231535388160983e-01
1.580447765919143488e-01 5.034056306727768826e-01 5.159558533581567463e-01
1.556925602346709570e-01 5.042895136618935625e-01 5.096855603565535464e-01
1.532232234951459993e-01 5.051724842338358723e-01 5.034020703832871035e-01
1.506341170297759435e-01 5.060552360210837097e-01 4.970950392165090492e-01
1.479251349064662491e-01 5.069381242363323326e-01 4.907540521606342243e-01
1.450989169821367275e-01 5.078211657594990935e-01 4.843686938332767977e-01
1.421610804711274656e-01 5.087040420049282474e-01 4.779286166356862542e-01
1.391204868117350002e-01 5.095861045039652426e-01 4.714236063620425576e-01
1.359895499317730627e-01 5.104663830048005169e-01 4.648436434523450234e-01
1.327845912322146715e-01 5.113435957613400307e-01 4.581789585308634893e-01
1.295262445313372923e-01 5.122161615629505782e-01 4.514200810916398354e-01
1.262399101195911255e-01 5.130822129532254250e-01 4.445578804842563181e-01
1.229555631631381596e-01 5.139397889872834302e-01 4.375821217122248252e-01
1.197103614797271343e-01 5.147863054564670859e-01 4.304858794552847701e-01
1.165471638531235532e-01 5.156190876837056791e-01 4.232614510628861515e-01
1.135154897666945939e-01 5.164352388982441644e-01 4.159014793262920673e-01
1.106718285206290120e-01 5.172316368920498730e-01 4.083990629543874373e-01
1.080796164815457927e-01 5.180049408378688547e-01 4.007477580736906742e-01
1.058087845385388415e-01 5.187515945922266392e-01 3.929415772123373007e-01
1.039347151626626997e-01 5.194678258709898300e-01 3.849749882596104622e-01
1.025359694379740783e-01 5.201499324230609567e-01 3.768393027599689638e-01
1.016935492416111697e-01 5.207935695646238594e-01 3.685305536972853790e-01
1.014861636581259607e-01 5.213941640446873027e-01 3.600451225619441531e-01
1.019864741263540875e-01 5.219469371865886886e-01 3.513789710355462725e-01
1.032573879600174327e-01 5.224467943632038480e-01 3.425286138069295050e-01
1.053486489718466201e-01 5.228882869876148032e-01 3.334912134806259920e-01
1.082944758843455979e-01 5.232655680556403954e-01 3.242647166459451946e-01
1.121126465502676472e-01 5.235723416856985502e-01 3.148480421610331259e-01
1.168080773882203272e-01 5.238019157845009710e-01 3.052357609469624200e-01
1.223851778388554135e-01 5.239470202573554003e-01 2.954047708042208953e-01
1.288112281029662742e-01 5.239987576746931719e-01 2.853827900944919116e-01
1.360562084538863015e-01 5.239482011263842942e-01 2.751760306766347641e-01
1.441051441567454128e-01 5.237851259274813875e-01 2.647689845875647596e-01
1.529337458740115119e-01 5.234977984931444839e-01 2.541648013875384415e-01
1.624797709886042218e-01 5.230749435633894606e-01 2.434144547473032216e-01
1.727628202070501939e-01 5.225005508142598343e-01 2.324837423146093873e-01
1.836989343202152236e-01 5.217627536881318528e-01 2.214609776488171700e-01
1.952977171051294425e-01 5.208432064756289837e-01 2.103468003936530550e-01
2.074630833256147500e-01 5.197310268323230842e-01 1.992559173422548402e-01
2.201804866465039145e-01 5.184088461339138032e-01 1.882360891676425896e-01
2.333258146567612767e-01 5.168714604143809233e-01 1.774406118077297145e-01
2.467881037119215581e-01 5.151165689450147855e-01 1.670111936474235337e-01
2.604357846291780465e-01 5.131495394379788078e-01 1.571028965071679384e-01
2.740979295754801814e-01 5.109884830669607636e-01 1.478869433321773896e-01
2.876188391824692214e-01 5.086577674346561828e-01 1.395008019535278476e-01
3.008697529112238644e-01 5.061850478527918362e-01 1.320370980173991848e-01
3.137048873182943232e-01 5.036086579898598758e-01 1.255636295957552051e-01
3.260953799514536011e-01 5.009490735705195430e-01 1.200560783434582723e-01
3.379942414546169838e-01 4.982324891216243778e-01 1.154774348588594179e-01
3.493848946613353212e-01 4.954799141616709757e-01 1.117578580398684851e-01
3.602962205045775468e-01 4.927011375065428189e-01 1.088017440281668080e-01
3.707557362021778880e-01 4.899046162256747716e-01 1.065118926400932531e-01
3.807631298544894016e-01 4.871051876936182690e-01 1.047958357127047602e-01
3.903962780488588469e-01 4.842942108846688964e-01 1.035513612645985182e-01
3.996816423057455991e-01 4.814750885349466381e-01 1.026934005733520505e-01
4.086118231999155692e-01 4.786607467967835539e-01 1.021446830610131806e-01
4.172551136119685977e-01 4.758391896276191191e-01 1.018301287096897711e-01
4.256394099183998248e-01 4.730091323278173276e-01 1.016864378549218739e-01
4.337892663754943090e-01 4.701691892705471276e-01 1.016574178896643976e-01
4.417279681251333190e-01 4.673173285719916525e-01 1.016932874387579400e-01
4.494774245867361739e-01 4.644510027208682001e-01 1.017499507788620638e-01
4.570581245596661568e-01 4.615672570884231107e-01 1.017881641465586406e-01
4.644891347214414323e-01 4.586628194897153787e-01 1.017726540180271533e-01
4.717881280670202515e-01 4.557341737961458361e-01 1.016712192056216102e-01
4.789714327661078186e-01 4.527776202342287948e-01 1.014538260631735189e-01
4.860540949604466587e-01 4.497893245680554819e-01 1.010916883277608036e-01
4.930499512220878899e-01 4.467653579561612909e-01 1.005563084901952564e-01
4.999717077592871761e-01 4.437017289996793012e-01 9.981844399132347445e-02
5.068310239546541807e-01 4.405944094751862994e-01 9.884694646914210514e-02
5.136385973536989891e-01 4.374393556236486624e-01 9.760740246030363831e-02
5.204042456108760595e-01 4.342325278539159172e-01 9.606047462512840029e-02
5.271369778483794288e-01 4.309699136110089479e-01 9.415979602623703038e-02
5.338450429270207298e-01 4.276475613850701518e-01 9.184919276569031288e-02
5.405359345634145329e-01 4.242616390225810163e-01 8.905887705893447692e-02
5.472163220206374135e-01 4.208085375489562519e-01 8.570001281033132190e-02
5.538918588348998862e-01 4.172850538650798025e-01 8.165659961692531277e-02
5.605667990303634385e-01 4.136887034926055584e-01 7.677269383879334330e-02
5.672433191540513508e-01 4.100182395068523444e-01 7.083093311283020221e-02];
%% Interpolate if necessary: 
% Because this map is cyclic, getting even spacing for the wrap-around requires
% adding an N for interpolation, then discarding the last value.   
if N~=256
   R = interp1(1:256,M(:,1),linspace(1,256,N+1)'); 
   G = interp1(1:256,M(:,2),linspace(1,256,N+1)'); 
   B = interp1(1:256,M(:,3),linspace(1,256,N+1)'); 
   M = [R(1:end-1) G(1:end-1) B(1:end-1)]; 
end
   
end
