function generateMidgetRGCmosaicComponents
    
    % Operations
    operations = {...
        'generateMosaic' ...
        'generateRTVFTobjList' ...
        'generateCenterSurroundRFs' ...
        'computeSTF' ...
        'visualizeResponses' ...
        };

    % Operation to compute
    operations = {'visualizeResponses'};
    %coneContrasts = [0.12 -0.12 0];
    coneContrasts = [1 1 0];

    eccSizeDegsExamined = [...
         0  0.3; ...
        -1  0.4; ...
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

    for iEcc = 3:3 % 3:size(eccSizeDegsExamined,1)
        fprintf('Generating components for mosaic %d of %d\n', iEcc, size(eccSizeDegsExamined,1));
        eccDegs  = eccSizeDegsExamined(iEcc,1) * [1 0];
        sizeDegs = eccSizeDegsExamined(iEcc,2) * [1 1];
        doIt(operations, eccDegs, sizeDegs, coneContrasts);
    end

end

function doIt(operations, eccDegs, sizeDegs, coneContrasts)

    fName = sprintf('mRGCmosaicComponents_eccDegs_%2.2f.mat', eccDegs(1));

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
                    surroundToCenterRcRatioGrid(iGridPosition) = 6.7;
    
                    % Also from Croner & Kaplan '95 (Figure 10b)
                    % "These mean ratios for P and M cells are not significantly different
                    % (Student's t-test: P = 0.482). The overall mean ratio is 0.55.
                    surroundToCenterIntegratedSensitivityRatioGrid(iGridPosition) = 0.55;
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

                    % Generate the C/S spatial RF
                    theMidgetRGCmosaic.generateCenterSurroundSpatialPoolingRF(RTVFTobjList, ...
                        eccDegsGrid, conesNumPooledByTheRFcenterGrid, ...
                        surroundToCenterRcRatioGrid, surroundToCenterIntegratedSensitivityRatioGrid);

                    save(fName, 'theMidgetRGCmosaic', '-append');

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

            case 'visualizeResponses'
                load(fName, 'theMidgetRGCmosaic', 'RTVFTobjList');
                RTVFTobjList{1}
                RTVFTobjList{1}.rfComputeStruct
                
                RTVFTobjList{1}.targetVisualRFDoGparams
                
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
                fNameVideo = strrep(fName, '.mat', '_Video');
                videoOBJ = VideoWriter(fNameVideo, 'MPEG-4');
                videoOBJ.FrameRate = 10;
                videoOBJ.Quality = 100;
                videoOBJ.open();

                for iii = 1:numel(sortedRGCindices)
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

                    % Visualize the examined retinal RF
                    ax = subplot(numel(spatialFrequenciesTested),2, [2 4 6 8]);
                    theMidgetRGCmosaic.visualizeSingleRetinalRF(iRGC, ...
                        'plotTitle', sprintf('RGC: %d of %d', iii, numel(sortedRGCindices)), ...
                        'figureHandle', hFig, ...
                        'axesHandle', ax);

                    % Visualize the computed STF
                    ax = subplot(numel(spatialFrequenciesTested),2, ((6:numel(spatialFrequenciesTested))-1)*2+2);
                    normalizedTargetSTF = RTVFTobjList{1}.rfComputeStruct.theSTF.target;
                    normVal = max(normalizedTargetSTF);
                    normalizedTargetSTF = normalizedTargetSTF / normVal;

                    % Plot the target STF
                    plot(ax, theRTVFTobj.rfComputeStruct.theSTF.support, normalizedTargetSTF, 'r-', 'Color', [1 0.5 0.5], 'LineWidth', 3.0);
                    hold(ax, 'on')
                    p1 = plot(ax, theRTVFTobj.rfComputeStruct.theSTF.support, normalizedTargetSTF, 'r-', 'LineWidth', 1.5);
                    
                    % Plot the measured STF
                    p2 = plot(ax,spatialFrequenciesTested, theResponseModulation/max(theResponseModulation), 'ko', ...
                        'MarkerSize', 14, 'MarkerFaceColor', [0.2 0.9 0.9], 'MarkerEdgeColor', [0 0.4 1], ...
                        'LineWidth', 1.0);

                    hold(ax, 'on');

                    hold(ax, 'off');
                    legend([p1, p2], {'target', 'measured'}, 'Location', 'SouthWest');
                    title(ax, sprintf('stim LMS contrast = < %2.1f, %2.1f, %2.1f >', coneContrasts(1), coneContrasts(2), coneContrasts(3)));
                    set(ax, 'XLim', [0.3 70], 'XTick', [0.1 0.3 1 3 10 30 100], 'YLim', [0 1.02], 'YTick', 0:0.1:1.0);
                    xlabel(ax, 'spatial frequency (c/deg)');
                    ylabel(ax, 'STF');
                    grid(ax, 'on');
                    set(ax, 'XLim', [0.1 100], 'XScale', 'Log', 'FontSize', 16);
                    
                    drawnow;
                    videoOBJ.writeVideo(getframe(hFig));
                end
                videoOBJ.close();


            otherwise
                error('Unknown operation: ''%s''.', operations{iOp});
        end % Switch
    end
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
        'wavefrontSpatialSamples', 701, ...
        'psfUpsampleFactor', 2 ...
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

        % Update targetVisualRFDoGparams conesNum for this grid position
        theGridTargetVisualRFDoGparams.conesNumPooledByTheRFcenter = conesNumPooledByTheRFcenterGrid(iGridPosition);
        
        % Update targetVisualRFDoGparams surroundToCenterRcRatio for this grid position
        theGridTargetVisualRFDoGparams.surroundToCenterRcRatio = surroundToCenterRcRatioGrid(iGridPosition);

        % Update targetVisualRFDoGparams surroundToCenterIntegratedSensitivityRatio
        theGridTargetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio = surroundToCenterIntegratedSensitivityRatioGrid(iGridPosition);

        % Compute the RetinaToVisualFieldTransformer for this grid position
        tic
        multiStartsNum = 16;
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
