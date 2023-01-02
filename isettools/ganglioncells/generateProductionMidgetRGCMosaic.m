function generateProductionMidgetRGCMosaic()

    mosaicCenterParams = struct(...
        'positionDegs',[0 0], ...
        'sizeDegs',  [3 3], ...        
        'whichEye', 'right eye');

    H1cellIndex = 1;

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = generateMosaicFileName(mosaicCenterParams);
    

    % Actions to perform
    actionToPerform = 'generateCenterConnectedMosaic';
    actionToPerform = 'generateR2VFTobjects';
    actionToPerform = 'reGenerateR2VFTobjectsAtSpecificPosition'
    %actionToPerform = 'inspectR2VFTobjects';
    %actionToPerform = 'generateCenterSurroundRFs';
   
    %actionToPerform = 'visualizeFittedSpatialRFs';
    %actionToPerform = 'computeMosaicSTFs';
    %actionToPerform = 'fitMosaicSTFs';

    switch (actionToPerform)

        case 'visualizeFittedSpatialRFs'
            rgcIndicesToAnalyze = generateListOfTargetRGCs(mosaicFileName);
            visualizeFittedSpatialRFmaps(mosaicFileName, rgcIndicesToAnalyze);

        case 'fitMosaicSTFs'
            responsesFileName = responsesFileNameForMosaicFileName(mosaicFileName, H1cellIndex);
            rgcIndicesToAnalyze = generateListOfTargetRGCs(mosaicFileName);
            fitMosaicSTFs(mosaicFileName, responsesFileName, rgcIndicesToAnalyze);

        case 'computeMosaicSTFs'
            responsesFileName = responsesFileNameForMosaicFileName(mosaicFileName, H1cellIndex);
            computeMosaicAchromaticSTFs(mosaicFileName, responsesFileName);

        case 'generateCenterSurroundRFs'
            R2VFTobjFileName = R2VFTobjFileNameForMosaicFileName(mosaicFileName, H1cellIndex);
            generateCenterSurroundRFs(mosaicFileName, R2VFTobjFileName);

        case 'inspectR2VFTobjects'
            dropboxDir = localDropboxPath();
            [file,path] = uigetfile(fullfile(localDropboxPath(), '*.mat'), ...
                        'Select an RTVF file');
            inspectSavedRTVFfile(fullfile(path,file));

        case {'reGenerateR2VFTobjectsAtSpecificPosition','generateR2VFTobjects'}

            mosaicSurroundParams = struct(...
                'eccentricitySamplingGridHalfSamplesNum', 1, ...                         % generate R2VFTobjects at 2*gridHalfSamplesNum + 1 spatial positions
                'centerConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ...   % cone types that can connect to the RF center
                'surroundConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ... % cone types that can connect to the RF surround
                'coneWeightsCompensateForVariationsInConeEfficiency', true, ...          % cone weight compensationfor eccentricity-dependent variations in cone efficiency
                'visualRFmodel', 'gaussian center, gaussian surround', ...
                'retinalConePoolingModel', sprintf('arbitrary center cone weights, double exponential surround from H1 cell with index %d', H1cellIndex),...
                'H1cellIndex', H1cellIndex, ...
                'targetSTFmatchMode', 'STFDoGparams' ...
            );
        
            opticsParams = struct(...
                'ZernikeDataBase', 'Polans2015', ...
                'subjectRankOrder', 6, ...
                'pupilDiameterMM', 3.0 ...
            );

            if (strcmp(actionToPerform, 'generateR2VFTobjects'))
                generateR2VFTobjects(mosaicCenterParams, mosaicSurroundParams, opticsParams, mosaicFileName, mosaicDirectory);
            else
                targetPosition = input('Enter position for which to update the RTVF object ([x y]): ');
                targetRFcenterConesNum = input('Enter center cones num for which to update the RTVF object (e.g, 1, 2): ');
                generateR2VFTobjects(mosaicCenterParams, mosaicSurroundParams, opticsParams, mosaicFileName, mosaicDirectory, ...
                    'updateRTVFobjectAtPosition', targetPosition, ...
                    'updateRTVFobjectWithCenterConesNum', targetRFcenterConesNum);
            end

        case 'generateCenterConnectedMosaic'
            generateCenterConnectedMosaic(mosaicCenterParams, mosaicFileName);

        otherwise
            error('Unknown action: ''%s''.', actionToPerform);
    end
end


function rgcIndicesToAnalyze = generateListOfTargetRGCs(mosaicFileName)
    % Load the fully-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    rgcsNum = size(theMidgetRGCmosaic.rgcRFsurroundConePoolingMatrix,2);
    maxVisualizedRFs = 1000;
    skippedRGCs = max([1 round(rgcsNum/maxVisualizedRFs)]);
    
    % Sort RGCs with respect to their distance from the mRGC mosaic center
    mRGCmosaicCenterDegs = mean(theMidgetRGCmosaic.rgcRFpositionsDegs,1);
    d = sum((bsxfun(@minus, theMidgetRGCmosaic.rgcRFpositionsDegs, mRGCmosaicCenterDegs)).^2,2);
    [~,sortedRGCindices] = sort(d, 'ascend');

    rgcIndicesToAnalyze = sortedRGCindices(1:skippedRGCs:rgcsNum);
end


function fitMosaicSTFs(mosaicFileName, responsesFileName, rgcIndicesToAnalyze)

    % Load the fully-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    % Load the mosaic responses
    load(responsesFileName, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts');

    % Allocate memory
    fittedSTFs = cell(1, numel(rgcIndicesToAnalyze));

    for iRGC = 1:numel(rgcIndicesToAnalyze)
        % Target RGC
        theRGCindex = rgcIndicesToAnalyze(iRGC);
        fprintf('Fitting RGC %d of %d, located at (%2.2f,%2.2f degs)\n', ...
            iRGC, numel(rgcIndicesToAnalyze), theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,1), theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,2));

        % Obtain the responses
        theSingleMidgetRGCResponses = squeeze(theMidgetRGCMosaicResponses(:, :, :, theRGCindex));

        % Select STF to fit (orientation that results in max extension into high SFs)
        [theMeasuredSTF, allMeasuredSTFs] = STFtoFit(theSingleMidgetRGCResponses, spatialFrequenciesTested, orientationsTested);

        % Fit the DoG model to the measured STF
        multiStartsNum = 128;
        [theFittedDoGmodelParams, theDoGmodelFitOfTheMeasuredSTF] = ...
            RetinaToVisualFieldTransformer.fitDoGmodelToMeasuredSTF(spatialFrequenciesTested, ...
                    theMeasuredSTF, ...
                    retinalRFcenterRcDegsInitialEstimate(theMidgetRGCmosaic, theRGCindex), ...
                    multiStartsNum);


        [triangulatingRTVFobjIndices, triangulatingRTVFobjWeights] = ...
             theMidgetRGCmosaic.triangulatingRTVFobjectIndicesAndWeights(theRGCindex);

        pipelineScaleFactorBasedOnLowestSF = 0;
        triangulatingRTVFobjSTFdata = cell(1, numel(triangulatingRTVFobjIndices));

        for iObj = 1:numel(triangulatingRTVFobjIndices)
            theRTVFobjIndex = triangulatingRTVFobjIndices(iObj);
            theRTVFTobj = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{theRTVFobjIndex};

            % Normalize theRTVFTobj.fitted to the theRTVFTobj.targetSTF
            scaleFactorBasedOnLowestSF = theRTVFTobj.rfComputeStruct.theSTF.target(1)/theRTVFTobj.rfComputeStruct.theSTF.fitted(1);
            pipelineScaleFactorBasedOnLowestSF = pipelineScaleFactorBasedOnLowestSF + theRTVFTobj.rfComputeStruct.theSTF.target(1)*triangulatingRTVFobjWeights(iObj);
            
            % The spatial support
            triangulatingRTVFobjSTFdata{iObj} = struct();
            triangulatingRTVFobjSTFdata{iObj}.spatialFrequencySupport = theRTVFTobj.rfComputeStruct.theSTF.support(:);
            % The model-achieved STF
            triangulatingRTVFobjSTFdata{iObj}.fittedSTF = theRTVFTobj.rfComputeStruct.theSTF.fitted(:)*scaleFactorBasedOnLowestSF;
            % The model-target STF
            triangulatingRTVFobjSTFdata{iObj}.targetSTF = theRTVFTobj.rfComputeStruct.theSTF.target(:);
        end

        % Normalize theRTVFTobj.fitted to the theRTVFTobj.targetSTF
        pipelineScaleFactorBasedOnLowestSF = pipelineScaleFactorBasedOnLowestSF /theDoGmodelFitOfTheMeasuredSTF.compositeSTFHiRes(1);

        figure(10000); clf;
        for iObj = 1:numel(triangulatingRTVFobjIndices) 

            RTVFobjPosition = theMidgetRGCmosaic.theSamplingPositionGrid(triangulatingRTVFobjIndices(iObj),:);

            subplot(1,numel(triangulatingRTVFobjIndices), iObj)
            plot(triangulatingRTVFobjSTFdata{iObj}.spatialFrequencySupport, triangulatingRTVFobjSTFdata{iObj}.targetSTF, ...
                'k-', 'LineWidth', 1.5);
            hold on;
            plot(triangulatingRTVFobjSTFdata{iObj}.spatialFrequencySupport, triangulatingRTVFobjSTFdata{iObj}.fittedSTF, ...
                'r-', 'LineWidth', 1.5);
            plot(spatialFrequenciesTested, theMeasuredSTF*pipelineScaleFactorBasedOnLowestSF, 'co-', 'MarkerSize', 14, 'LineWidth', 1.5, 'MarkerFaceColor', [0 1 1], 'Color', [0 0.6 0.6]);
            plot(theDoGmodelFitOfTheMeasuredSTF.sfHiRes, theDoGmodelFitOfTheMeasuredSTF.compositeSTFHiRes*pipelineScaleFactorBasedOnLowestSF, 'b-', 'LineWidth', 1.5);
            legend({'RTVF: target', 'RTVF: fitted', 'pipeline: measured', 'pipeline: DoGmodel'}, 'Location', 'SouthWest');
            set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], 'XLim', [0.01 100]);
            grid on
            set(gca, 'FontSize', 16);
            xlabel('spatial frequency (c/deg)')
            title(sprintf('component RTVF model\nposition:(%2.2f,%2.2f degs), weight: %2.3f', ...
                RTVFobjPosition(1), RTVFobjPosition(2), triangulatingRTVFobjWeights(iObj)));
        end
        

        fittedSTFs{iRGC} = struct(...
            'targetRGC', theRGCindex, ...
            'theMultiFocalRTVFmodelSTFdata', triangulatingRTVFobjSTFdata, ...
            'theMultiFocalRTVFmodelWeights', triangulatingRTVFobjWeights, ...
            'spatialFrequencySupport', spatialFrequenciesTested, ...
            'theMeasuredSTF', theMeasuredSTF, ... 
            'allMeasuredSTFs', allMeasuredSTFs, ...
            'theDoGmodelFitOfTheMeasuredSTF', theDoGmodelFitOfTheMeasuredSTF, ...
            'theFittedDoGmodelParams', theFittedDoGmodelParams);
            
    end % iRGC

    % Append the fittedSTFs structs
    save(responsesFileName, 'fittedSTFs', 'rgcIndicesToAnalyze', '-append');

end

function RcDegs = retinalRFcenterRcDegsInitialEstimate(theMidgetRGCmosaic, theRGCindex)
    connectivityVector = full(squeeze(theMidgetRGCmosaic.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
    indicesOfCenterCones = find(abs(connectivityVector) > 0.0001);
    conesNumPooledByTheRFcenter = numel(indicesOfCenterCones);
    coneRcDegs = mean(theMidgetRGCmosaic.inputConeMosaic.coneApertureDiametersDegs(indicesOfCenterCones)) * ...
                 theMidgetRGCmosaic.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;
    RcDegs = sqrt(conesNumPooledByTheRFcenter)*coneRcDegs;
end

function [theMeasuredSTFtoFit, theMeasuredSTFs] = STFtoFit(theMidgetRGCMosaicResponses, spatialFrequenciesTested, orientationsTested)

    % Allocate memory
    theMeasuredSTFs = zeros(numel(orientationsTested),numel(spatialFrequenciesTested));

    for iSF = 1:numel(spatialFrequenciesTested)
        for iOri = 1:numel(orientationsTested)
            % Retrieve the mRGC response time-series
            theResponseTimeSeries = squeeze(theMidgetRGCMosaicResponses(iOri, iSF, :));
    
            % Compute the response modulation for this SF
            theMeasuredSTFs(iOri, iSF) = max(theResponseTimeSeries)-min(theResponseTimeSeries);
        end % iORI
    end % iSF

    theMeasuredSTFs = theMeasuredSTFs / max(theMeasuredSTFs(:));

    % Determine the orientation that maximizes the STF extension to high spatial frequencies
    maxSF = nan(1,numel(orientationsTested));
    spatialFrequenciesInterpolated = linspace(spatialFrequenciesTested(1),spatialFrequenciesTested(end), 50);

    for iOri = 1:numel(orientationsTested)
        % Find spatial frequency at which STF drops to 20% of max
        theSTFatThisOri = squeeze(theMeasuredSTFs(iOri,:));
        theSTFatThisOriInterpolated = interp1(spatialFrequenciesTested, theSTFatThisOri, spatialFrequenciesInterpolated);
        [mag, iSFpeak] = max(theSTFatThisOri);
        thresholdSTF = mag * 0.2;

        ii = iSFpeak;
        keepGoing = true; iStop = [];
        while (ii < numel(spatialFrequenciesInterpolated)-1)&&(keepGoing)
            ii = ii + 1;
            if (theSTFatThisOriInterpolated(ii)>=thresholdSTF) && (theSTFatThisOriInterpolated(ii+1)<thresholdSTF)
                keepGoing = false;
                iStop = ii;
            end
        end % while
        if (~isempty(iStop))
            maxSF(iOri) = spatialFrequenciesInterpolated(iStop);
        end
    end % iOri

    % Best orientation
    if (any(isnan(maxSF)))
        theSTFatTheHighestSF = squeeze(theMeasuredSTFs(:,end));
        [~, iBestOri] = max(theSTFatTheHighestSF(:));
    else
        [~, iBestOri] = max(maxSF);
    end
    theMeasuredSTFtoFit = squeeze(theMeasuredSTFs(iBestOri,:));
end

function computeMosaicAchromaticSTFs(mosaicFileName, responsesFileName)
    % Load the fully-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    viewingDistanceMeters = 4;
    stimulusPixelsNum = 512*2;
    coneContrasts = [1 1 0];
    deltaOri = 15;
    orientationsTested = 0:deltaOri:(180-deltaOri);
    spatialFrequenciesTested = [0.25 0.5 1 2 4 6 8 12 16 20 24 32 48 64];

    % Generate a presentation display with a desired resolution
    sceneFOVdegs = theMidgetRGCmosaic.inputConeMosaic.sizeDegs;
    retinalImageResolutionDegs = max(sceneFOVdegs)/stimulusPixelsNum;

    % At least 6 samples / period
    maxSF = 1/(2*3*retinalImageResolutionDegs);
    if (max(spatialFrequenciesTested) > maxSF)
        fprintf('Max SF examined (%2.2f c/deg) is too high for this FOV (%2.2f degs) and pixels num (%d). (SFmax: %2.2f c/deg)\n', ...
            max(spatialFrequenciesTested), max(sceneFOVdegs), stimulusPixelsNum, maxSF);
        idx = find(spatialFrequenciesTested <= maxSF);
        spatialFrequenciesTested = spatialFrequenciesTested(idx);
        if (maxSF > max(spatialFrequenciesTested))
            spatialFrequenciesTested(numel(spatialFrequenciesTested)+1) = maxSF;
        end

        fprintf('Will only measure the STF up to %2.2f c/deg.\n', max(spatialFrequenciesTested));
    end

    
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

    % Allocate memory
    stimParams.orientationDegs = 0;
    stimParams.spatialFrequencyCPD = spatialFrequenciesTested(1);
    [~, spatialPhasesDegs] = rfMappingStimulusGenerator.driftingGratingFrames(stimParams);
    rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsMicrons,1);
    theMidgetRGCMosaicResponses = ...
        zeros(numel(orientationsTested), numel(spatialFrequenciesTested), numel(spatialPhasesDegs), rgcsNum);
   
    disp('Allocated memory');
    
    % No optics. We will get it back during the first call to
    % mRGCMosaic.compute()
    theOptics = [];

    % Go through all stimulus orientations
    for iOri = 1:numel(orientationsTested)
        stimParams.orientationDegs = orientationsTested(iOri);

        fprintf('Computing STF for the %d degs orientation patterns.\n', ...
                stimParams.orientationDegs);

        for iFreq = 1:numel(spatialFrequenciesTested)

            theStimParams = stimParams;
            theStimParams.spatialFrequencyCPD = spatialFrequenciesTested(iFreq);
            
            % Generate spatial modulation patterns for each stimulus frame
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

                fprintf('Computing mRGC mosaic response to frame (%d/%d) of the %2.2f c/deg stimulus.\n', ...
                    iFrame, numel(spatialPhasesDegs), theStimParams.spatialFrequencyCPD);

                % Get scene corresponding to this stimulus frame
                theScene = theDriftingGratingFrameScenes{iFrame};

                if (isempty(theOptics))
                    % Compute the mosaic's response to this stimulus frame
                    % and also retrieve the optics so we can pass it along
                    % in subsequent calls (avoid recomputing it)
                    [r,~,theConeMosaicActivation, theOptics] = theMidgetRGCmosaic.compute(...
                        theScene, ...
                        'nTrials', 1, ...
                        'theNullScene', theNullStimulusScene, ...
                        'normalizeConeResponsesWithRespectToNullScene', true);
                else
                    % Compute the mosaic's response to this stimulus frame
                    % using the returned optics
                    [r,~,theConeMosaicActivation] = theMidgetRGCmosaic.compute(...
                        theScene, ...
                        'nTrials', 1, ...
                        'theNullScene', theNullStimulusScene, ...
                        'withOptics', theOptics, ...
                        'normalizeConeResponsesWithRespectToNullScene', true);
                end


                % Store the mosaic's responses
                theFrameResponses(iFrame,:) = r;

            end % iFrame

            % Save memory
            theDriftingGratingFrameScenes = [];
            theNullStimulusScene = [];

            theMidgetRGCMosaicResponses(iOri, iFreq,:,:) = theFrameResponses;

        end % iFreq
    end % iOri

    % Save all response data to disk
    save(responsesFileName, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts', '-v7.3');
end



function visualizeFittedSpatialRFmaps(mosaicFileName, rgcIndicesToAnalyze)
    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    fprintf(['Center/Surround RFs of theMidgetRGCmosaic were generated based on %d RTVF objects\n' ...
             'for %d different cases of center cones num\n'], ...
             numel(theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList), ...
             numel(unique(theMidgetRGCmosaic.theConesNumPooledByTheRFcenterGrid)));

    videoFileName = 'RFs';
    videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    
    % Visualize the PSF of the RTVF object that was closest to a particular mosaic position
    targetMosaicPosition = [0 0];
    [~, objIndexForVisualizingPSF] = min(sum((bsxfun(@minus, theMidgetRGCmosaic.theSamplingPositionGrid, targetMosaicPosition)).^2,2));
    visualizedPSF = 'target positions';

    % Or visualize the PSF of the RTVF object that is closest to each RGC
    visualizedPSF = 'nearest';


    for iRGC = 1:numel(rgcIndicesToAnalyze)
        % Retrieve theTargetRGCindex
        theTargetRGCindex = rgcIndicesToAnalyze(iRGC);

        [triangulatingRTVFobjIndices, triangulatingRTVFobjWeights] = ...
             theMidgetRGCmosaic.triangulatingRTVFobjectIndicesAndWeights(theTargetRGCindex);
        [~,idx] = max(triangulatingRTVFobjWeights);
        nearestRTVFobjectIndex = triangulatingRTVFobjIndices(idx);

        % Extract the PSF data to visualize
        if (strcmp(visualizedPSF,'nearest'))
            % Visualize the PSF of the nearest RTVF object (i.e. max weight for this RGC)
            objIndexForVisualizingPSF = nearestRTVFobjectIndex;
        end
        theRTVFTobj = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{objIndexForVisualizingPSF};
        theVisualizedPSFData = theRTVFTobj.theVlambdaWeightedPSFData;
        theVisualizedPSFData.psfSupportXdegs = theVisualizedPSFData.psfSupportXdegs + ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(theTargetRGCindex,1);
        theVisualizedPSFData.psfSupportYdegs = theVisualizedPSFData.psfSupportYdegs + ...
            theMidgetRGCmosaic.rgcRFpositionsDegs(theTargetRGCindex,2);
       

        % Extract the STF of the nearest RTVFobject
        theNearestRTVFTobj = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{nearestRTVFobjectIndex};

        % Normalize theRTVFTobj.fitted to the theRTVFTobj.targetSTF
        scaleFactorBasedOnLowestSF = theNearestRTVFTobj.rfComputeStruct.theSTF.target(1)/theNearestRTVFTobj.rfComputeStruct.theSTF.fitted(1);
 
        % Assemble the nearest RTVF object STF data
        theNearestRTVFobjSTFdata = struct();
        theNearestRTVFobjSTFdata.spatialFrequencySupport = theNearestRTVFTobj.rfComputeStruct.theSTF.support(:);
        % The model-achieved STF
        theNearestRTVFobjSTFdata.fittedSTF = theNearestRTVFTobj.rfComputeStruct.theSTF.fitted(:)*scaleFactorBasedOnLowestSF;
        % The model-target STF
        theNearestRTVFobjSTFdata.targetSTF = theNearestRTVFTobj.rfComputeStruct.theSTF.target(:);


        [hFig, allAxes] = theMidgetRGCmosaic.visualizeSpatialRFs(...
                'onlyForRGCwithIndex', theTargetRGCindex, ...
                'visualizedRFspatialExtent', 0.2, ...
                'withPSFData', theVisualizedPSFData, ...
                'generateVideo', false, ...
                'withEccentricityCrossHairs', true, ...
                'fontSize', 16);

        % replace graphic is (1,1) with the STFs of the nearest RTVF
        cla(allAxes{1,1});
        plot(allAxes{1,1}, theNearestRTVFobjSTFdata.spatialFrequencySupport, theNearestRTVFobjSTFdata.targetSTF, 'k-', 'LineWidth', 1.5);
        hold(allAxes{1,1}, 'on');
        plot(allAxes{1,1}, theNearestRTVFobjSTFdata.spatialFrequencySupport, theNearestRTVFobjSTFdata.fittedSTF, 'r-', 'LineWidth', 1.5);
        legend(allAxes{1,1},{'RTVF: target', 'RTVF: fitted'}, 'Location', 'SouthWest');
        set(allAxes{1,1}, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], 'XLim', [0.01 100]);
        grid(allAxes{1,1}, 'on')
        xlabel(allAxes{1,1}, 'spatial frequency (c/deg)')
        set(allAxes{1,1}, 'FontSize', 16);

        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
    end % iRGC

    videoOBJ.close();
    fprintf('spatial RFs video saved at %s\n', videoFileName);

end


function generateCenterSurroundRFs(mosaicFileName, R2VFTobjFileName)

    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    % Load the computed R2VFTobjects
    load(R2VFTobjFileName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');

    % Generate the center/surround RFs
    theMidgetRGCmosaic.generateCenterSurroundSpatialPoolingRFs(theRTFVTobjList, ...
                        theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid, ...
                        theVisualSTFSurroundToCenterRcRatioGrid, theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid);

    % Save the updated midgetRGCmosaic which now includes  the computed
    % RTVFTobjList as well as the different grids:
    %  - 'theOpticsPositionGrid'
    %  - 'theConesNumPooledByTheRFcenterGrid'
    %  - 'theVisualSTFSurroundToCenterRcRatioGrid'
    %  - 'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid'
    save(mosaicFileName, 'theMidgetRGCmosaic', '-v7.3');
end


function inspectSavedRTVFfile(fName)
    load(fName, 'obj', 'iRTVobjIndex', ...
        'theConesNumPooledByTheRFcenterGrid', ...
        'theSamplingPositionGrid', ...                         
        'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid', ...                                     
        'theVisualSTFSurroundToCenterRcRatioGrid');

    theRTVFTobj = obj;
    clear 'obj';
    
    % Target and achieved ratios
    targetRsRcRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterRcRatio;
    targetSCintSensRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;
    fittedRsRcRatio = theRTVFTobj.rfComputeStruct.theSTF.fittedRsRcRatio;
    fittedSCintSensRatio = theRTVFTobj.rfComputeStruct.theSTF.fittedSCIntSensRatio;

    fprintf('RTVFobj at position (degs): %2.2f %2.2f\n', theSamplingPositionGrid(iRTVobjIndex,1), theSamplingPositionGrid(iRTVobjIndex,2));
    fprintf('Target Rs/Rc ratio: %2.2f, achieved: %2.2f\n', targetRsRcRatio, fittedRsRcRatio);
    fprintf('Target S/C int. sens. ratio: %2.3f, achieved: %2.3f\n', targetSCintSensRatio, fittedSCintSensRatio);

    hFig = figure(999); clf;
    set(hFig, 'Position', [10 10 900 350], ...
        'Name', sprintf('RTVF at position (degs): %2.2f, %2.2f', ...
                         theSamplingPositionGrid(iRTVobjIndex,1), ...
                         theSamplingPositionGrid(iRTVobjIndex,2)));
    ax = subplot(1,2,1);
    XLims = [1 15]; YLims = [1 15];
    plotTargetAndAchievedParam(ax, targetRsRcRatio, fittedRsRcRatio, XLims, YLims, 'Rs/Rc ratio');

    ax = subplot(1,2,2);
    XLims = [0 1]; YLims = [0 2];
    plotTargetAndAchievedParam(ax, targetSCintSensRatio, fittedSCintSensRatio, XLims, YLims, 'int. sens. S/C ratio');

end

function plotTargetAndAchievedParam(ax, targetVal, fittedVal, XLims, YLims, titleString)
    
    plot(ax, XLims, YLims, 'k-', 'LineWidth', 1.0); hold(ax,'on');
    plot(ax, targetVal, fittedVal, 'ro', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.5 0.5]);
    axis(ax, 'square')
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 16);
    title(ax, titleString);
    xlabel(ax,'target');
    ylabel(ax,'achieved');
end


function generateR2VFTobjects(mosaicCenterParams, mosaicSurroundParams, opticsParams, mosaicFileName, mosaicDirectory, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('updateRTVFobjectAtPosition', [], @(x)(isempty(x) || (numel(x)==2)));
    p.addParameter('updateRTVFobjectWithCenterConesNum', [], @(x)(isempty(x) || (numel(x)==1)));
    p.parse(varargin{:});

    updateRTVFobjectAtPosition = p.Results.updateRTVFobjectAtPosition;
    updateRTVFobjectWithCenterConesNum = p.Results.updateRTVFobjectWithCenterConesNum;

    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');
    
    % Generate the eccentricitySamplingGrid
    eccentricitySamplingGrid = midgetRGCMosaic.eccentricitySamplingGridCoords(...
        mosaicCenterParams.positionDegs, mosaicCenterParams.sizeDegs, ...
        mosaicSurroundParams.eccentricitySamplingGridHalfSamplesNum, ...
        'hexagonal', true);

    if (~isempty(updateRTVFobjectAtPosition))
        targetPosition = [updateRTVFobjectAtPosition(1) updateRTVFobjectAtPosition(2)];
        d = sum((bsxfun(@minus, eccentricitySamplingGrid, targetPosition)).^2,2);
        [~,updatedRTVFobjectIndex] = min(d);
        eccentricitySamplingGrid = eccentricitySamplingGrid(updatedRTVFobjectIndex,:);
        fprintf(2, 'Will update the RTVFobjList at index %d (pos (degs) = %2.2f,%2.2f)', ...
            updatedRTVFobjectIndex, eccentricitySamplingGrid(1), eccentricitySamplingGrid(2));
    end


    % Visualize the mosaic
    theMidgetRGCmosaic.visualize( ...
        'eccentricitySamplingGrid', eccentricitySamplingGrid, ...
        'inputPoolingVisualization', 'centerOnly');

    % multiStartsNum: select from:
    % - 1 (Single start run, fastest results), 
    % - some number (Multi-start), or 
    % - inf (Global search)
    fitParams = struct();
    fitParams.multiStartsNumDoGFit = 128;

    % Where to save the fitted RTVFobjects
    fitParams.exportsDirectory = mosaicDirectory;

    tStart = cputime;

    % Generate list of RTVT objects
    [theRTFVTobjList, theOpticsPositionGrid, ...
     theConesNumPooledByTheRFcenterGrid, ...
     theVisualSTFSurroundToCenterRcRatioGrid, ...
     theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid] = midgetRGCMosaic.R2VFTobjects(...
                theMidgetRGCmosaic, eccentricitySamplingGrid, ...
                mosaicSurroundParams, opticsParams, fitParams);

    timeLapsedMinutes = (cputime - tStart)/60;
    fprintf('\n\n midgetRGCMosaic.R2VFTobjects were generated in %d positions and fitting took %f minutes\n', ...
        size(eccentricitySamplingGrid,1), timeLapsedMinutes);
    
    % Save the computed list of RTVFTobj and the various grids to the mosaic mat file
    R2VFTobjFileName = R2VFTobjFileNameForMosaicFileName(mosaicFileName, mosaicSurroundParams.H1cellIndex);

    if (~isempty(updateRTVFobjectAtPosition))

        % Compute sourceRTVFobjectIndex
        sourceRTVFobjectIndex = find(...
            (theConesNumPooledByTheRFcenterGrid == updateRTVFobjectWithCenterConesNum));

        tmp = theRTFVTobjList{sourceRTVFobjectIndex};
        clear 'theRTFVTobjList';
        
        % Load previously generated theRTFVTobjList
        load(R2VFTobjFileName, 'theRTFVTobjList', ...
            'theOpticsPositionGrid', 'theConesNumPooledByTheRFcenterGrid');
        
        % Compute destinationRTVFobjectIndex
        % Retrieve the indices of the fitted RTVF objects that have the
        % same # of center cones
        centerConeMatchObjIndices = find(theConesNumPooledByTheRFcenterGrid == updateRTVFobjectWithCenterConesNum);

        % Compute distance based weights for this RGC and the fitted RTVF objects
        distancesToSamplingGridPositions = sqrt(sum((bsxfun(@minus, theOpticsPositionGrid(centerConeMatchObjIndices,:), targetPosition)).^2,2));
        [~, idx] = sort(distancesToSamplingGridPositions, 'ascend');
        destinationRTVFobjectIndex = idx(1);

        fprintf('Will override the RTVF with %d cones at position (degs): %2.2f %2.2f\n', ...
            theConesNumPooledByTheRFcenterGrid(destinationRTVFobjectIndex), ...
            theOpticsPositionGrid(destinationRTVFobjectIndex,1), theOpticsPositionGrid(destinationRTVFobjectIndex,2));
        pause

        % Overwrite list at the updatedRTVFobjectIndex
        theRTFVTobjList{destinationRTVFobjectIndex} = tmp;

        % Save the updated RTVFT list
        save(R2VFTobjFileName, 'theRTFVTobjList', '-append');
    else

        save(R2VFTobjFileName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                    'theConesNumPooledByTheRFcenterGrid', ...
                    'theVisualSTFSurroundToCenterRcRatioGrid', ...
                    'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');
        fprintf('Computed R2VFTobjects saved in: %s\n', R2VFTobjFileName);
    end


end

function generateCenterConnectedMosaic(mosaicParams, mosaicFileName)
    % Generate mRGC mosaic
    theMidgetRGCmosaic = midgetRGCMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'whichEye', mosaicParams.whichEye, ...
                'eccentricityDegs', mosaicParams.positionDegs, ...
                'sizeDegs', mosaicParams.sizeDegs ...
                );
    % Save the center-connected mosaic
    save(mosaicFileName, 'theMidgetRGCmosaic', '-v7.3');
end

function  responsesFileName = responsesFileNameForMosaicFileName(mosaicFileName, H1cellIndex)
    responsesFileName = strrep(mosaicFileName, '.mat', sprintf('_H1cellIndex%d_Responses.mat', H1cellIndex));
end


function R2VFTobjFileName = R2VFTobjFileNameForMosaicFileName(mosaicFileName, H1cellIndex)
    R2VFTobjFileName = strrep(mosaicFileName, '.mat', sprintf('_R2VFTobjectsForH1cellIndex_%d.mat', H1cellIndex));
end

function [mosaicFileName, mosaicDirectoryPath] = generateMosaicFileName(mosaicParams)
    dropboxDir = localDropboxPath();
    mosaicDirectoryPath = sprintf('%s/productionMidgetRGCMosaics', dropboxDir);
    mosaicFileName = sprintf('MRGCmosaic_Ecc_%2.1f_%2.1f_sizeDegs_%2.1f_%2.1f.mat', ...
        mosaicParams.positionDegs(1), mosaicParams.positionDegs(2), ...
        mosaicParams.sizeDegs(1), mosaicParams.sizeDegs(2));

    mosaicFileName = fullfile(mosaicDirectoryPath, mosaicFileName);
end


function dropboxDir = localDropboxPath()
    % Get dropboxDir & intermediate data files location
    computerInfo = GetComputerInfo();
    switch (computerInfo.localHostName)
        case 'Ithaka'
            dropboxDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

        case 'Crete'
            dropboxDir = '/Volumes/Dropbox/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

        case 'Santorini'
            dropboxDir = '/Users/nicolas/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

        otherwise
            if (contains(computerInfo.networkName, 'leviathan'))
                dropboxDir = '/media/dropbox_disk/Aguirre-Brainard Lab Dropbox/isetbio isetbio';
            else
                error('Could not establish dropbox location')
            end
    end
end
