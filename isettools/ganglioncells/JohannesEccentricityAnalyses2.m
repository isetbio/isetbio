function JohannesEccentricityAnalyses2

    % Get dropboxDir & intermediate data files location
    computerInfo = GetComputerInfo();
    switch (computerInfo.localHostName)
        case 'Ithaka'
            dropboxDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris/midgetRGCMosaics';
            mappedRFsDir = '/Volumes/SSDdisk/MATLAB/toolboxes/isetbio/isettools/ganglioncells/JohannesAnalysesData';
   
        case 'Crete'
            dropboxDir = '/Volumes/Dropbox/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris/midgetRGCMosaics';
            mappedRFsDir = '/Volumes/MATLAB/toolboxes/isetbio/isettools/ganglioncells/JohannesAnalysesData';

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

    % 3 degrees in the temporal retina (negative horizontal ecc)
    mosaicEccDegs = [-3 0];
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

    % multiStartsNum: select from:
    % - 1 (Single start run, fastest results), 
    % - some number (Multi-start), or 
    % - inf (Global search)
    multiStartsNum = 3;

    generateRTVobjects = ~true;
    if (generateRTVobjects)
        for iEcc = 1:size(mosaicEccDegs,1)
            tic
            fprintf('Generating RTVF objects for mosaic %d of %d ... \n', iEcc, size(mosaicEccDegs,1));
            generateTheRTVFobjects(mosaicEccDegs(iEcc,:), mappedRFsDir, ...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM, multiStartsNum);
            fprintf('Finished generating RTVF objects for mosaic %d of %d in %2.1f hours.\n', iEcc, size(mosaicEccDegs,1), toc/(60*60));
        end
    end

    generateCenterSurroundRFstructure = true;
    if (generateCenterSurroundRFstructure)
        for iEcc = 1:size(mosaicEccDegs,1)
            generateTheCenterSurroundRFs(mosaicEccDegs(iEcc,:), mappedRFsDir, ...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM);
        end
    end

    % STF stimulus parameters
    % L+M contrast for gratings used to measure the STFs 
    coneContrasts = [1 1 0];
    stimulusPixelsNum = 300;
    deltaOri = 15;
    orientationsTested = 0:deltaOri:(180-deltaOri);
    spatialFrequenciesTested = [0.25 0.5 1 2 4 6 8 12 16 20 24 32 48 64];

%     stimulusPixelsNum = 256;
%     orientationsTested = 0:30:150;
%     spatialFrequenciesTested = [0.25 0.5 1 2 4 8 16 32 64];


    computeTheSTFs = ~true;
    if (computeTheSTFs)
        for iEcc = 1:size(mosaicEccDegs,1)
            computeTheMosaicSTFs(mosaicEccDegs(iEcc,:), mappedRFsDir, ...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                coneContrasts, stimulusPixelsNum, ...
                orientationsTested, spatialFrequenciesTested);
        end
    end


    fitTheSTFs = ~true;
    centerMostRGCsNumToAnalyze = 16;
    if (fitTheSTFs)
        for iEcc = 1:size(mosaicEccDegs,1)
            fitTheMosaicSTFs(mosaicEccDegs(iEcc,:), mappedRFsDir, ...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                coneContrasts, centerMostRGCsNumToAnalyze);
        end
    end

    inspectSyntheticRFsComponents = true;
    if (inspectSyntheticRFsComponents)
        for iEcc = 1:size(mosaicEccDegs,1)
            inspectSyntheticRFs(mosaicEccDegs(iEcc,:), mappedRFsDir, ...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                coneContrasts, centerMostRGCsNumToAnalyze);
        end
    end
end

function inspectSyntheticRFs(mosaicEccDegs, mappedRFsDir, ZernikeDataBase, subjectRankOrder, ...
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
    
    for iRGC = 1:numel(fittedSTFs)
        f = fittedSTFs{iRGC}
        conesNumPooledByTheRFcenter = f.targetVisualRFDoGparams.conesNumPooledByTheRFcenter;
        theRGCindex = f.targetRGC;
        theFittedSTF = f.theFittedSTF;
        theTargetSTFsupport = f.theTargetSTFdata(:,1);
        theTargetSTFDoGModelApproximation = f.theTargetSTFdata(:,2);
        theTargetSTFmeasured = f.theTargetSTFdata(:,3);

        theFittedSTFDoGparams = f.theFittedSTFDoGparams
        figure(1); clf;
        p1 = plot(spatialFrequenciesTested, f.allMeasuredSTFs, '.', ...
            'MarkerSize', 20, 'LineWidth', 1.0, 'MarkerEdgeColor', [0.5 0.5 0.5]);
        hold on;
        p2 = plot(spatialFrequenciesTested, f.theMeasuredSTFtoFit, 'ro', 'MarkerSize', 14, 'LineWidth', 1.5);
        p3 = plot(f.theFittedSTF.sfHiRes, f.theFittedSTF.compositeSTFHiRes, 'r-', 'LineWidth', 1.5);
        p4 = plot(theTargetSTFsupport, theTargetSTFDoGModelApproximation, 'b-', 'LineWidth', 1.5);
        p5 = plot(theTargetSTFsupport, theTargetSTFmeasured, 'c-', 'LineWidth', 1.5);
        legend([p1(1) p2 p3 p4 p5], ...
               {'measured (all orientations)','measured (opt. orientation)', 'measured (DoG model fit)', 'target (DoGmodel fit)', 'target (measured)'}, ...
               'Location', 'SouthWest');
        title(sprintf('RF center: %d input cones', conesNumPooledByTheRFcenter));
        set(gca, 'XScale', 'log');
        pause
    end


end


function fitTheMosaicSTFs(mosaicEccDegs, mappedRFsDir, ZernikeDataBase, subjectRankOrder, ...
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
         'spatialPhasesDegs');


    % Find the indices of the centerMostRGCsNumToAnalyze RGCs
    relativeRGCpositions = bsxfun(@minus, theMidgetRGCmosaic.rgcRFpositionsDegs, theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs);
    radii = sum(relativeRGCpositions.^2,2);
    [~,sortedRGCindices] = sort(radii, 'ascend');
    rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsDegs,1);
    rgcIndicesOfAnalyzedRFs = sortedRGCindices(1:min([centerMostRGCsNumToAnalyze, rgcsNum]));

    % Allocate memory
    fittedSTFs = cell(1, numel(rgcIndicesOfAnalyzedRFs));

    for iRGC = 1:numel(rgcIndicesOfAnalyzedRFs)

        % Retrieve the RGCindex
        theRGCindex = rgcIndicesOfAnalyzedRFs(iRGC);

        % Compute the STFs for all examined orientations
        theMeasuredSTFs = zeros(numel(orientationsTested),numel(spatialFrequenciesTested));
        
        for iSF = 1:numel(spatialFrequenciesTested)
            for iOri = 1:numel(orientationsTested)
                % Retrieve the mRGC response time-series
                theResponseTimeSeries = squeeze(theMidgetRGCMosaicResponses(iOri, iSF, :, theRGCindex));
    
                % Compute the response modulation for this SF
                theMeasuredSTFs(iOri, iSF) = max(theResponseTimeSeries)-min(theResponseTimeSeries);
            end % iORI
        end % iSF
        theMeasuredSTFs = theMeasuredSTFs / max(theMeasuredSTFs(:));
    
        % Determine the orientation that maximizes the STF extension to high spatial frequencies
        maxSF = nan(1,numel(orientationsTested));
        for iOri = 1:numel(orientationsTested)
            % Find spatial frequency at which STF drops to 20% of max
            theSTFatThisOri = squeeze(theMeasuredSTFs(iOri,:));

            spatialFrequenciesInterpolated = linspace(spatialFrequenciesTested(1),spatialFrequenciesTested(end), 50);
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

        

        % Initial estimate of the retinal RF center
        connectivityVector = full(squeeze(theMidgetRGCmosaic.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
        indicesOfCenterCones = find(abs(connectivityVector) > 0.0001);
        conesNumPooledByTheRFcenter = numel(indicesOfCenterCones);
        coneRcDegs = mean(theMidgetRGCmosaic.inputConeMosaic.coneApertureDiametersDegs(indicesOfCenterCones)) * ...
                     theMidgetRGCmosaic.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;
        retinalRFcenterRcDegsInitialParam = sqrt(conesNumPooledByTheRFcenter)*coneRcDegs;

        % Fit the DoG model to the measured STF
        multiStartsNum = 64;
        [DoGparams, theFittedSTF] = ...
            fitDoGmodelToMeasuredSTF(spatialFrequenciesTested, ...
                    theMeasuredSTFtoFit, ...
                    retinalRFcenterRcDegsInitialParam, ...
                    multiStartsNum);

        % Save the fit results

        % Retrieve the correct RTVFTobj based on this cells position and
        % #of center cones. For now only checking the centerConesNum
        iObj = find(...
             (theMidgetRGCmosaic.theConesNumPooledByTheRFcenterGrid == numel(indicesOfCenterCones)) ...  % match the conesNum in the center
        );
        theRTVFTobj = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList{iObj};

        % Extract the targetSTF data
        maxTargetSTF = max([max(theRTVFTobj.rfComputeStruct.theSTF.fitted(:)) max(theRTVFTobj.rfComputeStruct.theSTF.target(:))]);
        theTargetSTFdata = [ ...
            theRTVFTobj.rfComputeStruct.theSTF.support(:) ...
            theRTVFTobj.rfComputeStruct.theSTF.fitted(:)/maxTargetSTF ...
            theRTVFTobj.rfComputeStruct.theSTF.target(:)/maxTargetSTF];

        fittedSTFs{iRGC} = struct(...
               'targetRGC', theRGCindex, ...
               'targetRGCeccentricityDegs', theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,:), ...
               'targetVisualRFDoGparams', theRTVFTobj.targetVisualRFDoGparams, ...
               'spatialFrequenciesTested', spatialFrequenciesTested, ...
               'allMeasuredSTFs', theMeasuredSTFs, ...
               'theMeasuredSTFtoFit',theMeasuredSTFtoFit, ...
               'theFittedSTF', theFittedSTF, ...
               'theTargetSTFdata', theTargetSTFdata, ...
               'theFittedSTFDoGparams', DoGparams ...
               );

        figure(iRGC); clf;
        plot(spatialFrequenciesTested, theMeasuredSTFtoFit, 'ko'); hold on
        plot(theFittedSTF.sfHiRes, theFittedSTF.compositeSTFHiRes, 'r-');
        plot(theTargetSTFdata(:,1), theTargetSTFdata(:,2), 'b--');
        plot(theTargetSTFdata(:,1), theTargetSTFdata(:,3), 'm--');
        set(gca, 'XScale', 'log', 'XLim', [0.1 100], 'XTick', [0.1 0.3 1 3 10 30 100], 'YLim', [0 2]);
        drawnow;
        
    end % for iRGC

    % Append the fittedSTFs structs
    save(fNameResponses, 'fittedSTFs', '-append');
    fprintf('Appended the measured and fitted STF data for the center-most %d RGCs to %s\n', numel(rgcIndicesOfAnalyzedRFs), fName);

end % fitTheSTFs


function computeTheMosaicSTFs(mosaicEccDegs, mappedRFsDir, ZernikeDataBase, subjectRankOrder, ...
    pupilDiameterMM, coneContrasts, stimulusPixelsNum, ...
    orientationsTested, spatialFrequenciesTested)

    % Load the midgetRGCMosaic
    fName = fullfile(mappedRFsDir, sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2)));
    load(fName, 'theMidgetRGCmosaic');


    % Assemble the responses filename
    responsesFileNamePostFix = sprintf('_STFresponses_cLMS_%2.2f_%2.2f_%2.2f.mat', coneContrasts(1), coneContrasts(2), coneContrasts(3));
    fNameResponses = strrep(fName, '.mat', responsesFileNamePostFix);
    fprintf('STF responses will be saved to %s\n', fNameResponses);

   
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

    % Allocate memory
    stimParams.orientationDegs = 0;
    stimParams.spatialFrequencyCPD = spatialFrequenciesTested(1);
    [~, spatialPhasesDegs] = rfMappingStimulusGenerator.driftingGratingFrames(stimParams);
    rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsMicrons,1);
    theMidgetRGCMosaicResponses = ...
        zeros(numel(orientationsTested), numel(spatialFrequenciesTested), numel(spatialPhasesDegs), rgcsNum);
   
    % Go through all stimulus orientations
    for iOri = 1:numel(orientationsTested)
        stimParams.orientationDegs = orientationsTested(iOri);

        fprintf('Computing STF for the %d degs orientation patterns.\n', ...
                stimParams.orientationDegs);
    
        parfor iFreq = 1:numel(spatialFrequenciesTested)
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

    % Save all response data to disk
    save(fNameResponses, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts', '-v7.3');

end

function generateTheCenterSurroundRFs(mosaicEccDegs, mappedRFsDir, ZernikeDataBase, subjectRankOrder, pupilDiameterMM)
    fName = fullfile(mappedRFsDir, sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2)));
    
    load(fName, 'theMidgetRGCmosaic', 'thePSFData', 'theRTFVTobjList', ...
                'theOpticsPositionGrid', ...
                'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');

    % If we re-run this step a second time, the
    % 'theRTFVTobjList' variable does not exist in the file (it has been
    % embedded in the  midgetRGCMosaic object, so extract it from there)
    if (~exist('theRTFVTobjList', 'var'))
        if (~isempty(theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList))
            fprintf('Loaded theRTFVTobjList from the mosaic itself\n');
            theRTFVTobjList = theMidgetRGCmosaic.theRetinaToVisualFieldTransformerOBJList;
        else
            error('theRTFVTobjList does not exist neither in the file (''%s''), nor in the midgetRGCMosaic object ', fName);
        end
    else
        fprintf('Loaded theRTFVTobjList from the saved file\n');
    end
    

    % Generate C/S spatial RFs for all cells in the  midgetRGCmosaic
    theMidgetRGCmosaic.generateCenterSurroundSpatialPoolingRF(theRTFVTobjList, ...
                        theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid, ...
                        theVisualSTFSurroundToCenterRcRatioGrid, theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid);

    % theRTFVTobjList is now embedded in theMidgetRGCmosaic, so clear it
    clear 'theRTFVTobjList';

    % Save the updated midgetRGCmosaic which now includes  the computed RTVFTobjList
    save(fName, 'theMidgetRGCmosaic', 'thePSFData', ...
                'theOpticsPositionGrid', ...
                'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid', ...
                '-v7.3');


end

function generateTheRTVFobjects(mosaicEccDegs, mappedRFsDir, ZernikeDataBase, subjectRankOrder, pupilDiameterMM, multiStartsNum)
   
    fName = fullfile(mappedRFsDir, sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2)));
    load(fName, 'theMidgetRGCmosaic');

    % Eccentricity grid for RTVTobjects: just a single optical position
    % within the mosaic. For large mosaics, (e.g. 4x4) we should have an array of
    % positions
    eccXY = mean(theMidgetRGCmosaic.rgcRFpositionsDegs,1);
    eccentricitySamplingGrid(1,:) = [eccXY(1) eccXY(2)];

    % Start timer
    tic

    % Generate list of RTVT objects
    [theRTFVTobjList, theOpticsPositionGrid, ...
     theConesNumPooledByTheRFcenterGrid, ...
     theVisualSTFSurroundToCenterRcRatioGrid, ...
     theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid] = midgetRGCMosaic.generateRTVFobjects(...
                ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                theMidgetRGCmosaic, eccentricitySamplingGrid, ...
                multiStartsNum);

    % Measure lapsed time
    timeLapsed = toc;
    
    fprintf('Finished with RTV generation in %f minutes\n', timeLapsed/60);

    % Append the computed list of RTVFTobj and the various grids to the mosaic mat file
    save(fName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid', ...
                '-append');
    fprintf('Appended the computed RTVFTobjList to %s\n', fName);

end


function [DoGparams, theFittedSTF] = fitDoGmodelToMeasuredSTF(sf, theMeasuredSTF, retinalRFcenterRcDegs, multiStartsNum)
    % DoG param initial values and limits: center gain, kc
    Kc = struct(...    
        'low', 1e-1, ...
        'high', 1e4, ...
        'initial', 1);

    % DoG param initial values and limits: Ks/Kc ratio
    KsToKc = struct(...
        'low', 1e-6, ...
        'high', 1.0, ...
        'initial', 0.1);

    % DoG param initial values and limits: RsToRc ratio
    RsToRc = struct(...
        'low', 1.5, ...
        'high', 30, ...
        'initial', 5);

    % DoG param initial values and limits: RcDegs
    RcDegs = struct(...
        'low', retinalRFcenterRcDegs/sqrt(2.0), ...
        'high', retinalRFcenterRcDegs*200, ...
        'initial', retinalRFcenterRcDegs*5);
    
     %                          Kc           kS/kC             RsToRc            RcDegs    
     DoGparams.initialValues = [Kc.initial   KsToKc.initial    RsToRc.initial    RcDegs.initial];
     DoGparams.lowerBounds   = [Kc.low       KsToKc.low        RsToRc.low        RcDegs.low];
     DoGparams.upperBounds   = [Kc.high      KsToKc.high       RsToRc.high       RcDegs.high];
     DoGparams.names         = {'Kc',        'kS/kC',         'RsToRc',         'RcDegs'};
     DoGparams.scaling       = {'log',       'log',           'linear',         'linear'};
     
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
     DoGparams.finalValues = run(ms, problem, multiStartsNum);

     theFittedSTF.compositeSTF = DoGSTF(DoGparams.finalValues, sf);
     theFittedSTF.centerSTF = DoGparams.finalValues(1) * ( pi * DoGparams.finalValues(4)^2 * exp(-(pi*DoGparams.finalValues(4)*sf).^2) );
     theFittedSTF.surroundSTF = DoGparams.finalValues(1)*DoGparams.finalValues(2) * ( pi * (DoGparams.finalValues(4)*DoGparams.finalValues(3))^2 * exp(-(pi*DoGparams.finalValues(4)*DoGparams.finalValues(3)*sf).^2) );
     
     sfHiRes = logspace(log10(0.1), log10(100), 64);
     theFittedSTF.sfHiRes = sfHiRes;
     theFittedSTF.compositeSTFHiRes = DoGSTF(DoGparams.finalValues, sfHiRes);
     theFittedSTF.centerSTFHiRes = DoGparams.finalValues(1) * ( pi * DoGparams.finalValues(4)^2 * exp(-(pi*DoGparams.finalValues(4)*sfHiRes).^2) );
     theFittedSTF.surroundSTFHiRes = DoGparams.finalValues(1)*DoGparams.finalValues(2) * ( pi * (DoGparams.finalValues(4)*DoGparams.finalValues(3))^2 * exp(-(pi*DoGparams.finalValues(4)*DoGparams.finalValues(3)*sfHiRes).^2) );
     
end
