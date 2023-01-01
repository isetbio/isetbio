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
    actionToPerform = 'generateSurroundParameterEccentricityMaps';
    actionToPerform = 'computeMosaicSTFs';

    switch (actionToPerform)

        case 'computeMosaicSTFs'
            responsesFileName = responsesFileNameForMosaicFileName(mosaicFileName, H1cellIndex);
            computeMosaicAchromaticSTFs(mosaicFileName, responsesFileName);


        case 'generateSurroundParameterEccentricityMaps'
            R2VFTobjFileName = R2VFTobjFileNameForMosaicFileName(mosaicFileName, H1cellIndex);
            generateSurroundParameterEccentricityMaps(mosaicFileName, R2VFTobjFileName);

        case 'generateR2VFTobjects' 
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

            generateR2VFTobjects(mosaicCenterParams, mosaicSurroundParams, opticsParams, mosaicFileName, mosaicDirectory);

        case 'generateCenterConnectedMosaic'
            generateCenterConnectedMosaic(mosaicCenterParams, mosaicFileName);

        otherwise
            error('Unknown action: ''%s''.', actionToPerform);
    end
end


function computeMosaicAchromaticSTFs(mosaicFileName, responsesFileName)
    % Load the fully-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    viewingDistanceMeters = 4;
    stimulusPixelsNum = 512;
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
   
    disp('Allocated memory')
    

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

                % Compute the mosaic's response to this stimulus frame
                [r,~,theConeMosaicActivation] = theMidgetRGCmosaic.compute(...
                    theScene, ...
                    'nTrials', 1, ...
                    'theNullScene', theNullStimulusScene, ...
                    'normalizeConeResponsesWithRespectToNullScene', true);

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


function generateSurroundParameterEccentricityMaps(mosaicFileName, R2VFTobjFileName)

    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    % Load the computed R2VFTobjects
    load(R2VFTobjFileName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');


    visualizeFittedModels = ~true;

    
    if (visualizeFittedModels)
        theConesNumPooledList = unique(theConesNumPooledByTheRFcenterGrid);
        
        xLims(1) = min(theRTFVTobjList{1}.theConeMosaic.coneRFpositionsDegs(:,1));
        xLims(2) = max(theRTFVTobjList{1}.theConeMosaic.coneRFpositionsDegs(:,1));
        yLims(1) = min(theRTFVTobjList{1}.theConeMosaic.coneRFpositionsDegs(:,2));
        yLims(2) = max(theRTFVTobjList{1}.theConeMosaic.coneRFpositionsDegs(:,2));

        for theConesNumPooledIndex = 1:numel(theConesNumPooledList)
            iRTVobjIndices = find(theConesNumPooledByTheRFcenterGrid == theConesNumPooledList(theConesNumPooledIndex));
            spatialPositions = theOpticsPositionGrid(iRTVobjIndices,:);
            
            figure(theConesNumPooledIndex); clf;
            ax = subplot('Position', [0.01 0.01 0.99 0.99]);
            set(ax, 'Color', [0 0 0]);

            for idx = 1:numel(iRTVobjIndices)
                iRTVobjIndex = iRTVobjIndices(idx)
                theRTVFobj = theRTFVTobjList{iRTVobjIndex};
    
                [~, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, ...
                    pooledConeIndicesAndWeights] = RetinaToVisualFieldTransformer.visualRFfromRetinalConePooling(...
                     theRTVFobj.rfComputeStruct.modelConstants, ...
                     theRTVFobj.rfComputeStruct.retinalConePoolingParams.finalValues);
    
                spatialPositions(idx,:)
                xDegs = theRTVFobj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1) + spatialPositions(idx,1);
                yDegs = theRTVFobj.rfComputeStruct.modelConstants.spatialSupportDegs(:,2) + spatialPositions(idx,2);
                imagesc(ax, xDegs, yDegs, theRetinalRFsurroundConeMap/max(theRetinalRFsurroundConeMap(:)));
                hold(ax, 'on')
                %plot(ax, spatialPositions(idx,1), spatialPositions(idx,2), 'r+');
                axis 'equal'
                colormap(ax, gray(1024))
                

            end
    
            pause
            
            %title('spatial positions fitted (%d center cones)', theConesNumPooledList(theConesNumPooledIndex)); 
            %axis 'equal'
            %set(gca, 'XLim', xLims, 'yLim', yLims);
        end
    end

    theMidgetRGCmosaic.generateCenterSurroundSpatialPoolingRF(theRTFVTobjList, ...
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


function generateR2VFTobjects(mosaicCenterParams, mosaicSurroundParams, opticsParams, mosaicFileName, mosaicDirectory)

    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    
    % Generate the eccentricitySamplingGrid
    eccentricitySamplingGrid = midgetRGCMosaic.eccentricitySamplingGridCoords(...
        mosaicCenterParams.positionDegs, mosaicCenterParams.sizeDegs, ...
        mosaicSurroundParams.eccentricitySamplingGridHalfSamplesNum, ...
        'hexagonal', true);


    %centerConesNum = full(sum(theMidgetRGCmosaic.rgcRFcenterConeConnectivityMatrix,1))
    %[max(centerConesNum) min(centerConesNum) numel(find(centerConesNum ==max(centerConesNum))) numel(find(centerConesNum ==min(centerConesNum)))]

    % Visualize the mosaic
    theMidgetRGCmosaic.visualize( ...
        'eccentricitySamplingGrid', eccentricitySamplingGrid, ...
        'inputPoolingVisualization', 'centerOnly');


    % multiStartsNum: select from:
    % - 1 (Single start run, fastest results), 
    % - some number (Multi-start), or 
    % - inf (Global search)
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
    fprintf('\n\n midgetRGCMosaic.R2VFTobjects were generated in %d positions and fitting took %f minutes\n', timeLapsedMinutes);
    

    % Save the computed list of RTVFTobj and the various grids to the mosaic mat file
    R2VFTobjFileName = R2VFTobjFileNameForMosaicFileName(mosaicFileName, mosaicSurroundParams.H1cellIndex);

    save(R2VFTobjFileName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');
    fprintf('Computed R2VFTobjects saved in: %s\n', R2VFTobjFileName);

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
