function generateMidgetRGCmosaicComponents
    
    % Operations
    operations = {...
        'generateMosaic' ...
        'generateRTVFTobjList' ...
        'generateCenterSurroundRFs' ...
        'computeConeMosaicSTF' ...
        'visualizeResponses' ...
        };

    % Operation to compute
    operations = operations(1:3);

    eccSizeDegsExamined = [...
         0  0.5; ...
        -1  0.6; ...
        -2  0.75; ...
        -3  0.8; ...
        -4  1.0; ...
        -6  1.25; ...
        -8  1.5; ...
        -10 2.0; ...
        -12 2.5; ...
        -14 3.0; ...
        -16 3.5; ...
        -20 4.0; ...
        -24 5.0 ...
       ];

    for iEcc = 1:numel(eccSizeDegsExamined)
        eccDegs  = eccSizeDegsExamined(iEcc,1) * [1 0];
        sizeDegs = eccSizeDegsExamined(iEcc,2) * [1 1];
        doIt(operations, eccDegs, sizeDegs);
    end

end

function doIt(operations, eccDegs, sizeDegs)

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

                % Generate the list of RTVFT objects for different gridPositions
                iGridPosition = 1;

                eccDegsGrid(iGridPosition,:) = [0 0];
                conesNumPooledByTheRFcenterGrid(iGridPosition) = 1;

                % From Croner & Kaplan '95 (Figure 4c and text)
                % "P surrounds were on average 6.7 times wider than the centers of
                % the same cells, or about 45 times larger in area".
                surroundToCenterRcRatioGrid(iGridPosition) = 6.7;

                % Also from Croner & Kaplan '95 (Figure 10b)
                % "These mean ratios for P and M cells are not significantly different
                % (Student's t-test: P = 0.482). The overall mean ratio is 0.55.
                surroundToCenterIntegratedSensitivityRatioGrid(iGridPosition) = 0.55;
                

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

            case 'computeConeMosaicSTF'
                load(fName, 'theMidgetRGCmosaic');
                [theMidgetRGCMosaicResponses, spatialFrequenciesTested, spatialPhasesDegs] = ...
                    computeTheSTF(theMidgetRGCmosaic);

                % Save the responses to a separate file
                fNameResponses = strrep(fName, '.mat', '_Responses.mat');
                save(fNameResponses, ...
                    'theMidgetRGCMosaicResponses', ...
                    'spatialFrequenciesTested', ...
                    'spatialPhasesDegs', ...
                    '-v7.3');

            case 'visualizeResponses'
                % Load the responses to a separate file
                fNameResponses = strrep(fName, '.mat', '_Responses.mat');
                load(fNameResponses, ...
                    'theMidgetRGCMosaicResponses', ...
                    'spatialFrequenciesTested', ...
                    'spatialPhasesDegs');

                figure(55);
                m1 = min(theMidgetRGCMosaicResponses(:));
                m2 = max(theMidgetRGCMosaicResponses(:));
                for iFreq = 1:numel(spatialFrequenciesTested)
                    m = squeeze(theMidgetRGCMosaicResponses(iFreq,:,:));
                    mm = (m'-m1)/(m2-m1);
                    imagesc(mm);
                    set(gca, 'CLim', [0 1]);
                    title(sprintf('max = %2.3f', max(mm(:))));
                    colormap(gray(1024));
                    pause
                end

            otherwise
                error('Unknown operation: ''%s''.', operations{iOp});
        end % Switch
    end
end



function [theMidgetRGCMosaicResponses, spatialFrequenciesTested, spatialPhasesDegs] = computeTheSTF(theMidgetRGCmosaic)

    sceneFOVdegs = theMidgetRGCmosaic.inputConeMosaic.sizeDegs;

    % Generate a presentation display with a desired resolution
    pixelsNum = 256;
    retinalImageResolutionDegs = max(sceneFOVdegs)/pixelsNum;
    viewingDistanceMeters = 4;
    theDisplay = rfMappingStimulusGenerator.presentationDisplay(...
            theMidgetRGCmosaic.inputConeMosaic.wave, retinalImageResolutionDegs, ...
            viewingDistanceMeters);

    % Stim params for the STF mapping
    stimParams = struct(...
            'backgroundLuminanceCdM2', 50.0, ...
            'backgroundChromaticity', [0.301 0.301], ...
            'coneContrasts', [1 1 1], ...
            'contrast', 0.75, ...
            'spatialFrequencyCPD', [], ...
            'orientationDegs', 0, ...
            'spatialPhaseIncrementDegs', 30, ...
            'pixelSizeDegs', retinalImageResolutionDegs, ...
            'stimSizeDegs', max(sceneFOVdegs), ...
            'wavelengthSupport', displayGet(theDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(theDisplay, 'viewing distance') ...
            );

    
    spatialFrequenciesTested = [0.5 1 2 4 8 16 32 64];

    rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsMicrons,1); 
    
    for iFreq = 1:numel(spatialFrequenciesTested)
   
        stimParams.spatialFrequencyCPD = spatialFrequenciesTested(iFreq);
        fprintf('Generating scenes for the frames of the %2.1f c/deg pattern.\n', stimParams.spatialFrequencyCPD);

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
   
    
        parfor iFrame = 1:numel(theDriftingGratingFrameScenes)
            theScene = theDriftingGratingFrameScenes{iFrame};
            theMidgetRGCMosaicResponses(iFreq, iFrame,:) = theMidgetRGCmosaic.compute(theScene);
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
        'wavefrontSpatialSamples', 401, ...
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
