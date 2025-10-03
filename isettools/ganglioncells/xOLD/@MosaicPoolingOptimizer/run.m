function run(restartParPool)

    if (nargin == 0)
        restartParPool = false;
    end

    [rawFiguresRoot, scaledFiguresRoot] = initializeRun(restartParPool);

    % Ask user which mosaic to use
    [mosaicHorizontalEccentricityDegs,  mosaicEccsForSummaryStatistics] = ...
        MosaicPoolingOptimizer.chooseMosaicToUse();

    % Get mosaic params
    mosaicParams = MosaicPoolingOptimizer.getMosaicParams(mosaicHorizontalEccentricityDegs);

    % This controls the visualized spatial support range
    tickSeparationArcMin = 6;
    normalizedPeakSurroundSensitivity = 0.4;
    visualizedSpatialFrequencyRange = [0.1 100];

    %MosaicPoolingOptimizer.performComputeVisuallyProjectedAnatomicalConeRcsByFittingRFmap([]);
    
    % Get operation to perform
    operationSetToPerformContains = MosaicPoolingOptimizer.operationsMenu(mosaicParams);
  
    % Perform the generateRGCMosaic operation
    if (operationSetToPerformContains.generateCenterConnectedRGCMosaic)
        if (isfield(mosaicParams, 'maxConeInputsPerRGCToConsiderTransferToNearbyRGCs'))
            MosaicPoolingOptimizer.performGenerateCenterConnectedRGCMosaicOp(mosaicParams, ...
                'maxConeInputsPerRGCToConsiderTransferToNearbyRGCs', mosaicParams.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs);
        else
            MosaicPoolingOptimizer.performGenerateCenterConnectedRGCMosaicOp(mosaicParams);
        end
        return;
    end

    % Perform the visualizeRGCMosaic operation
    if (operationSetToPerformContains.visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs)
        
        xMin = mosaicParams.eccDegs(1) - 1.0;
        xMax = mosaicParams.eccDegs(1) + 1.0;
        yMin = mosaicParams.eccDegs(2) - 0.35;
        yMax = mosaicParams.eccDegs(2) + 0.35;


        identifyPooledCones = true;
        identifyInputCones = ~true;
        plotRFoutlines = true;
        labelRetinalMeridians  = ~true;
        backgroundColor = [1 1 1];
        domainVisualizationLimits = [xMin xMax yMin yMax];
        domainVisualizationTicks = struct('x', -30:0.2:30, 'y', -20:0.2:20);

        if (mosaicParams.eccDegs(1)<0)
            reverseXDir = true;
        else
            reverseXDir = false;
        end

        MosaicPoolingOptimizer.performVisualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCsOp(mosaicParams, ...
            'identifyPooledCones', identifyPooledCones, ...
            'identifyInputCones', identifyInputCones, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'plotRFoutlines', plotRFoutlines, ...
            'labelRetinalMeridians', labelRetinalMeridians, ...
            'backgroundColor', backgroundColor, ...
            'reverseXDir', reverseXDir);
        return;
    end

    % Perform the visualizePSFsWithinRGCMosaic operation
    if (operationSetToPerformContains.visualizePSFsWithinRGCMosaic)
        MosaicPoolingOptimizer.performVisualizePSFsWithinRGCMosaicOp(mosaicParams, tickSeparationArcMin);
        return;
    end

     % Perform the visualizeConePSFsAtLocationWithinRGCMosaic operation
    if (operationSetToPerformContains.visualizeConePSFsAtLocationWithinRGCMosaic)
        MosaicPoolingOptimizer.performVisualizeConePSFsWithinRGCMosaicOp(mosaicParams, 3.0);
        return;
    end

    % Perform the computeInputConeMosaicSTFresponses operation
    if (operationSetToPerformContains.computeInputConeMosaicSTFresponses)
        MosaicPoolingOptimizer.performComputeInputConeMosaicSTFresponsesOp(mosaicParams);
        return;
    end

    % Perform the computeInputConeMosaicVisuallyProjectedCharacteristicRadii operation
    if (operationSetToPerformContains.computeInputConeMosaicVisuallyProjectedCharacteristicRadii)
        MosaicPoolingOptimizer.performComputeInputConeMosaicVisuallyProjectedRcDegsOp(mosaicParams);
        return;
    end

    % Perform the optimizeSurroundConePoolingModels operation
    if (operationSetToPerformContains.optimizeSurroundConePoolingModels)
        
        % Fitting options
        multiStartsNumRetinalPooling = 20; % 12
        multiStartsNumDoGFit = 128;

        % More weight for matching the Rs/Rc ratio
        rmseWeightForRsRcResidual = 2.0;
        rmseWeightForSCintSensResidual = 1.0;

        MosaicPoolingOptimizer.performOptimizeSurroundConePoolingModelsOp(mosaicParams, ...
            'multiStartsNumRetinalPooling', multiStartsNumRetinalPooling, ...
            'multiStartsNumDoGFit', multiStartsNumDoGFit, ...
            'rmseWeightForRsRcResidual', rmseWeightForRsRcResidual , ...
            'rmseWeightForSCintSensResidual', rmseWeightForSCintSensResidual);

        return;
    end


    % Perfom the inspectOptimizedSurroundConePoolingModels operation
    if (operationSetToPerformContains.inspectOptimizedSurroundConePoolingModels)
        MosaicPoolingOptimizer.performInspectOptimizedSurroundConePoolingModelsOp(...
            mosaicParams, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
            'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange, ...
            'gridlessLineWeightingFunctions', true);

        return;
    end

    % Perfom the summarizeOptimizedSurroundConePoolingModels operation
    if (operationSetToPerformContains.summarizeOptimizedSurroundConePoolingModels)
        MosaicPoolingOptimizer.performSummarizeOptimizedSurroundConePoolingModelsOp(...
            mosaicEccsForSummaryStatistics, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
            'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange);

        return;
    end

    % Perform the generateComputeReadyMidgetRGCMosaic operation
    if (operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.performGenerateComputeReadyMidgetRGCMosaicOp(mosaicParams);
        return;
    end

    % Perform the computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic operation
    if (operationSetToPerformContains.computeVisualSTFsOfTheComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.performComputeVisualSTFsOfTheComputeReadyMidgetRGCMosaicOp(mosaicParams);
        return;
    end

    % Perform the fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic operation
    if (operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.performFitVisualSTFsOfTheComputeReadyMidgetRGCMosaicOp(mosaicParams);
        return;
    end

    % Perform the adjustGainOfComputeReadyMRGCMosaicBasedOnVisualSTFdata
    if (operationSetToPerformContains.adjustGainOfComputeReadyMidgetRGCMosaicBaseOnFittedVisualSTFs)
        MosaicPoolingOptimizer.performAdjustGainOfComputeReadyMRGCMosaicBasedOnVisualSTFdata(mosaicParams)
        return;
    end

    % Perform the resetGainOfComputeReadyMRGCMosaic
    if (operationSetToPerformContains.resetGainOfComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.performResetGainOfComputeReadyMRGCMosaic(mosaicParams);
        return;
    end

    % Perform the visualizeConePoolingRFmapAndVisualSTFforTargetRGC operation
    if (operationSetToPerformContains.visualizeConePoolingRFmapAndVisualSTFforTargetRGC)
        
        if (mosaicParams.eccDegs(1)<0)
            reverseXDir = true;
        else
            reverseXDir = false;
        end
        reverseXDir = false;

        MosaicPoolingOptimizer.performVisualizeConePoolingRFmapAndVisualSTFforTargetRGC(mosaicParams, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
            'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange, ...
            'reverseXDir', reverseXDir, ...
            'gridlessLineWeightingFunctions', true);
        return;
    end

    % Perform the visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic operation
    if (operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic)
        MosaicPoolingOptimizer.performVisualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic(mosaicParams);
        return;
    end

    % Perform the visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic operation
    if (operationSetToPerformContains.visualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaics)
        MosaicPoolingOptimizer.performVisualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic(mosaicEccsForSummaryStatistics);
        return;
    end

    % Perform the computeVisualRFsUsingSubSpaceMapping operation
    if (operationSetToPerformContains.computeVisualRFsUsingSubspaceMapping)

        reComputeInputConeMosaicResponses = ~true;
        reComputeMRGCMosaicResponses = ~true;
        recomputeRFs = ~true;
        visualizeOptimallyMappedRFmapLocations = true;
        stimPositionDegs = [0.65 0.72];
        stimSizeDegs = [0.35 0.35];

        % Upper SF limit for the subspace stimuli
        maxSFLimit = 90;
        % or leave it empty to use the optimal max SF
        % maxSFLimit = []; 

        % Upper SF limit for the stimuli used to generate the RF
        % This can be empty or lower than the maxSFLimit
        maxSFanalyzed = []; %input('Enter max SF to include in subspace analysis. Hit enter for full range: ');


        rfMappingPixelMagnificationFactor = 1.0;

        msequencePixelSizeArcMin = 0.25*4;
        msequencePixelSizeDegs = msequencePixelSizeArcMin/60;

        visualizedRGCindex = input('Enter RGC index to visualize its subspace RF. Hit enter to visualize all:');

        MosaicPoolingOptimizer.performComputeVisualRFsUsingSubSpaceMapping(mosaicParams, ...
            'stimSizeDegs', stimSizeDegs, ....
            'stimPositionDegs', stimPositionDegs, ...
            'maxSFLimit', maxSFLimit, ...
            'maxSFToBeAnalyzed', maxSFanalyzed, ...
            'rfMappingPixelMagnificationFactor', rfMappingPixelMagnificationFactor, ...
            'reComputeInputConeMosaicResponses', reComputeInputConeMosaicResponses, ...
            'reComputeMRGCMosaicResponses', reComputeMRGCMosaicResponses, ...
            'reComputeRFs', recomputeRFs, ...
            'visualizeOptimallyMappedRFmapLocations', visualizeOptimallyMappedRFmapLocations, ...
            'visualizedRGCindex', visualizedRGCindex, ...
            'msequencePixelSizeDegs', msequencePixelSizeDegs);
        return;
    end

    % Perform the computeVisualRFsUsingMSequenceMapping operation
    if (operationSetToPerformContains.computeVisualRFsUsingMSequenceMapping)

        reComputeInputConeMosaicResponses = true;
        reComputeMRGCMosaicResponses = true;
        recomputeRFs = true;
        visualizeOptimallyMappedRFmapLocations = true;

        % (0,0) mosaic
        stimPositionDegs = [0.65 0.72];
        stimSizeDegs = [0.2 0.2];

        % (2.5, 0 ) mosaic
        % stimPositionDegs = [3 -1];
        % stimPositionDegs = [3.5 -0.1];
        % stimPositionDegs = [2 0.4];
        % stimSizeDegs = [0.6 0.6];

        % Spatial sampling of RF (# of squares)
        rfPixelsAcross = 16;
        rfPixelsAcross = 44.0;

        % Bit length of m-sequence
        mSequenceBitLength = 12;

        % Use ternary instead of binary m-sequence
        ternaryInsteadOfBinaryMsequence = ~true;

        MosaicPoolingOptimizer.performComputeVisualRFsUsingMSequenceMapping(mosaicParams, ...
            'stimSizeDegs', stimSizeDegs, ....
            'stimPositionDegs', stimPositionDegs, ...
            'rfPixelsAcross', rfPixelsAcross, ...
            'mSequenceBitLength', mSequenceBitLength, ...
            'ternaryInsteadOfBinaryMsequence', ternaryInsteadOfBinaryMsequence, ...
            'reComputeInputConeMosaicResponses', reComputeInputConeMosaicResponses, ...
            'reComputeMRGCMosaicResponses', reComputeMRGCMosaicResponses, ...
            'reComputeRFs', recomputeRFs, ...
            'visualizeOptimallyMappedRFmapLocations', visualizeOptimallyMappedRFmapLocations ...
            );     
    end

 
    % Perform the ContrastSTFsAcrossDifferentOpticsOrChromaticities operation
    if (operationSetToPerformContains.contrastSTFsAcrossDifferentOpticsOrChromaticities)
        % 0.20 0.31 0.47 0.64 0.75
        targetRangeForSurroundConeMix = 0.20 + [0.00 0.05];
        %targetRangeForSurroundConeMix = 0.31 + [0.00 0.05];
        targetRangeForSurroundConeMix = 0.475 + [0.00 0.05];
        %targetRangeForSurroundConeMix = 0.64 + [0.00 0.05];
        %targetRangeForSurroundConeMix = 0.75 + [0.00 0.05];

        % Only cells with RF centers receiving a single type of cone inputs
        targetRangeForCenterMix = [1 1];

        MosaicPoolingOptimizer.performContrastSTFsAcrossDifferentChromaticities(...
            mosaicParams, rawFiguresRoot, scaledFiguresRoot, ...
            'performSurroundAnalysisForConesExclusiveToTheSurround', true, ...
            'targetRangeForSurroundConeMix', targetRangeForSurroundConeMix, ...
            'targetRangeForCenterMix', targetRangeForCenterMix);
    end

    % Perform the ContrastSTFsAcrossDifferentOpticsOrChromaticities operation
    if (operationSetToPerformContains.contrastMSequenceRFsAcrossDifferentOpticsOrChromaticities)
        % 0.20 0.31 0.47 0.64 0.75
        %targetRangeForSurroundConeMix = 0.20 + [0.00 0.05];
        %targetRangeForSurroundConeMix = 0.31 + [0.00 0.05];
        targetRangeForSurroundConeMix = 0.475 + [0.00 0.05];
        %targetRangeForSurroundConeMix = 0.64 + [0.00 0.05];
        %targetRangeForSurroundConeMix = 0.75 + [0.00 0.05];

        % Only cells with RF centers receiving a single type of cone inputs
        targetRangeForCenterMix = [1 1];

        % (0,0) mosaic
        stimPositionDegs = [0.65 0.72];
        stimSizeDegs = [0.2 0.2];

        % (2.5, 0 ) mosaic
        %stimPositionDegs = [3 -1];
        %stimPositionDegs = [2 0.4];
        %stimSizeDegs = [0.6 0.6];

        % Spatial sampling of RF (# of squares)
        rfPixelsAcross = 16;
        rfPixelsAcross = 32;
        rfPixelsAcross = 44;

        % Bit length of m-sequence
        mSequenceBitLength = 12;

        % Use ternary instead of binary m-sequence
        ternaryInsteadOfBinaryMsequence = ~true;

        MosaicPoolingOptimizer.performContrastMSequenceRFsAcrossDifferentChromaticities(...
            mosaicParams, ...
            stimPositionDegs, ternaryInsteadOfBinaryMsequence, ...
            mSequenceBitLength, rfPixelsAcross, ...
            rawFiguresRoot, scaledFiguresRoot, ...
            'performSurroundAnalysisForConesExclusiveToTheSurround', true, ...
            'targetRangeForSurroundConeMix', targetRangeForSurroundConeMix, ...
            'targetRangeForCenterMix', targetRangeForCenterMix);
    end

    % Perform the computeVisualRFcenterMapsViaDirectConvolutionWithPSF  operation
    if (operationSetToPerformContains.computeVisualRFcenterMapsViaDirectConvolutionWithPSF)
        MosaicPoolingOptimizer.performComputeVisualRFcenterMapsViaDirectConvolutionWithPSF()
    end

end

function [rawFiguresRoot, scaledFiguresRoot] = initializeRun(restartParPool)

    if (restartParPool)
        computerInfo = GetComputerInfo;
        switch (lower(computerInfo.localHostName))
            case 'ithaka'
                maxNumWorkers = 12;
            case 'crete'
                maxNumWorkers = 12;
            otherwise
                maxNumWorkers = [];
        end
    
        delete(gcp('nocreate'))
        c = parcluster('Processes');
        numWorkers = input('Enter number of parallel workers: ');
        if (isempty(numWorkers))
            numWorkers = maxNumWorkers;
        end
        if (~isempty(numWorkers))
            c.NumWorkers = numWorkers;
        end
        parpool(c);
    end

    % Generate path to rawFigures root directory.
    % These are figures generated for manuscripts/presentations, formatted
    % using a common figureFormat
    p = getpref('isetbio');
    rawFiguresRoot = strrep(p.alternateFullDataDir, ...
        'ISETBioValidationFiles/gradleFiles/validationFull', ...
        'ManuscriptSupportMaterials/PLOS2024/figures/raw');

    % Where scaled figures are stored
    scaledFiguresRoot = '/Users/nicolas/Documents/Manuscripts/PLOS2024/figures';

end