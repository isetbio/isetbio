

function t_mRGCMosaicBasicValidationAgainstCronerAndKaplan

    % Initialize
    clearvars;
    close all;

    % Use 6 CPU cores
    AppleSiliconParPoolManager(4);

    % Load an 11x11 deg mosaic that was synthesized at 7 degrees along the nasal meridian
    mosaicParams.eccDegs  = [7 0];
    mosaicParams.sizeDegs = [11 11];


    % (B) Surround optimization method
    mosaicParams.spatialCompactnessSpectralPurityTradeoff = 1;
    mosaicParams.surroundOptimizationSubString = 'PackerDacey2002H1freeLowH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0';


    % (C) Optics under which the mosaic was optimized
    opticsParams.ZernikeDataBase = 'Polans2015';
    opticsParams.subjectRankOrder = 9; 
    opticsParams.type = 'nativeOptics';
    opticsParams.refractiveErrorDiopters = [];
    opticsParams.visualizePSFonTopOfConeMosaic = ~true;


    % Subdirectory for exporting the generated PDFs
    exportVisualizationPDFdirectory = 'mosaicValidationPDFs';

    
    positionDegs = [5 0];
    sizeDegs = [1 1];

    % Params for the spatial transfer function
    minSF = 0.01; maxSF = 30; sfSamplesNum = 8;
    stfParams.sfSupport = logspace(log10(minSF), log10(maxSF), sfSamplesNum);
    stfParams.orientationDeltaDegs = 45;
    stfParams.spatialPhaseIncrementDegs = 30;

    theMRGCmosaic = computeSTF(mosaicParams, opticsParams, stfParams, positionDegs, sizeDegs, ...
        exportVisualizationPDFdirectory);

end

function theMRGCmosaic = computeSTF(mosaicParams, opticsParams, stfParams, positionDegs, sizeDegs, exportVisualizationPDFdirectory)

    % Crop a small patch (2x2 degs) of the mosaic, centered at 7 degrees
    mosaicParams.cropParams = struct(...
        'sizeDegs', sizeDegs, ...
        'eccentricityDegs', positionDegs ...
        );

    % Load the desired mRGCmosaic and generated the optics for the computation
    [theMRGCmosaic, opticsForSTFresponses, thePSF] = mRGCMosaic.loadPrebakedMosaic(mosaicParams, opticsParams);
    
    % Visualize the mosaic
    visualizeTheMosaic(theMRGCmosaic, exportVisualizationPDFdirectory);


    % Determine the stimulus pixel resolution to be a fraction of the minimum cone aperture or cone spacing in the mosaic
    % here, half of the cone spacing
    theMetric = 'cone aperture';  % choose from {'cone aperture' or cone spacing'}
    if ((strcmp(opticsForSTFresponses, 'adaptiveOptics6MM')) || (strcmp(opticsForSTFresponses, 'adaptiveOptics6MMwithLCA')))
        theFraction = 0.1;
    else
        theFraction = 0.25;
    end
    targetRGCindices =  1:theMRGCmosaic.rgcsNum;
    stimulusResolutionDegs = RGCMosaicConstructor.helper.simulateExperiment.stimulusResolutionFromConeApertureOrConeSpacing(...
                theMRGCmosaic, targetRGCindices, theFraction, theMetric);


    % Chromaticity for validation STF responses
    chromaticityForSTFresponses = 'Achromatic';   % Choose between 'LconeIsolating', 'MconeIsolating', 'Achromatic'
    coneFundamentalsOptimizedForStimPosition = true;

    % Params struct for the STF
    STFparamsStruct = struct(...
        'backgroundChromaticity', [0.301 0.301], ...
        'backgroundLuminanceCdM2', 40.0, ...
        'chromaticity', 'Achromatic', ...
        'coneFundamentalsOptimizedForStimPosition', true, ...
        'resolutionDegs', stimulusResolutionDegs, ...                       % to be determined separately for each optimization position
        'sfSupport', stfParams.sfSupport, ...
        'orientationDeltaDegs', stfParams.orientationDeltaDegs, ...
        'spatialPhaseIncrementDegs', stfParams.spatialPhaseIncrementDegs, ...
        'positionDegs', theMRGCmosaic.eccentricityDegs, ...  
        'sizeDegs', 1.1*max(theMRGCmosaic.inputConeMosaic.sizeDegs));

    % Filenames for intermediate responses
    p = getpref('isetbio');
    intermediateDataDir = p.rgcResources.intermediateDataDir;

    theInputConeMosaicSTFResponsesFileName = fullfile('scratchspace', 'inputConeMosaicSTFresponses.mat');
    theInputConeMosaicSTFResponsesFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
        intermediateDataDir, theInputConeMosaicSTFResponsesFileName, ...
        'generateMissingSubDirs', true);

    theMRGCMosaicSTFResponsesFileName = fullfile('scratchspace', 'mRGCMosaicSTFresponses.mat');
    theMRGCMosaicSTFResponsesFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
        intermediateDataDir, theMRGCMosaicSTFResponsesFileName, ...
        'generateMissingSubDirs', true);

    theCronerKaplanAnalysisFileName = fullfile('scratchspace', 'CronerKaplanAnalysis.mat');
    theCronerKaplanAnalysisFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
        intermediateDataDir, theCronerKaplanAnalysisFileName, ...
        'generateMissingSubDirs', true);

    % Whether to visualize the computed input cone mosaic STF responses
    visualizeInputConeMosaicResponses = ~true;

    % Compute the input cone mosaic STF responses and the mRGC mosaic STF
    % responses
    computeInputConeMosaicResponses = ~true;
    computeMRGCMosaicResponses = ~true;
    
    RGCMosaicConstructor.helper.simulateExperiment.spatialTransferFunction(...
        theMRGCmosaic, opticsForSTFresponses, ...
        STFparamsStruct, ...
        theInputConeMosaicSTFResponsesFullFileName, ...
        theMRGCMosaicSTFResponsesFullFileName, ...
        'computeInputConeMosaicResponses', computeInputConeMosaicResponses, ...
        'computeMRGCMosaicResponses', computeMRGCMosaicResponses , ...
        'visualizeResponse', visualizeInputConeMosaicResponses);


    reAnalyzeSTFData = true;
    aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs = false;
    aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities = false;
    onlyDisplayCronerKaplanData = false;

    % Choose orientation to analyze
    % 0 deg orientation 
    fixedOptimalOrientation = 0;

    % Random orientation.
    % We use this option, because according to my correspondense with Lisa Croner (email-June 11, 2024), 
    % they did not optimize orientation in their 1995 paper
    fixedOptimalOrientation = nan;

    % Orientation for which the STF at half max extends to the highest SF 
    fixedOptimalOrientation = [];

    % Target RGCs with surround specificity around 0.5 (i.e., 50/50 L/M cone net weight in the surround)
    targetedSurroundPurityRange = [];           % [0.4 0.6];
    targetedRadialEccentricityRange = [];       % Empty means full range
    targetedCenterConeNumerosityRange = [];     % [1 1], etc, empty means all numerosities
    targetedCenterPurityRange = [];

    [RsToRcVarianceCK, intStoCsensVarianceCK, RsToRcVariance, intStoCsensVariance] = ...
        RGCMosaicAnalyzer.compute.CronerAndKaplanSTFanalysis(...
            theMRGCMosaicSTFResponsesFullFileName, ...
            theCronerKaplanAnalysisFullFileName, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities, ...
            targetedSurroundPurityRange, ...
            targetedRadialEccentricityRange, ...
            targetedCenterConeNumerosityRange, ...
            targetedCenterPurityRange, ...
            reAnalyzeSTFData, ...
            'onlyDisplayCronerKaplanData', onlyDisplayCronerKaplanData, ...
            'fixedOptimalOrientation', fixedOptimalOrientation, ...
            'deltaThresholdForLimitingFittedSTFtoPrimaryPeak', 0.01, ...
            'limitVisualizedCKdataToTheEccentricititesOfSyntheticCells', ~true);

end


function visualizeTheMosaic(theMRGCmosaic, exportVisualizationPDFdirectory)
    % The full path where the generated PDFs will be stored can be found like so:
    % p = getpref('isetbio');
    % fullPath = fullfile(p.rgcResources.figurePDFsDir,exportVisualizationPDFdirectory);

    % Visualize the mosaic of mRGC RF centers
    % identifying cones that are pooled by the RF center mechanism with
    % a weight >= mRGCMosaic.sensitivityAtPointOfOverlap;
    % This representation is like the representation used in visualizing 
    % mosaics of RGCs in typical in-vitro experiments (e.g. by the Chichilnisky lab)
    minCenterConeWeight = mRGCMosaic.sensitivityAtPointOfOverlap;

    % Visualization limits and ticks
    visualizedWidthDegs = theMRGCmosaic.inputConeMosaic.sizeDegs(1);
    visualizedHeightDegs = theMRGCmosaic.inputConeMosaic.sizeDegs(2);
    domainVisualizationLimits(1:2) = theMRGCmosaic.eccentricityDegs(1) + 0.5 * visualizedWidthDegs * [-1 1];
    domainVisualizationLimits(3:4) = theMRGCmosaic.eccentricityDegs(2) + 0.5 * visualizedHeightDegs * [-1 1];
    domainVisualizationTicks = struct(...
        'x', theMRGCmosaic.eccentricityDegs(1) + 0.5 * visualizedWidthDegs * [-1 -0.5 0 0.5 1], ...
        'y', theMRGCmosaic.eccentricityDegs(2) + 0.5 * visualizedHeightDegs * [-1 -0.5 0 0.5 1]);

    ff = PublicationReadyPlotLib.figureComponents('1x1 giant rectangular-wide mosaic');

    hFig = figure(1); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    theMRGCmosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'identifyInputCones', true, ...
        'identifyPooledCones', true, ...
        'inputConesAlpha', 0.5, ...
        'identifiedConeAperture', 'lightCollectingArea4sigma', ...
        'identifiedConeApertureThetaSamples', 16, ...
        'minConeWeightVisualized', minCenterConeWeight, ...
        'centerSubregionContourSamples', 32, ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'plotTitle', sprintf('min center weight visualized: %2.3f', minCenterConeWeight), ...
        'visualizationPDFfileName', sprintf('mRGCmosaic_Ecc_%2.1f_%2.1f', theMRGCmosaic.eccentricityDegs(1), theMRGCmosaic.eccentricityDegs(2)), ...
        'exportVisualizationPDF', true, ...
        'exportVisualizationPDFdirectory', exportVisualizationPDFdirectory);

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

end
