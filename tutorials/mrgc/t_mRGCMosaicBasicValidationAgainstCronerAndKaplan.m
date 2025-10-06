

function t_mRGCMosaicBasicValidationAgainstCronerAndKaplan

    % Initialize
    clearvars;
    close all;

    hPositions = [2 3 4 5 6 7 8 9 10 11 12];
    vPositions = [0];
    [X,Y] = meshgrid(hPositions, vPositions);
    horizontalPosDegs = X(:);
    verticalPosDegs = Y(:);


    % (A) Prebaked mRGCmosaic: 11x11 deg synthesized at 7 degrees along the nasal meridian
    mosaicParams.sizeDegs = [11 11];
    mosaicParams.eccDegs  = [7 0];

    % (B) Synthesis params
    % RF center: maximize RF center spatial homogenity
    mosaicParams.spatialCompactnessSpectralPurityTradeoff = 1;

    % RF surround: narror tolerance for visual STF, H1 params in the bottom half of the macaque H1 range
    mosaicParams.surroundOptimizationSubString = 'PackerDacey2002H1freeLowH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0';


    % (C) Optics under which the mosaic was optimized
    opticsParams.ZernikeDataBase = 'Polans2015';
    opticsParams.subjectRankOrder = 9; 
    opticsParams.type = 'nativeOptics';
    opticsParams.refractiveErrorDiopters = [];
    opticsParams.visualizePSFonTopOfConeMosaic = ~true;


    % (D) STF params 
    minSF = 0.01; maxSF = 70; sfSamplesNum = 15;

    stfParams = struct();
    stfParams.sfSupport = logspace(log10(minSF), log10(maxSF), sfSamplesNum);
    stfParams.orientationDeltaDegs = 45;
    stfParams.spatialPhaseIncrementDegs = 30;

    stfParams.backgroundChromaticity = [0.301 0.301];
    stfParams.backgroundLuminanceCdM2 = 40.0;
    stfParams.chromaticity = 'Achromatic';   % Choose between 'LconeIsolating', 'MconeIsolating', 'Achromatic'
    stfParams.coneFundamentalsOptimizedForStimPosition = true;

    
    % Which computation stage to perform

    % 1. Compute the input cone mosaic STF responses 
    computeInputConeMosaicResponses = ~true;

    % 2. Compute the mRGC mosaic STFresponses
    computeMRGCMosaicResponses = ~true;

    % 3. Fit DoG models to the STF responses
    reAnalyzeSTFData = ~true;

    % 4. Analyze previously analyzed mRGC STF responses
    aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs = false;
    aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities = true;




    % Set the stage for aggregating analyed STF data from multiple mosaics
    onlyReturnAggegatedFilenames = false;
    if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities) 
        computeInputConeMosaicResponses = false;
        computeMRGCMosaicResponses = false;
        reAnalyzeSTFData = false;
        onlyReturnAggegatedFilenames = true;
    end
    

    % Subdirectory for exporting the generated PDFs
    exportVisualizationPDFdirectory = 'mosaicValidationPDFs';


    if (computeInputConeMosaicResponses || computeMRGCMosaicResponses)
        % Use a smaller number of the available CPU cores
        AppleSiliconParPoolManager('conservative');
    else
        AppleSiliconParPoolManager('extreme');
    end

    if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
        allMRGCMosaicSTFResponsesFullFileNames = cell(1, numel(horizontalPosDegs));
        allCronerKaplanAnalysisFullFileNames = cell(1, numel(horizontalPosDegs));
    end

    % Analyze all mosaic patches
    for iPos = 1:numel(horizontalPosDegs)

        % Update mosaic params with current crop params
        mosaicParams.cropParams = struct(...
            'sizeDegs', [1 1], ...
            'eccentricityDegs', [horizontalPosDegs(iPos) verticalPosDegs(iPos)]);
    
        % Go !
        [theMRGCMosaicSTFResponsesFullFileName, ...
         theCronerKaplanAnalysisFullFileName] = runValidation(mosaicParams, opticsParams, ...
                stfParams, ...
                computeInputConeMosaicResponses, ...
                computeMRGCMosaicResponses, ...
                reAnalyzeSTFData, ...
                aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
                false, ...
                onlyReturnAggegatedFilenames, ...
                [], [], ...
                exportVisualizationPDFdirectory);

        if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
            allMRGCMosaicSTFResponsesFullFileNames{iPos} = theMRGCMosaicSTFResponsesFullFileName;
            allCronerKaplanAnalysisFullFileNames{iPos} = theCronerKaplanAnalysisFullFileName;
        end
    end % iPos

    % Aggregate the previously analyzed data
    if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
        runValidation(mosaicParams, opticsParams, ...
                stfParams, ...
                computeInputConeMosaicResponses, ...
                computeMRGCMosaicResponses, ...
                reAnalyzeSTFData, ...
                aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
                aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities, ...
                false, ...
                allMRGCMosaicSTFResponsesFullFileNames, ...
                allCronerKaplanAnalysisFullFileNames, ...
                exportVisualizationPDFdirectory);

        return;
    end % if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)



end

function [theMRGCMosaicSTFResponsesFullFileName, ...
          theCronerKaplanAnalysisFullFileName] = runValidation(...
            mosaicParams, opticsParams, ...
            stfParams, ...
            computeInputConeMosaicResponses, ...
            computeMRGCMosaicResponses, ...
            reAnalyzeSTFData, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities, ...
            onlyReturnFilenames, ...
            theAggregatedMRGCMosaicSTFResponsesFullFileNames, ...
            theAggregatedCronerKaplanAnalysisFullFileNames, ...
            exportVisualizationPDFdirectory)


    % Load the desired mRGCmosaic and generated the optics for the computation
    [theMRGCmosaic, opticsForSTFresponses, thePSF, ...
        prebakedMRGCMosaicDir, prebakedMRGCMosaicFilename] = mRGCMosaic.loadPrebakedMosaic(mosaicParams, opticsParams);
    
    postFix = sprintf('%s_Ecc_%2.1f_%2.1f_Size_%2.1f_%2.1f', prebakedMRGCMosaicFilename, theMRGCmosaic.eccentricityDegs(1), theMRGCmosaic.eccentricityDegs(2), theMRGCmosaic.sizeDegs(1), theMRGCmosaic.sizeDegs(2));
    
    fprintf('\n Analyzing mosaic at [%2.1f, %2.1f] ...\n', theMRGCmosaic.eccentricityDegs(1), theMRGCmosaic.eccentricityDegs(2));
       
    % Filenames for intermediate responses
    p = getpref('isetbio');
    intermediateDataDir = p.rgcResources.intermediateDataDir;

    
    theInputConeMosaicSTFResponsesFileName = fullfile('scratchspace', sprintf('inputConeMosaicSTFresponses_%s.mat', postFix));
    theInputConeMosaicSTFResponsesFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
        intermediateDataDir, theInputConeMosaicSTFResponsesFileName, ...
        'generateMissingSubDirs', true);

    theMRGCMosaicSTFResponsesFileName = fullfile('scratchspace', sprintf('mRGCMosaicSTFresponses_%s.mat', postFix));
    theMRGCMosaicSTFResponsesFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
        intermediateDataDir, theMRGCMosaicSTFResponsesFileName, ...
        'generateMissingSubDirs', true);

    theCronerKaplanAnalysisFileName = fullfile('scratchspace', sprintf('CronerKaplanAnalysis_%s.mat', postFix));
    theCronerKaplanAnalysisFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
        intermediateDataDir, theCronerKaplanAnalysisFileName, ...
        'generateMissingSubDirs', true);

    if (onlyReturnFilenames)
        return;
    end

    if (computeInputConeMosaicResponses || computeMRGCMosaicResponses)
        % Generate the PSD data struct (with vLambda-weighted PSF) for visualization 
        thePSFData.data = RGCMosaicAnalyzer.compute.vLambdaWeightedPSF(thePSF);
        thePSFData.supportXdegs = thePSF.supportX/60;
        thePSFData.supportYdegs = thePSF.supportY/60;

        % Visualize the mosaic
        visualizeTheMosaic(theMRGCmosaic, thePSFData, exportVisualizationPDFdirectory, postFix);

        
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


        % Update stfParams adding 
        % - the stimulus resolution
        % - the stimulus position
        % - the stimulus size
        stfParams.resolutionDegs = stimulusResolutionDegs;
        stfParams.positionDegs = theMRGCmosaic.eccentricityDegs;
        stfParams.sizeDegs = 1.1*max(theMRGCmosaic.inputConeMosaic.sizeDegs);


        % Whether to visualize the computed input cone mosaic STF responses
        visualizeInputConeMosaicResponses = ~true;

        % Go !
        RGCMosaicConstructor.helper.simulateExperiment.spatialTransferFunction(...
            theMRGCmosaic, opticsForSTFresponses, ...
            stfParams, ...
            theInputConeMosaicSTFResponsesFullFileName, ...
            theMRGCMosaicSTFResponsesFullFileName, ...
            'computeInputConeMosaicResponses', computeInputConeMosaicResponses, ...
            'computeMRGCMosaicResponses', computeMRGCMosaicResponses , ...
            'visualizeResponse', visualizeInputConeMosaicResponses);

        return;
    end % if (computeInputConeMosaicResponses || computeMRGCMosaicResponses)

    
    % Choose orientation to analyze (slice through the 2D STF)
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


    % Visualization options
    % Whether to visualize the computed STF together with the cone weights
    % map and the 1D profile of cone weights map along with the STF-derived RF map
    % separately for each mRGC in the mosaic
    visualizeSTFwithConeWeightsMap = ~true;

    % Whether visualize the full 2D STF data computed along with the chosen STF slice
    visualizeFullAndMaximalExcursionSTF = false;

    % Whether to visualize the DoG model fit to the chosen STF slice
    visualizeSTFfits = false;

    % Whether to visualize the sinusoidal fit to the time-series response of the mRGCs
    visualizeSinusodalFitsToResponseTimeSeries = false;

    % Whether to only display the macaque STF data from the Croner & Kaplan study
    onlyDisplayCronerKaplanData = false;

    
    % If we are aggregating data over multiple runs, detemine the data
    % files to be included in the aggregate data

    if (~isempty(theAggregatedCronerKaplanAnalysisFullFileNames))
        theMRGCMosaicSTFResponsesFullFileName = theAggregatedMRGCMosaicSTFResponsesFullFileNames;
        theCronerKaplanAnalysisFullFileName = theAggregatedCronerKaplanAnalysisFullFileNames;
    end

    % Do it !
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
            exportVisualizationPDFdirectory, ...
            'visualizeFullAndMaximalExcursionSTF', visualizeFullAndMaximalExcursionSTF, ...
            'visualizeSTFfits', visualizeSTFfits, ...
            'visualizeSTFwithConeWeightsMap', visualizeSTFwithConeWeightsMap, ...
            'visualizeModelFitting', visualizeSinusodalFitsToResponseTimeSeries, ...
            'onlyDisplayCronerKaplanData', onlyDisplayCronerKaplanData, ...
            'fixedOptimalOrientation', fixedOptimalOrientation, ...
            'deltaThresholdForLimitingFittedSTFtoPrimaryPeak', 0.01, ...
            'limitVisualizedCKdataToTheEccentricititesOfSyntheticCells', ~true);

end


function visualizeTheMosaic(theMRGCmosaic, thePSFData, exportVisualizationPDFdirectory, postFix)
    
    % Visualize the mosaic of mRGC RF centers
    % identifying cones that are pooled by the RF center mechanism with
    % a weight >= mRGCMosaic.sensitivityAtPointOfOverlap;
    % This representation is like the representation used in visualizing 
    % mosaics of RGCs in typical in-vitro experiments (e.g. by the Chichilnisky lab)
    minCenterConeWeight = mRGCMosaic.sensitivityAtPointOfOverlap;

    % Visualization limits and ticks
    visualizedWidthDegs = theMRGCmosaic.sizeDegs(1);
    visualizedHeightDegs = theMRGCmosaic.sizeDegs(2);
    domainVisualizationLimits(1:2) = theMRGCmosaic.eccentricityDegs(1) + 0.51 * visualizedWidthDegs * [-1 1];
    domainVisualizationLimits(3:4) = theMRGCmosaic.eccentricityDegs(2) + 0.51 * visualizedHeightDegs * [-1 1];
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
        'withSuperimposedPSF', thePSFData, ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'plotTitle', sprintf('min center weight visualized: %2.3f', minCenterConeWeight));

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);


    % The full path where the generated PDFs will be stored can be found like so:
    p = getpref('isetbio');
    pdfExportRootDir = p.rgcResources.figurePDFsDir;

    theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('Mosaic_%s.pdf', postFix));
    % Generate the path if we need to
    thePDFFullFileName  = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);

    NicePlot.exportFigToPDF(thePDFFullFileName, hFig, 300);

end
