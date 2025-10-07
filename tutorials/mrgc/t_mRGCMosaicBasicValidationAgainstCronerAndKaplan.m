function t_mRGCMosaicBasicValidationAgainstCronerAndKaplan(options)
% Validate synthetic mRGC mosaic against the in-vivo data of Croner&Kaplan '95
%
% Syntax:
%   t_mRGCMosaicBasicValidationAgainstCronerAndKaplan()
%
% Description:
%   Demonstrates how to Validate synthetic mRGC mosaic against the in-vivo data of Croner&Kaplan '95
%   by computing visual space STFs for all neurons in the examined mRGCs
%   mosaics, fitting DoG models to each of the computed STF, and comparing
%   the fitted DoG model parameters of the synthetic mRGCs against the
%   population of macaque mRGCs from Croner & Kaplan '95
%
%  This is set up with key/value pairs that demonstate how to select different
%  options. Different choices are illustrated in the examples
%  in the source code.
%
% Optional key/value pairs
%    See source code arguments block for a list of key/value pairs.

% History:
%    08/28/25  NPC  Wrote it.

% Examples:
%{

    
    % Aggregate previously analyzed data over specific XY eccentricities
    t_mRGCMosaicBasicValidationAgainstCronerAndKaplan(...
        'validationEccentritiesDegsX', [2 12], ...
        'validationEccentritiesDegsY', [-2 -1 0 1 2], ...
        'aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities', true);
%}


arguments
    % ---- Name encoding properties of the rgcMosaic, such as its eccentricity ---
    % See RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters
    % for what is available and to add new mosaics
    options.rgcMosaicName (1,:) char = 'PLOSpaperNasal7DegsMosaic';

    % ---- Which species to employ ----
    % Choose between {'macaque', 'human'}. If 'macaque' is chosen, the input
    % cone mosaic has a 1:1 L/M cone ratio.
    options.coneMosaicSpecies  (1,:) char {mustBeMember(options.coneMosaicSpecies,{'human','macaque'})} = 'human';


    % ----- Which subject optics to employ -----
    options.opticsSubjectName (1,:) char = 'PLOSpaperDefaultSubject';


    % ------ targetVisualSTF options ----
    % Options are : {'default', 'x1.3 RsRcRatio'}
    % These are with respect to the macaque data of the Croner & Kaplan '95 study
    % 'default': target the mean Rs/Rc, and the mean Ks/Kc (Rs/Rc)^2
    % 'x1.3 RsRcRatio': target Rs/Rc ratio that is 1.3 x mean, and target Ks/Kc (Rs/Rc)^2: mean
    options.targetVisualSTFdescriptor (1,:) char = 'default';


    % Whether to close previously open figures
    options.closeOpenFigures (1,1) logical = true;

    % ------- Validation options ------
    % Validation positions within the prebaked mosaic specified by options.rgcMosaicName 
    options.validationEccentritiesDegsX (1,:) double = [2  12];
    options.validationEccentritiesDegsY (1,:) double = [0];
    options.validationSizeDegs (1,2) double = [1 1];

    % STF params for the validation
    options.validationSFsupportMin (1,1) double = 0.01;
    options.validationSFsupportMax (1,1) double = 70;
    options.validationSFsupportMaxSamplesNum (1,1) double = 15;

    options.validationOrientationDeltaDegs (1,1) double = 45;
    options.validationSpatialPhaseIncrementDegs (1,1) double = 30;

    options.validationBackgroundChromaticity (1,2) double = [0.301 0.301];
    options.validationBackgroundLuminanceCdM2 (1,1) double = 40.0;
    options.validationChromaticity  (1,:) char =  'Achromatic';   % Choose between 'LconeIsolating', 'MconeIsolating', 'Achromatic'
    options.validationConeFundamentalsOptimizedForStimPosition (1,1) logical = true;


    % ------ Choices of actions to perform -------
    % 1. Compute the input cone mosaic STF responses 
    options.computeInputConeMosaicResponses (1,1) logical = false

    % 2. Compute the mRGC mosaic STFresponses
    options.computeMRGCMosaicResponses (1,1) logical = false

    % 3. Fit DoG models to the STF responses
    options.reAnalyzeSTFData (1,1) logical = false

    % 4. Aggregate previously analyzed mRGC STF responses
    options.aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs (1,1) logical = false
    options.aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities (1,1) logical = false



end % arguments


% Close previously open figures
closePreviouslyOpenFigures = options.closeOpenFigures;

if (closePreviouslyOpenFigures)
    % Close any stray figs
    close all;
end


% Set flags from key/value pairs
rgcMosaicName = options.rgcMosaicName;
coneMosaicSpecies = options.coneMosaicSpecies;
opticsSubjectName = options.opticsSubjectName;
targetVisualSTFdescriptor = options.targetVisualSTFdescriptor;


% Generate pStruct with synthesized mosaic params
pStruct = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
    coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor);

% Extract mosaic and optics params for the original prebaked mosaic
[mosaicParams, opticsParams] = RGCMosaicConstructor.helper.utils.extractSynthesizedMosaicAndOpticsParams(...
    pStruct, targetVisualSTFdescriptor);

% Validation eccentricities and patch size within the prebaked mosaic
validationEccentritiesDegsX = options.validationEccentritiesDegsX;
validationEccentritiesDegsY = options.validationEccentritiesDegsY;
validationSizeDegs = options.validationSizeDegs; 


%STF params 
stfParams = struct();
stfParams.sfSupport = logspace(log10(options.validationSFsupportMin), log10(options.validationSFsupportMax), options.validationSFsupportMaxSamplesNum);
stfParams.orientationDeltaDegs = options.validationOrientationDeltaDegs;
stfParams.spatialPhaseIncrementDegs = options.validationSpatialPhaseIncrementDegs;

stfParams.backgroundChromaticity = options.validationBackgroundChromaticity;
stfParams.backgroundLuminanceCdM2 = options.validationBackgroundLuminanceCdM2;
stfParams.chromaticity = options.validationChromaticity;
stfParams.coneFundamentalsOptimizedForStimPosition = options.validationConeFundamentalsOptimizedForStimPosition;


% Actions to perform
computeInputConeMosaicResponses = options.computeInputConeMosaicResponses;
computeMRGCMosaicResponses = options.computeMRGCMosaicResponses;
reAnalyzeSTFData = options.reAnalyzeSTFData;
aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs = options.aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs;
aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities = options.aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities;


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

% Grid it
[X,Y] = meshgrid(validationEccentritiesDegsX, validationEccentritiesDegsY);
validationEccentritiesDegsXgrid = X(:); validationEccentritiesDegsYgrid = Y(:);

if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
    allMRGCMosaicSTFResponsesFullFileNames = cell(1, numel(validationEccentritiesDegsXgrid));
    allCronerKaplanAnalysisFullFileNames = cell(1, numel(validationEccentritiesDegsYgrid));
end

    validationEccentritiesDegsXgrid
% Analyze all validation positions
for iPos = 1:numel(validationEccentritiesDegsXgrid)

    % Update mosaic params with current crop params
    mosaicParams.cropParams = struct(...
        'sizeDegs', validationSizeDegs, ...
        'eccentricityDegs', [validationEccentritiesDegsXgrid(iPos) validationEccentritiesDegsYgrid(iPos)]);

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

% Aggregate the previously analyzed data over all validation positions
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
