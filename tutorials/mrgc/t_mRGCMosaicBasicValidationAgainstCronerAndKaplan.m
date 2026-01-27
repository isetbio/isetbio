function t_mRGCMosaicBasicValidationAgainstCronerAndKaplan(options)
% Validate a synthesized mRGC mosaic against the in-vivo data of Croner&Kaplan '95
%
% Syntax:
%   t_mRGCMosaicBasicValidationAgainstCronerAndKaplan()
%
% Description:
%   Demonstrates how to validate synthetic mRGC mosaic against the in-vivo data of Croner&Kaplan '95
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

% UTTBSkip

% Skip running these examples during autovalidation because they require material
% not included in the ISETBio repository.

%
% NOTE: To run any RGC-related ISETBio code, such as this tutorial, users must follow
% the directions discribed in:
%    https://github.com/isetbio/isetbio/wiki/Retinal-ganglion-cell-(RGC)-mosaics
% under section "Configuring ISETBio to access RGC resources and run RGC simulations"
%

    % Aggregate previously analyzed data over specific XY eccentricities
    t_mRGCMosaicBasicValidationAgainstCronerAndKaplan(...
        'validationEccentricitiesDegsX', [2 3 4 5 6 7 8 9 10 11 12], ...
        'validationEccentricitiesDegsY', [-2 -1 0 1 2], ...
        'aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities', true);


    % Compute the requisite inputConeMosaic and rgcMosaic STF responses
    % for a mosaic synthesized at 2 degrees
    t_mRGCMosaicBasicValidationAgainstCronerAndKaplan(...
        'rgcMosaicName', 'PLOSpaperNasal2DegsTinyMosaic', ...
        'validationEccentricitiesDegsX', 2, ...
        'validationEccentricitiesDegsY', 0, ...
        'computeInputConeMosaicResponses', true, ...
        'computeMRGCMosaicResponses', true);

    % Analyze the mRGCMosaic STF responses computed previously
    % by fitting a DoG model to each RGC STF, and validate against the
    % Croner&Kaplan data
    t_mRGCMosaicBasicValidationAgainstCronerAndKaplan(...
        'rgcMosaicName', 'PLOSpaperNasal2DegsTinyMosaic', ...
        'validationEccentricitiesDegsX', 2, ...
        'validationEccentricitiesDegsY', 0, ...
        'reAnalyzeSTFData', true);


    % Just plot previously analyzed and fitted STF validations
    t_mRGCMosaicBasicValidationAgainstCronerAndKaplan(...
        'rgcMosaicName', 'PLOSpaperNasal2DegsTinyMosaic', ...
        'validationEccentricitiesDegsX', 2, ...
        'validationEccentricitiesDegsY', 0, ...
        'reAnalyzeSTFData', ~true);

    % Reanalyze the STF data for a random orientation (nan), 
    % optimal, [], or 90 degs
    analyzedSTFsliceOrientation = [];
    %analyzedSTFsliceOrientation = nan;
    %%analyzedSTFsliceOrientation = 90;

    % Re-analyze
    t_mRGCMosaicBasicValidationAgainstCronerAndKaplan(...
        'reAnalyzeSTFData', true, ...
        'validationEccentricitiesDegsX', [2 3 4 5 6 7 8 9 10 11 12], ...
        'validationEccentricitiesDegsY', [0], ...
        'analyzedSTFsliceOrientation', analyzedSTFsliceOrientation, ...
        'aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities', ~true);

    % Plot only
    t_mRGCMosaicBasicValidationAgainstCronerAndKaplan(...
        'reAnalyzeSTFData', ~true, ...
        'validationEccentricitiesDegsX', [2 3 4 5 6 7 8 9 10 11 12], ...
        'validationEccentricitiesDegsY', [0], ...
        'analyzedSTFsliceOrientation', analyzedSTFsliceOrientation, ...
        'aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities', true);

%}



arguments
    % ---- Mosaic specifiers for selecting a prebaked mRGC mosaic ---
    
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



    % ------- Validation options ------
    % Validation positions within the prebaked mosaic specified by options.rgcMosaicName 
    options.validationEccentricitiesDegsX (1,:) double = [2 3 4 5 6 7 8 9 10 11 12];
    options.validationEccentricitiesDegsY (1,:) double = [0];
    options.validationSizeDegs (1,2) double = [1 1];

    % STF stimulus params for the validation
    options.validationSFsupportMin (1,1) double = 0.01;
    options.validationSFsupportMax (1,1) double = 70;
    options.validationSFsupportMaxSamplesNum (1,1) double = 15;

    options.validationSpatialPhaseIncrementDegs (1,1) double = 30;
    options.validationOrientationDeltaDegs (1,1) double = 45;

    options.validationBackgroundChromaticity (1,2) double = [0.301 0.301];
    options.validationBackgroundLuminanceCdM2 (1,1) double = 40.0;
    options.validationChromaticity  (1,:) char =  'Achromatic';   % Choose between 'LconeIsolating', 'MconeIsolating', 'Achromatic'
    options.validationConeFundamentalsOptimizedForStimPosition (1,1) logical = true;

    % STF analysis params
    options.analyzedSTFsliceOrientation = [];   
    % Option 1: nan.  A random orientation for each cell
    % Using this option may simulate the experimental conditions of Croner&Kaplan as, 
    % according to NPC's correspondense with Lisa Croner (email-June 11, 2024), 
    % they did not optimize orientation in their 1995 paper
    %
    % Option 2: []. The orientation at which the cell's STF has a maximum extension at the high SF regime (default)
    %
    % Option 3: A specific value in 0:validationOrientationDeltaDegs:(180-validationOrientationDeltaDegs)}

    % Criteria for selecting RGCs for which to analyze their STFs 
    options.analyzedCenterPurityRange = [];          % Empty means all, or select a specific range, e.g., [0.4 0.6]
    options.analyzedSurroundPurityRange = [];        % Empty means all, or select a specific range, e.g., [0.4 0.6]
    options.analyzedRadialEccentricityRange = [];    % Empty means full range, or select an eccentricity range, e.g., [1 2]
    options.analyzedCenterConeNumerosityRange = [];  % Empty means all numerositities, or select a numerosity range, e.g., [1 1], [1 3] etc.
    

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



    % ----- Visualization options -----
    % Whether to visualize the mRGC mosaic and PSF at each validation position
    options.visualizeMosaicsAndPSFsAtValidationEccentricities (1,1) logical = false

    % Whether to visualize and contrast the custom cone fundamentals 
    % to the SS2 cone fundamentals. This occurs only if
    % validationConeFundamentalsOptimizedForStimPosition is set to true
    options.visualizeCustomConeFundamentals (1,1) logical = false

    % Whether to visualize the computed input cone mosaic STF responses
    options.visualizeInputConeMosaicResponses (1,1) logical = false

    % Whether to visualize the computed STF together with the cone weights
    % map and the 1D profile of cone weights map along with the STF-derived RF map
    % separately for each mRGC in the mosaic
    options.visualizeSTFwithConeWeightsMap (1,1) logical = false

    % Whether visualize the full 2D STF data computed along with the chosen STF slice
    options.visualizeFullAndMaximalExcursionSTF (1,1) logical = false

    % Whether to visualize the DoG model fit to the chosen STF slice
    options.visualizeSTFfits (1,1) logical = false

    % Whether to visualize the sinusoidal fit to the time-series response of the mRGCs
    options.visualizeSinusodalFitsToResponseTimeSeries (1,1) logical = false

    % Whether to only display the macaque STF data from the Croner & Kaplan study
    options.onlyDisplayCronerKaplanData (1,1) logical = false

    % Whether to close previously open figures
    options.closeOpenFigures (1,1) logical = true;

end % arguments


% Close previously open figures
closePreviouslyOpenFigures = options.closeOpenFigures;

if (closePreviouslyOpenFigures)
    % Close any stray figs
    close all;
end


% Set flags from key/value pairs

% Mosaic specifiers for selecting a prebaked mRGC mosaic
rgcMosaicName = options.rgcMosaicName;
coneMosaicSpecies = options.coneMosaicSpecies;
opticsSubjectName = options.opticsSubjectName;
targetVisualSTFdescriptor = options.targetVisualSTFdescriptor;

% Validation eccentricities and patch size within the prebaked mosaic
validationEccentricitiesDegsX = options.validationEccentricitiesDegsX;
validationEccentricitiesDegsY = options.validationEccentricitiesDegsY;
validationSizeDegs = options.validationSizeDegs; 


% STF stimulus params 
stfParams = struct();
stfParams.sfSupport = logspace(log10(options.validationSFsupportMin), log10(options.validationSFsupportMax), options.validationSFsupportMaxSamplesNum);
stfParams.spatialPhaseIncrementDegs = options.validationSpatialPhaseIncrementDegs;
stfParams.orientationDeltaDegs = options.validationOrientationDeltaDegs;


stfParams.backgroundChromaticity = options.validationBackgroundChromaticity;
stfParams.backgroundLuminanceCdM2 = options.validationBackgroundLuminanceCdM2;
stfParams.chromaticity = options.validationChromaticity;
stfParams.coneFundamentalsOptimizedForStimPosition = options.validationConeFundamentalsOptimizedForStimPosition;

% STF analysis params
stfAnalysisParams.analyzedSTFsliceOrientation = options.analyzedSTFsliceOrientation;
stfAnalysisParams.analyzedCenterPurityRange =  options.analyzedCenterPurityRange;
stfAnalysisParams.analyzedSurroundPurityRange = options.analyzedSurroundPurityRange;
stfAnalysisParams.analyzedRadialEccentricityRange = options.analyzedRadialEccentricityRange;
stfAnalysisParams.analyzedCenterConeNumerosityRange =  options.analyzedCenterConeNumerosityRange;
    

% Actions to perform
computeInputConeMosaicResponses = options.computeInputConeMosaicResponses;
computeMRGCMosaicResponses = options.computeMRGCMosaicResponses;
reAnalyzeSTFData = options.reAnalyzeSTFData;
aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs = options.aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs;
aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities = options.aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities;


% Visualization options
visualizeMosaicsAndPSFsAtValidationEccentricities = options.visualizeMosaicsAndPSFsAtValidationEccentricities;
visualizeCustomConeFundamentals = options.visualizeCustomConeFundamentals;
visualizeInputConeMosaicResponses = options.visualizeInputConeMosaicResponses;
visualizeSTFwithConeWeightsMap = options.visualizeSTFwithConeWeightsMap;
visualizeFullAndMaximalExcursionSTF = options.visualizeFullAndMaximalExcursionSTF;
visualizeSTFfits = options.visualizeSTFfits;
visualizeSinusodalFitsToResponseTimeSeries = options.visualizeSinusodalFitsToResponseTimeSeries;
onlyDisplayCronerKaplanData = options.onlyDisplayCronerKaplanData;


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


% Generate validation position grid
[X,Y] = meshgrid(validationEccentricitiesDegsX, validationEccentricitiesDegsY);
validationEccentricitiesDegsXgrid = X(:); validationEccentricitiesDegsYgrid = Y(:);

if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
    allMRGCMosaicSTFResponsesFullFileNames = cell(1, numel(validationEccentricitiesDegsXgrid));
    allCronerKaplanAnalysisFullFileNames = cell(1, numel(validationEccentricitiesDegsYgrid));
end

% Analyze all validation positions
for iPos = 1:numel(validationEccentricitiesDegsXgrid)

    % Update mosaic params with current crop params
    cropParams = struct(...
        'sizeDegs', validationSizeDegs, ...
        'eccentricityDegs', [validationEccentricitiesDegsXgrid(iPos) validationEccentricitiesDegsYgrid(iPos)]);

    % Go !
    [theMRGCMosaicSTFResponsesFullFileName, ...
     theCronerKaplanAnalysisFullFileName] = runValidation(...
            coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor, ...
            cropParams, ...
            stfParams, ...
            stfAnalysisParams, ...
            computeInputConeMosaicResponses, ...
            computeMRGCMosaicResponses, ...
            reAnalyzeSTFData, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
            false, ...
            onlyReturnAggegatedFilenames, ...
            [], [], ...
            visualizeMosaicsAndPSFsAtValidationEccentricities, ...
            visualizeCustomConeFundamentals, ...
            visualizeInputConeMosaicResponses, ...
            visualizeSTFwithConeWeightsMap, ...
            visualizeFullAndMaximalExcursionSTF, ...
            visualizeSTFfits, ...
            visualizeSinusodalFitsToResponseTimeSeries, ...
            onlyDisplayCronerKaplanData, ...
            exportVisualizationPDFdirectory);

    if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
        allMRGCMosaicSTFResponsesFullFileNames{iPos} = theMRGCMosaicSTFResponsesFullFileName;
        allCronerKaplanAnalysisFullFileNames{iPos} = theCronerKaplanAnalysisFullFileName;
    end
end % iPos

% Aggregate the previously analyzed data over all validation positions
if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
    runValidation(coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor, ...
        cropParams, ...
        stfParams, ...
        stfAnalysisParams, ...
        computeInputConeMosaicResponses, ...
        computeMRGCMosaicResponses, ...
        reAnalyzeSTFData, ...
        aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
        aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities, ...
        false, ...
        allMRGCMosaicSTFResponsesFullFileNames, ...
        allCronerKaplanAnalysisFullFileNames, ...
        visualizeMosaicsAndPSFsAtValidationEccentricities, ...
        visualizeCustomConeFundamentals, ...
        visualizeInputConeMosaicResponses, ...
        visualizeSTFwithConeWeightsMap, ...
        visualizeFullAndMaximalExcursionSTF, ...
        visualizeSTFfits, ...
        visualizeSinusodalFitsToResponseTimeSeries, ...
        onlyDisplayCronerKaplanData, ...
        exportVisualizationPDFdirectory);

    return;
end % if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)

end % t_mRGCMosaicBasicValidationAgainstCronerAndKaplan



function [theMRGCMosaicSTFResponsesFullFileName, ...
          theCronerKaplanAnalysisFullFileName] = runValidation(...
            coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor, ...
            cropParams, ...
            stfParams, ...
            stfAnalysisParams, ...
            computeInputConeMosaicResponses, ...
            computeMRGCMosaicResponses, ...
            reAnalyzeSTFData, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities, ...
            onlyReturnFilenames, ...
            theAggregatedMRGCMosaicSTFResponsesFullFileNames, ...
            theAggregatedCronerKaplanAnalysisFullFileNames, ...
            visualizeMosaicsAndPSFsAtValidationEccentricities, ...
            visualizeCustomConeFundamentals, ...
            visualizeInputConeMosaicResponses, ...
            visualizeSTFwithConeWeightsMap, ...
            visualizeFullAndMaximalExcursionSTF, ...
            visualizeSTFfits, ...
            visualizeSinusodalFitsToResponseTimeSeries, ...
            onlyDisplayCronerKaplanData, ...
            exportVisualizationPDFdirectory)


    % Load the mRGCmosaic specified by the passed parameters:
    %   coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor
    % and generate the optics that were used to synthesize the mosaic
    [theMRGCmosaic, opticsForSTFresponses, thePSF, ...
        prebakedMRGCMosaicDir, prebakedMRGCMosaicFilename] = mRGCMosaic.loadPrebakedMosaic(...
            coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor, ...
            'computeTheMosaicOptics', computeInputConeMosaicResponses||computeMRGCMosaicResponses, ...
            'cropParams', cropParams);

    postFix = sprintf('%s_Ecc_%2.1f_%2.1f_Size_%2.1f_%2.1f', prebakedMRGCMosaicFilename, theMRGCmosaic.eccentricityDegs(1), theMRGCmosaic.eccentricityDegs(2), theMRGCmosaic.sizeDegs(1), theMRGCmosaic.sizeDegs(2));
    
    if (computeInputConeMosaicResponses || computeMRGCMosaicResponses)
        fprintf('\nComputing STF responses at [%2.1f, %2.1f] ...\n', theMRGCmosaic.eccentricityDegs(1), theMRGCmosaic.eccentricityDegs(2));
    else
        fprintf('\nAnalyzing STF responses at [%2.1f, %2.1f] ...\n', theMRGCmosaic.eccentricityDegs(1), theMRGCmosaic.eccentricityDegs(2));
    end

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
        if (visualizeMosaicsAndPSFsAtValidationEccentricities)
            visualizeMosaicAndEmployedPSF(theMRGCmosaic, thePSFData, exportVisualizationPDFdirectory, postFix);
        end

        
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

        % Update stfParams for the current validation position, adding 
        % - the stimulus resolution
        % - the stimulus position = mosaic position
        % - the stimulus size = 1.1 x input cone mosaic size
        stfParams.resolutionDegs = stimulusResolutionDegs;
        stfParams.positionDegs = theMRGCmosaic.eccentricityDegs;
        stfParams.sizeDegs = 1.1*max(theMRGCmosaic.inputConeMosaic.sizeDegs);

        % Compute STF responses for the current validation position
        RGCMosaicConstructor.helper.simulateExperiment.spatialTransferFunction(...
            theMRGCmosaic, opticsForSTFresponses, ...
            stfParams, ...
            theInputConeMosaicSTFResponsesFullFileName, ...
            theMRGCMosaicSTFResponsesFullFileName, ...
            'computeInputConeMosaicResponses', computeInputConeMosaicResponses, ...
            'computeMRGCMosaicResponses', computeMRGCMosaicResponses , ...
            'visualizeResponse', visualizeInputConeMosaicResponses, ...
            'visualizeCustomConeFundamentals', visualizeCustomConeFundamentals && stfParams.coneFundamentalsOptimizedForStimPosition);

        return;
    end % if (computeInputConeMosaicResponses || computeMRGCMosaicResponses)

    % If we are aggregating data over multiple runs, detemine the data
    % files to be included in the aggregate data

    if (~isempty(theAggregatedCronerKaplanAnalysisFullFileNames))
        theMRGCMosaicSTFResponsesFullFileName = theAggregatedMRGCMosaicSTFResponsesFullFileNames;
        theCronerKaplanAnalysisFullFileName = theAggregatedCronerKaplanAnalysisFullFileNames;
    end

    % Analyze the computed STFs by fitting DoG models to extracting the
    % Croner & Kaplan parameters
    [RsToRcVarianceCK, intStoCsensVarianceCK, RsToRcVariance, intStoCsensVariance] = ...
        RGCMosaicAnalyzer.compute.CronerAndKaplanSTFanalysis(...
            theMRGCMosaicSTFResponsesFullFileName, ...
            theCronerKaplanAnalysisFullFileName, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities, ...
            stfAnalysisParams.analyzedSurroundPurityRange, ...
            stfAnalysisParams.analyzedRadialEccentricityRange, ...
            stfAnalysisParams.analyzedCenterConeNumerosityRange, ...
            stfAnalysisParams.analyzedCenterPurityRange, ...
            reAnalyzeSTFData, ...
            exportVisualizationPDFdirectory, ...
            'visualizeMosaic', visualizeMosaicsAndPSFsAtValidationEccentricities, ...
            'visualizeFullAndMaximalExcursionSTF', visualizeFullAndMaximalExcursionSTF, ...
            'visualizeSTFfits', visualizeSTFfits, ...
            'visualizeSTFwithConeWeightsMap', visualizeSTFwithConeWeightsMap, ...
            'visualizeModelFitting', visualizeSinusodalFitsToResponseTimeSeries, ...
            'onlyDisplayCronerKaplanData', onlyDisplayCronerKaplanData, ...
            'fixedOptimalOrientation', stfAnalysisParams.analyzedSTFsliceOrientation, ...
            'deltaThresholdForLimitingFittedSTFtoPrimaryPeak', 0.01, ...
            'limitVisualizedCKdataToTheEccentricititesOfSyntheticCells', ~true);

end


function visualizeMosaicAndEmployedPSF(theMRGCmosaic, thePSFData, exportVisualizationPDFdirectory, postFix)
    
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


    % The full path where the generated PDFs will be stored
    pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir;
    theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('Mosaic_%s_withPSF.pdf', postFix));

    % Generate the path if we need to
    thePDFFullFileName  = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);

    NicePlot.exportFigToPDF(thePDFFullFileName, hFig, 300, 'beVerbose');

end
