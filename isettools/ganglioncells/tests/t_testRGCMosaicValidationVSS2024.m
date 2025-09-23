 % Script to measure STF responses of macaque compute-ready mosaics 
% This script is also used to generate materials for the validation
% figures for the VSS2024 paper
%
% Usage:
%{
    t_testRGCMosaicValidationVSS2024    
%}

% Initialize session
close all; clear all;

% Actions to perform: compute responses
computeInputConeMosaicResponses = ~true;
computeMRGCMosaicResponses = ~true;
visualizeInputConeMosaicResponses = false;
visualizePSFonTopOfConeMosaic = false;

% Actions to perform: analyze computed responses
performCronerKaplanAnalysis = true;
performLeeShapleyAnalysis = ~true;

% Whether to realize the STF data or load previously analyzed data
reAnalyzeSTFData = true;

if (performCronerKaplanAnalysis && performLeeShapleyAnalysis)
    error('Choose EITHER the Croner&Kaplan analysis OR the Lee&Shapley analysis.');
end


% Aggregation over multiple runs 
aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs = ~true;
aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities = false;

% Aggregated targetVisualSTFdescriptors
aggregatedTargetVisualSTFdescriptors = { ...
           'default' ...
           'x1.3 RsRcRatio' ...
        };

% Aggregated eccentricities
aggregatedRGCmosaicNames = {...
            'VSS2024TalkTemporal3DegsMosaic' ...
            'VSS2024TalkTemporal7DegsMosaic' };

% Macaque mosaic 
coneMosaicSpecies = 'macaque';

% Target eccentricity
rgcMosaicName = 'VSS2024TalkTemporal3DegsMosaic';  % For Lee&Shapley analyses
rgcMosaicName = 'VSS2024TalkTemporal7DegsMosaic';  % For Reid&Shapley analyses (their RFs were in the range 3-13 degs)

% Optics
opticsSubjectName = 'VSS2024TalkFirstSubject';
opticsSubjectName = 'VSS2024TalkSecondSubject';


% Target visualizeSTFdescriptor. Options are:
% - 'default' ...
% - 'x1.3 RsRcRatio' ...
% - 'x1.5 RsRcRatio' ...
% - 'x0.75 RsRcRatio' ...
targetVisualSTFdescriptorToOptimizeFor = 'default';
%targetVisualSTFdescriptorToOptimizeFor = 'x1.3 RsRcRatio';

% Default pStruct
allConditionPstructs{1} = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
            coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptorToOptimizeFor);


% Chromaticity for validation STF responses
chromaticityForSTFresponses = 'Achromatic';   % Choose between 'LconeIsolating', 'MconeIsolating', 'Achromatic'
coneFundamentalsOptimizedForStimPosition = false;

% Optics for validation STF responses
% Option 1. Native optics (what was used to optimize the surround),
opticsForSTFresponses = 'nativeOptics';
residualWithRespectToNativeOpticsDefocusDiopters = [];

% Option 2. residual defocus with respect to native optics
%opticsForSTFresponses = 'refractionResidualWithRespectToNativeOptics';
%residualWithRespectToNativeOpticsDefocusDiopters = 0.75;
%residualWithRespectToNativeOpticsDefocusDiopters = 0.50;
%residualWithRespectToNativeOpticsDefocusDiopters = 0.25;
%residualWithRespectToNativeOpticsDefocusDiopters =-0.25;

% Option 3. adaptive optics, 6 mm pupil, no LCA
%opticsForSTFresponses = 'adaptiveOptics6MM';

% Option 4.  adaptive optics with LCA
% opticsForSTFresponses = 'adaptiveOptics6MMwithLCA';
 
% Option 5. a custom refraction
% opticsForSTFresponses = 'customRefraction';
% customRefractionDiopters = -4;


if (performLeeShapleyAnalysis)
    % Pick the orientation for the STF analyses. 
    % This has an effect only if reAnalyzeSTFData == true
    
    % 0deg orientation 
    fixedOptimalOrientation = 0;

    % Random orientation.
    % We use this option, because according to my correspondense with Lisa Croner (email-June 11, 2024), 
    % they did not optimize orientation in their 1995 paper
    fixedOptimalOrientation = nan;

    % Orientation for which the STF at half max extends to the highest SF 
    fixedOptimalOrientation = [];

    % Target RGCs with surround specificity around 0.5 (i.e., 50/50 L/M cone net weight in the surround)
    targetedSurroundPurityRange = [0 1];            % [0.4 0.6];

    % Or, target RGCs with surround specificity < 0.25 (i.e., similar surround and center cone dominances)
    %targetedSurroundSpecificifyRange = [0.0 0.15];
    % Or, target RGCs with surround specificity > 0.75 (i.e., opposite surround and center cone dominances)
    %targetedSurroundSpecificifyRange = [0.85 1.0];

    targetedRadialEccentricityRange = [];           % Empty means full range
    targetedCenterConeNumerosityRange = [1 1];      % Only examine cells with a single cone in their RF center

    % Target RGCs with any center specificity 
    targetedCenterPurityRange = [ ];   
end

if (performCronerKaplanAnalysis)
    % Pick the orientation for the STF analyses. 
    % This has an effect only if reAnalyzeSTFData == true
    
    % 0deg orientation 
    fixedOptimalOrientation = 0;

    % Random orientation.
    % We use this option, because according to my correspondense with Lisa Croner (email-June 11, 2024), 
    % they did not optimize orientation in their 1995 paper
    fixedOptimalOrientation = nan;

    % Orientation for which the STF at half max extends to the highest SF 
    fixedOptimalOrientation = [];

    % Target RGCs with surround specificity around 0.5 (i.e., 50/50 L/M cone net weight in the surround)
    targetedSurroundPurityRange = [];           % [0.4 0.6];
    targetedRadialEccentricityRange = [4.9 5.5];       % Empty means full range
    targetedCenterConeNumerosityRange = [1 2];     % all center cone numerosities

    % Target RGCs with any center specificity 
    targetedCenterPurityRange = [1 1];   

    % Croner & Kaplan analysis for achromatic gratings always
    chromaticityForSTFresponses = 'Achromatic';
    coneFundamentalsOptimizedForStimPosition = false;

    % Croner & Kaplan analysis for native optics always 
    opticsForSTFresponses = 'nativeOptics';
    residualWithRespectToNativeOpticsDefocusDiopters = [];
end


% If we are aggregating data over multiple runs, assemble the aggregated conditions
if ((~reAnalyzeSTFData) && (aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs||aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities))
    % AGGREGATE analyzed data from different runs from multiple visual STFs
    % descriptors (e.g. with different target Rs/Rc)
    if (aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs)

        allConditionPstructs = cell(1, numel(aggregatedTargetVisualSTFdescriptors));
        for iRun = 1:numel(aggregatedTargetVisualSTFdescriptors)
            allConditionPstructs{iRun} = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
                            coneMosaicSpecies, opticsSubjectName, RGCmosaicName, aggregatedTargetVisualSTFdescriptors{iRun});
        end % iRun
    elseif (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
        % Target visual STF descriptor
        targetVisualSTFdescriptorToOptimizeFor = 'default'; % mean CK mean RsRcRatio, mean CK intSCratio (what was used in PLOSpaper)

        allConditionPstructs = cell(1, numel(aggregatedRGCmosaicNames));
        for iRun = 1:numel(aggregatedRGCmosaicNames)
            allConditionPstructs{iRun} = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
                            coneMosaicSpecies, opticsSubjectName, aggregatedRGCmosaicNames{iRun}, targetVisualSTFdescriptorToOptimizeFor);
        end % iRun
    end
end

% Configure a conservative parpool manager. This gives at least 8 GB RAM/core
%ASPPManager = AppleSiliconParPoolManager(12);
%
ASPPManager = AppleSiliconParPoolManager('half max');
%ASPPManager = AppleSiliconParPoolManager('conservative');

% Generate filenames
pStruct = allConditionPstructs{1};
if (computeInputConeMosaicResponses || computeMRGCMosaicResponses)

    [theSurroundConnectedMRGCMosaicFullFileName, theInputConeMosaicSTFResponsesFullFileName, ...
    theMRGCMosaicSTFResponsesFullFileName, surroundConnectedParamsStruct, targetVisualSTFmodifierStruct] = ...
        RGCMosaicAnalyzer.filepathFor.testRuns(...
            pStruct.whichEye, pStruct.whichZernikeDataBase, pStruct.whichSubjectID, ...
            pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
            pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
            pStruct.rgcMosaic.employRFCenterOverlappingMosaic, ...
            pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
            pStruct.customLMSconeDensities, ...
            pStruct.rgcMosaicSurroundOptimization.targetVisualSTFdescriptor, ...
            pStruct.rgcMosaicSurroundOptimization.optimizationStrategy, ...
            opticsForSTFresponses, residualWithRespectToNativeOpticsDefocusDiopters, ...
            chromaticityForSTFresponses, ...
            coneFundamentalsOptimizedForStimPosition, ...
            'STFresponses');
else
    
    if (performCronerKaplanAnalysis)
        [theSurroundConnectedMRGCMosaicFullFileName, theInputConeMosaicSTFResponsesFullFileName, ...
        theMRGCMosaicSTFResponsesFullFileName, surroundConnectedParamsStruct, targetVisualSTFmodifierStruct] = ...
            RGCMosaicAnalyzer.filepathFor.testRuns(...
                pStruct.whichEye, pStruct.whichZernikeDataBase, pStruct.whichSubjectID, ...
                pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
                pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
                pStruct.rgcMosaic.employRFCenterOverlappingMosaic, ...
                pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
                pStruct.customLMSconeDensities, ...
                pStruct.rgcMosaicSurroundOptimization.targetVisualSTFdescriptor, ...
                pStruct.rgcMosaicSurroundOptimization.optimizationStrategy, ...
                'nativeOptics', [], ...
                'Achromatic', ...
                coneFundamentalsOptimizedForStimPosition, ...
                'STFresponses');

    elseif (performLeeShapleyAnalysis)
        [theSurroundConnectedMRGCMosaicFullFileName, theInputConeMosaicSTFResponsesFullFileName, ...
        theMRGCMosaicSTFResponsesFullFileName, surroundConnectedParamsStruct, targetVisualSTFmodifierStruct] = ...
            RGCMosaicAnalyzer.filepathFor.testRuns(...
                pStruct.whichEye, pStruct.whichZernikeDataBase, pStruct.whichSubjectID, ...
                pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
                pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
                pStruct.rgcMosaic.employRFCenterOverlappingMosaic, ...
                pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
                pStruct.customLMSconeDensities, ...
                pStruct.rgcMosaicSurroundOptimization.targetVisualSTFdescriptor, ...
                pStruct.rgcMosaicSurroundOptimization.optimizationStrategy, ...
                opticsForSTFresponses, residualWithRespectToNativeOpticsDefocusDiopters, ...
                'Achromatic', ...
                coneFundamentalsOptimizedForStimPosition, ...
                'STFresponses');
    else
        error('One of ``performCronerKaplanAnalysis`` or ``performanceLeeShapleyAnalysis`` must be set to true')
    end
end

if (computeInputConeMosaicResponses) || (computeMRGCMosaicResponses)
    % Load theMRGCmosaic
    load(theSurroundConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');
    % Print info on how the surrounds where optimized
    %RGCMosaicConstructor.helper.utils.printStruct(theMRGCMosaic.rfSurroundConnectivityParams, 'rfSurroundConnectivityParams')
end


if (computeInputConeMosaicResponses)
    switch (opticsForSTFresponses) 
        case  'refractionResidualWithRespectToNativeOptics'
            RGCMosaicAnalyzer.compute.inputConeMosaicSTFsForStimulusChromaticityAndOptics(...
                theMRGCMosaic, opticsForSTFresponses, residualWithRespectToNativeOpticsDefocusDiopters, ...
                chromaticityForSTFresponses, ...
                coneFundamentalsOptimizedForStimPosition, ...
                theInputConeMosaicSTFResponsesFullFileName, theMRGCMosaicSTFResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic);

        case 'customRefraction'
            RGCMosaicAnalyzer.compute.inputConeMosaicSTFsForStimulusChromaticityAndOptics(...
                theMRGCMosaic, opticsForSTFresponses, customRefractionDiopters, ...
                chromaticityForSTFresponses, ...
                coneFundamentalsOptimizedForStimPosition, ...
                theInputConeMosaicSTFResponsesFullFileName, theMRGCMosaicSTFResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic);

        otherwise
            RGCMosaicAnalyzer.compute.inputConeMosaicSTFsForStimulusChromaticityAndOptics(...
                theMRGCMosaic, opticsForSTFresponses, [],  ...
                chromaticityForSTFresponses,  ...
                coneFundamentalsOptimizedForStimPosition, ...
                theInputConeMosaicSTFResponsesFullFileName, theMRGCMosaicSTFResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic);
    end % switch
end % if (computeInputConeMosaicResponses)

if (computeMRGCMosaicResponses)
    RGCMosaicConstructor.helper.simulateExperiment.spatialTransferFunction(...
        theMRGCMosaic, [], [],...
        theInputConeMosaicSTFResponsesFullFileName, ...
        theMRGCMosaicSTFResponsesFullFileName, ...
        'computeInputConeMosaicResponses', false, ...
        'computeMRGCMosaicResponses', true, ...
        'visualizeResponse', ~true);
    return;
end % if (computeMRGCMosaicResponses)

if (performLeeShapleyAnalysis)
    % Output file
    theLeeShapleyAnalysisFileName = strrep(theMRGCMosaicSTFResponsesFullFileName, 'STFresponses', 'demos/LeeShapleyAnalyses');
    % Compute randomized orientations
    load(theMRGCMosaicSTFResponsesFullFileName, 'theMRGCMosaic', 'stimParams');

    % Reset the random number generator so we get the same sequence for the same mosaic
    rngSeed = theMRGCMosaic.rgcsNum;
    rng(rngSeed);

    % Randomize orientations
    theRandomlySelectedAnalyzedOrientations = randi(numel(stimParams.orientationDegs), [1 theMRGCMosaic.rgcsNum]);

    % Do it !
    RGCMosaicAnalyzer.compute.LeeShapleySTFanalysis(...
            theMRGCMosaicSTFResponsesFullFileName, ...
            theLeeShapleyAnalysisFileName, ...
            targetedSurroundPurityRange, ...
            targetedRadialEccentricityRange, ...
            targetedCenterConeNumerosityRange, ...
            targetedCenterPurityRange, ...
            reAnalyzeSTFData, ...
            'fixedOptimalOrientation', theRandomlySelectedAnalyzedOrientations, ...
            'limitVisualizedCKdataToTheEccentricititesOfSyntheticCells', ~true);
    return;
end

if (performCronerKaplanAnalysis)
    % Output file
    theCronerKaplanAnalysisFileName = strrep(theMRGCMosaicSTFResponsesFullFileName, 'STFresponses', 'demos/CronerKaplanAnalyses');

    if (~aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs) && (~aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
        allMRGCMosaicSTFResponsesFullFileNames = theMRGCMosaicSTFResponsesFullFileName;
        allCronerKaplanAnalysisFileNames = theCronerKaplanAnalysisFileName;
    else
        % Generate STFResponse filenames for the aggregated runs
        allMRGCMosaicSTFResponsesFullFileNames = {};

        if (aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs)
            % Generate filenames for all aggregated runs
            for iRun = 1:numel(allConditionPstructs)
                pStructForAggregatedRun = allConditionPstructs{iRun};
                pStructForAggregatedRun.rgcMosaicSurroundOptimization
                pStructForAggregatedRun.rgcMosaicSurroundOptimization.targetVisualSTFdescriptor

                % Generate targetVisualSTFmodifierStruct
                targetVisualSTFmodifierStruct = ...
                    RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateTargetVisualSTFmodifiersStruct(pStructForAggregatedRun.rgcMosaicSurroundOptimization.targetVisualSTFdescriptor);

                theRunParams = [...
                    targetVisualSTFmodifierStruct.surroundToCenterRcRatioMultiplier ...
                    targetVisualSTFmodifierStruct.surroundToCenterIntegratedSensitivityRatioMultiplier];

                allMRGCMosaicSTFResponsesFullFileNames{iRun} = strrep(theMRGCMosaicSTFResponsesFullFileName, 'vSTF_1.0_1.0', ...
                    sprintf('vSTF_%1.1f_%1.1f', theRunParams(1), theRunParams(2)));
                allCronerKaplanAnalysisFileNames{iRun} = strrep(theCronerKaplanAnalysisFileName, 'vSTF_1.0_1.0', ...
                    sprintf('vSTF_%1.1f_%1.1f', theRunParams(1), theRunParams(2)));
            end % iRun
        end  % if (aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs)

        if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
            targetEccString = sprintf('Ecc%2.1f_%2.1f_Size%2.1fx%2.1f', ...
                pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs(1), ...
                pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs(2), ...
                pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs(1), ...
                pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs(2));

            for iRun = 1:numel(allConditionPstructs)
                pStructForAggregatedRun = allConditionPstructs{iRun};
                aggregatedEccString = sprintf('Ecc%2.1f_%2.1f_Size%2.1fx%2.1f', ...
                    pStructForAggregatedRun.rgcMosaicSurroundOptimization.mosaicEccDegs(1), ...
                    pStructForAggregatedRun.rgcMosaicSurroundOptimization.mosaicEccDegs(2), ...
                    pStructForAggregatedRun.rgcMosaicSurroundOptimization.mosaicSizeDegs(1), ...
                    pStructForAggregatedRun.rgcMosaicSurroundOptimization.mosaicSizeDegs(2));

                tmp1 = strrep(theMRGCMosaicSTFResponsesFullFileName, targetEccString, aggregatedEccString);
                tmp2 = strrep(theCronerKaplanAnalysisFileName, targetEccString, aggregatedEccString);

                allMRGCMosaicSTFResponsesFullFileNames{iRun} = strrep(tmp1, ...
                    pStruct.rgcMosaicSurroundOptimization.optimizationStrategy, ...
                    pStructForAggregatedRun.rgcMosaicSurroundOptimization.optimizationStrategy);
                allCronerKaplanAnalysisFileNames{iRun} = strrep(tmp2, ...
                    pStruct.rgcMosaicSurroundOptimization.optimizationStrategy, ...
                    pStructForAggregatedRun.rgcMosaicSurroundOptimization.optimizationStrategy);
            end % for iRun
        end % if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
    end

    % Do it !
    RGCMosaicAnalyzer.compute.CronerAndKaplanSTFanalysis(...
            allMRGCMosaicSTFResponsesFullFileNames, ...
            allCronerKaplanAnalysisFileNames, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities, ...
            targetedSurroundPurityRange, ...
            targetedRadialEccentricityRange, ...
            targetedCenterConeNumerosityRange, ...
            targetedCenterPurityRange, ...
            reAnalyzeSTFData, ...
            'fixedOptimalOrientation', fixedOptimalOrientation, ...
            'limitVisualizedCKdataToTheEccentricititesOfSyntheticCells', ~true);
    return;
end % if (performCronerKaplanAnalysis)