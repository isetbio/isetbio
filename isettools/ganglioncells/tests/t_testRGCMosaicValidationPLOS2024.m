% Script to examine the properties (visual space referred STFs)  of compute-ready mosaics 
% This script is also used to generate materials for the validation
% figures for the PLOS2024 paper
%
% Usage:
%{
    t_testRGCMosaicValidationPLOS2024    
%}

% Initialize session
close all; clear all;

% Actions to perform: compute STF responses
computeInputConeMosaicResponses = ~true;
computeMRGCMosaicResponses = ~true;
visualizeInputConeMosaicResponses = false;
visualizePSFonTopOfConeMosaic = ~true;

% Actions to perform: analyze computed STF responses
performCronerKaplanAnalysis = true;

% Whether to realize the STF data or load previously analyzed data
reAnalyzeSTFData = ~true;

% Aggregation over multiple runs 
aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs = ~true;
aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities = true;

onlyDisplayCronerKaplanData = true;

% Human mosaic 
coneMosaicSpecies = 'human';


% Optics subject for which the mRGC mosaics were optimized for
opticsSubjectName = 'PLOSpaperDefaultSubject';  % This is subject #2 (ranked #3)
%opticsSubjectName = 'VSS2024TalkFirstSubject';  % This is subject #3 (ranked #7)

%opticsSubjectName = 'PLOSpaperStrehlRatio_0.21';  % Still very high RcDegs at fovea
%opticsSubjectName = 'PLOSpaperStrehlRatio_0.19';  % Lets try this one

% Single condition pStruct
% Eccentricity: select from
%   'PLOSpaperFovealMosaic' ...
%   'PLOSpaperTemporal2DegsMosaic' ...
%   'PLOSpaperTemporal4DegsMosaic' ...
%   'PLOSpaperTemporal7DegsMosaic' ...
%   'PLOSpaperTemporal10DegsMosaic' ...
%   'PLOSpaperTemporal14DegsMosaic' ...
%   'PLOSpaperTemporal19DegsMosaic' ...
%   'PLOSpaperTemporal25DegsMosaic' ...
%   'PLOSpaperTemporal32DegsMosaic' ...

% Target eccentricity
% Run with mRGCMosaic.amplificationInCenterOnlySensitivityCausedByInactiveSurrounds set to 0.7 (from 0.85)
%RGCmosaicName = 'PLOSpaperTemporal32DegsMosaicLowerOverlap';
%RGCmosaicName = 'PLOSpaperTemporal25DegsMosaicLowerOverlap';
%RGCmosaicName = 'PLOSpaperTemporal19DegsMosaicLowerOverlap';
%RGCmosaicName = 'PLOSpaperTemporal14DegsMosaicLowerOverlap';
%RGCmosaicName = 'PLOSpaperTemporal10DegsMosaicLowerOverlap';
%RGCmosaicName = 'PLOSpaperTemporal7DegsMosaicLowerOverlap';
%RGCmosaicName = 'PLOSpaperTemporal4DegsMosaicLowerOverlap';
%RGCmosaicName = 'PLOSpaperTemporal2DegsMosaicLowerOverlap';
%RGCmosaicName = 'PLOSpaperFovealMosaicLowerOverlap';


% Run with
% mRGCMosaic.amplificationInCenterOnlySensitivityCausedByInactiveSurrounds 0.6 
% (best agreement in RF center overlap with Gauthier)
%RGCmosaicName = 'PLOSpaperTemporal32DegsMosaic';  % Fully Done. Also done with second  subject optics
%RGCmosaicName = 'PLOSpaperTemporal25DegsMosaic';  % Fully done. Also done with second  subject optics
%RGCmosaicName = 'PLOSpaperTemporal19DegsMosaic';  % Fully done Also done with second  subject optics
%RGCmosaicName = 'PLOSpaperTemporal14DegsMosaic';  % Fully done Also done with second  subject optics 
%RGCmosaicName = 'PLOSpaperTemporal10DegsMosaic';  % Fully done Also done with second  subject optics 
%RGCmosaicName = 'PLOSpaperTemporal7DegsMosaic';   % Fully done. Optimizing L cone on Crete now
%RGCmosaicName = 'PLOSpaperTemporal4DegsMosaic';   % Fully done. Also done with second  subject optics
RGCmosaicName = 'PLOSpaperTemporal2DegsMosaic';   % Fully done. Last one Computing input cone mosaic on Leviahtnan now
%RGCmosaicName = 'PLOSpaperFovealMosaic';           % Fully done. Also done with second  subject optics

% Default pStruct
targetVisualSTFdescriptorToOptimizeFor = 'default';
allConditionPstructs{1} = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
            coneMosaicSpecies, opticsSubjectName, RGCmosaicName, targetVisualSTFdescriptorToOptimizeFor);

% Chromaticity for validation STF responses
chromaticityForSTFresponses = 'Achromatic';   % Choose between 'LconeIsolating', 'MconeIsolating', 'Achromatic'
coneFundamentalsOptimizedForStimPosition = false;

% Optics for validation STF response
% Option 1. Native optics (what was used to optimize the surround)
opticsForSTFresponses = 'nativeOptics';
residualWithRespectToNativeOpticsDefocusDiopters = [];


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
    targetedRadialEccentricityRange = [];       % Empty means full range
    targetedCenterConeNumerosityRange = [];     % single cone centers only
    targetedCenterPurityRange = [];
end

% If we are aggregating data over multiple runs, assemble the aggregated conditions
if ((~reAnalyzeSTFData) && (aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs||aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities))
    % AGGREGATE analyzed data from different runs from multiple visual STFs
    % descriptors (e.g. with different target Rs/Rc)
    if (aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs)

        % Aggregated targetVisualSTFdescriptors
        aggregatedTargetVisualSTFdescriptors = { ...
           'default' ...
           'x1.3 RsRcRatio' ...
           'x1.5 RsRcRatio' ...
           'x0.75 RsRcRatio' ...
        };

        allConditionPstructs = cell(1, numel(aggregatedTargetVisualSTFdescriptors));
        for iRun = 1:numel(aggregatedTargetVisualSTFdescriptors)
            allConditionPstructs{iRun} = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
                            coneMosaicSpecies, opticsSubjectName, RGCmosaicName, aggregatedTargetVisualSTFdescriptors{iRun});
        end % iRun

    elseif (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
     
        % Target visual STF descriptor
        targetVisualSTFdescriptorToOptimizeFor = 'default'; % mean CK mean RsRcRatio, mean CK intSCratio (what was used in PLOSpaper)


        % Fully done for 'PLOSpaperDefaultSubject', which is subject #2(ranked #3), and 
        % also for 'VSS2024TalkFirstSubject', which is subject #3 (ranked #7)
        aggregatedRGCmosaicNames = {...
            'PLOSpaperTemporal32DegsMosaic' ...
            'PLOSpaperTemporal25DegsMosaic' ...
            'PLOSpaperTemporal19DegsMosaic' ...
            'PLOSpaperTemporal14DegsMosaic' ...
            'PLOSpaperTemporal10DegsMosaic' ...
            'PLOSpaperTemporal7DegsMosaic' ...
            'PLOSpaperTemporal4DegsMosaic' ...
            'PLOSpaperTemporal2DegsMosaic' ...
            'PLOSpaperFovealMosaic' ...
        };

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


if (computeInputConeMosaicResponses) || (computeMRGCMosaicResponses)
    % Load theMRGCmosaic
    load(theSurroundConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');
    % Print info on how the surrounds where optimized
    RGCMosaicConstructor.helper.utils.printStruct(theMRGCMosaic.rfSurroundConnectivityParams, 'rfSurroundConnectivityParams');
    RGCMosaicConstructor.helper.utils.printStruct(theMRGCMosaic.surroundVarianceInComputeReadyMosaic, 'surroundVarianceInComputeReadyMosaic');
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


if (performCronerKaplanAnalysis)
    % Output file
    theCronerKaplanAnalysisFileName = strrep(theMRGCMosaicSTFResponsesFullFileName, 'STFresponses', 'demos/CronerKaplanAnalyses');

    if (~aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs) && (~aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
        allMRGCMosaicSTFResponsesFullFileNames = theMRGCMosaicSTFResponsesFullFileName;
        allCronerKaplanAnalysisFileNames = theCronerKaplanAnalysisFileName;
    else
        % Accumulate all filenames
        allMRGCMosaicSTFResponsesFullFileNames = {};

        if (aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs)
            % Generate filenames for all aggregated runs
            for iRun = 1:numel(allConditionPstructs)
                pStructForAggregatedRun = allConditionPstructs{iRun};

                theRunParams = [...
                    pStructForAggregatedRun.rgcMosaicSurroundOptimization.targetVisualSTFdescriptor.surroundToCenterRcRatioMultiplier ...
                    pStructForAggregatedRun.rgcMosaicSurroundOptimization.targetVisualSTFdescriptor.surroundToCenterIntegratedSensitivityRatioMultiplier];

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
    [RsToRcVarianceCK, intStoCsensVarianceCK, RsToRcVariance, intStoCsensVariance] = ...
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
            'onlyDisplayCronerKaplanData', onlyDisplayCronerKaplanData, ...
            'fixedOptimalOrientation', fixedOptimalOrientation, ...
            'deltaThresholdForLimitingFittedSTFtoPrimaryPeak', 0.01, ...
            'limitVisualizedCKdataToTheEccentricititesOfSyntheticCells', ~true);

end % if (performCronerKaplanAnalysis)


