% This script is also used to generate materials for the validation
% figures for the PLOS2024 paper
% Usage:
%{
    t_validateAgainstInVitroData
%}

% Initialize session
close all; clear all;


theRunIndexToCompute = 56;
generateVisualRFandConePoolingMapComboPlots = true;

% Human mosaic 
coneMosaicSpecies = 'human';

targetVisualSTFdescriptorToOptimizeFor = 'default';

% Target eccentricity
rgcMosaicName = 'PLOSpaperTemporal32DegsMosaicLowerOverlap';

% Optics subject
opticsSubjectName = 'PLOSpaperDefaultSubject';

% Generate params struct
pStruct = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
    coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptorToOptimizeFor);

% RF mapping params
%adaptive optics, 6 mm pupil, no LCA
opticsForRFmapping = 'adaptiveOptics6MM';
residualWithRespectToNativeOpticsDefocusDiopters = [];


% Choose chromaticity to run
chromaticityForRFmapping = 'Achromatic';   % Choose between 'LconeIsolating', 'MconeIsolating', 'Achromatic'
coneFundamentalsOptimizedForStimPosition = false;

% Generate filenames
[theSurroundConnectedMRGCMosaicFullFileName, theInputConeMosaicRFmappingResponsesFullFileName, ...
    theMRGCMosaicRFmappingResponsesFullFileName, surroundConnectedParamsStruct, targetVisualSTFmodifierStruct] = ...
        RGCMosaicAnalyzer.filepathFor.testRuns(...
            pStruct.whichEye, ...
            pStruct.whichZernikeDataBase, ...
            pStruct.whichSubjectID, ...
            pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
            pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs,...
            pStruct.rgcMosaic.employRFCenterOverlappingMosaic, ...
            pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
            pStruct.customLMSconeDensities, ...
            pStruct.rgcMosaicSurroundOptimization.targetVisualSTFdescriptor, ...
            pStruct.rgcMosaicSurroundOptimization.optimizationStrategy, ...
            opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters, ...
            chromaticityForRFmapping, ...
            coneFundamentalsOptimizedForStimPosition, ...
            'mSequenceResponses');

theOriginalSurroundConnectedMRGCMosaicFullFileName = theSurroundConnectedMRGCMosaicFullFileName;
theOriginalMRGCMosaicRFmappingResponsesFullFileName = theMRGCMosaicRFmappingResponsesFullFileName;
theOriginalInputConeMosaicRFmappingResponsesFullFileName = theInputConeMosaicRFmappingResponsesFullFileName;

targetString = sprintf('Ecc%2.1f_%2.1f_Size%2.1fx%2.1f',...
    pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs(1), ...
    pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs(2), ...
    pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs(1), ...
    pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs(2));

if (~isempty(strfind(theOriginalMRGCMosaicRFmappingResponsesFullFileName, 'UpperH1paramsNarrowVisualSTFparamTolerance')))
    targetSurroundOptimizationStrategy = 'UpperH1paramsNarrowVisualSTFparamTolerance';
end

if (~isempty(strfind(theOriginalMRGCMosaicRFmappingResponsesFullFileName, 'MidH1paramsNarrowVisualSTFparamTolerance')))
    targetSurroundOptimizationStrategy  = 'MidH1paramsNarrowVisualSTFparamTolerance';
end

if (~isempty(strfind(theOriginalMRGCMosaicRFmappingResponsesFullFileName, 'LowH1paramsNarrowVisualSTFparamTolerance')))
    targetSurroundOptimizationStrategy  = 'LowH1paramsNarrowVisualSTFparamTolerance';
end


[surroundOptimizationStrategyToBeAggregated, ...
     targetNumberOfMappedCellsToBeAggregated, ...
     targetHorizontalEccentricitiesToBeAggregated, ...
     targetMosaicSizesToBeAggregated, ...
     targetMappedPositionDegsToBeAggregated ] = testsHelperFunctions.assembleEccentricityMatrices();

% Load the 32-deg mosaic
mosaicEccDegs(1) = targetHorizontalEccentricitiesToBeAggregated(theRunIndexToCompute);
mosaicEccDegs(2) = pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs(2);
aggregatedEccString = sprintf('Ecc%2.1f_%2.1f_Size%2.1fx%2.1f',mosaicEccDegs(1), mosaicEccDegs(2), ...
                targetMosaicSizesToBeAggregated(theRunIndexToCompute,1),  targetMosaicSizesToBeAggregated(theRunIndexToCompute,1));
     
theSurroundConnectedMRGCMosaicFullFileName = strrep(...
            strrep(theOriginalSurroundConnectedMRGCMosaicFullFileName, targetString, aggregatedEccString), ...
            targetSurroundOptimizationStrategy, surroundOptimizationStrategyToBeAggregated{theRunIndexToCompute});

theMRGCMosaicRFmappingResponsesFullFileName = strrep(...
            strrep(theOriginalMRGCMosaicRFmappingResponsesFullFileName, targetString, aggregatedEccString), ...
            targetSurroundOptimizationStrategy, surroundOptimizationStrategyToBeAggregated{theRunIndexToCompute});
    
theInputConeMosaicRFmappingResponsesFullFileName = strrep(...
            strrep(theOriginalInputConeMosaicRFmappingResponsesFullFileName, targetString, aggregatedEccString), ...
            targetSurroundOptimizationStrategy, surroundOptimizationStrategyToBeAggregated{theRunIndexToCompute});

% Load theMRGCmosaic
load(theSurroundConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');

% Add the mapped position info
mappedPositionDegs = targetMappedPositionDegsToBeAggregated(theRunIndexToCompute,:);
mappedSizeDegs = RGCMosaicAnalyzer.compute.mSequenceRFmapSizeDegsForMosaic(theMRGCMosaic, mappedPositionDegs, targetNumberOfMappedCellsToBeAggregated{theRunIndexToCompute});
        
mappedPositionString = sprintf('_@%2.1f_%2.1f.mat', mappedPositionDegs(1), mappedPositionDegs(2));
theMRGCMosaicRFmappingResponsesFullFileName = strrep(theMRGCMosaicRFmappingResponsesFullFileName, ...
            '.mat', mappedPositionString);

theInputConeMosaicRFmappingResponsesFullFileName = strrep(theInputConeMosaicRFmappingResponsesFullFileName, ...
            '.mat', mappedPositionString);


% Crop mosaic to the mapped region
fprintf('Trimming mRGCMosaic to %2.1f x %2.1f degs at (x,y) = (%2.1f,%2.1f)\n', ...
        mappedSizeDegs(1), mappedSizeDegs(2), mappedPositionDegs(1), mappedPositionDegs(2));

% Extra support for input cone mosaic supporting the RF surrounds
extraSupportDegsForInputConeMosaic = 0.3*mRGCMosaic.extraSupportDegsForMidgetRGCSurrounds(abs(mappedPositionDegs)+mappedSizeDegs, 0.0);
visualizeSpatialRelationshipToSourceMosaic = true;

theMRGCMosaic.cropToSizeAtEccentricity(mappedSizeDegs, mappedPositionDegs, ...
        'extraSupportDegsForInputConeMosaic', extraSupportDegsForInputConeMosaic, ...
        'visualizeSpatialRelationshipToSourceMosaic', visualizeSpatialRelationshipToSourceMosaic);

theMRGCMosaic.visualize();

theTargetVisualizedRGCindices = theMRGCMosaic.indicesOfRGCsWithinROI(...
	mappedPositionDegs, [mappedSizeDegs(1) mappedSizeDegs(2)]);

% Specify smoothing kernel params 
[smoothinKernelSizePixels, smoothinKernelSigmaPixels] = testsHelperFunctions.smoothinKernelParams();

coneFundamentalAndSpatialResolutionString = '';
profileGain = 1;
zLevelsNegative = -1/profileGain*(0.002:0.05:0.5);
zLevelsPositive =  1/profileGain*(0.002:0.05:0.5);

onlyPlotPreviouslyComputedRoRincRatios = false;
[smoothinKernelSizePixels, smoothinKernelSigmaPixels] = RGCMosaicAnalyzer.visualize.mSequenceRFmaps(...
        theMRGCMosaicRFmappingResponsesFullFileName, ...
        surroundConnectedParamsStruct, ...
        theMRGCMosaic, theTargetVisualizedRGCindices, ...
        profileGain, zLevelsNegative, zLevelsPositive, ...
        smoothinKernelSizePixels, smoothinKernelSigmaPixels, ...
        coneFundamentalAndSpatialResolutionString, ...
        chromaticityForRFmapping, opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters, ...
        generateVisualRFandConePoolingMapComboPlots, ...
        onlyPlotPreviouslyComputedRoRincRatios);
    fprintf('RF map smoothing kernel employed: size = %2.4f, sigma = %2.4f\n', smoothinKernelSizePixels, smoothinKernelSigmaPixels);

