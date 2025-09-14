%% Generate mRGC mosaics at different substages of stage 3 synthesis (cone-RFsurround connectivity)
%
% Description:
%   Demonstrates how to generate an mRGC mosaic at stage 3A, 3B or 3C of connectivity
%   At stage 3A we compute input cone mosaic STF responses
%   At stage 3B ...
%

% History:
%    08/28/25  NPC  Wrote it.


close all; clear all;

% Specify desired spatialChromaticUniformityTradeoff for the RGC RF centers
% A value of 1 corresponds to maximal spatial homogeneity
% A value of 0 corresponds to maximal spectral purity
spatialChromaticUniformityTradeoff = 1.0; 

% Name encoding rgcMosaic
rgcMosaicName = 'PLOSpaperNasal7DegsMosaic';

% Which optics to employ
opticsSubjectName = 'PLOSpaperDefaultSubject';

% Which species to employ 
% Choose between {'macaque', 'human'}. If 'macaque' is chosen, the input
% cone mosaic has a 1:1 L/M cone ratio.
coneMosaicSpecies = 'human';


% Default STF parameters: mean Rs/Rc, mean Ks/Kc (Rs/Rc)^2
targetVisualSTFdescriptorToOptimizeFor = 'default';

% STF with an Rs/Rc ratio that is 1.3 x mean, and mean Ks/Kc (Rs/Rc)^2
%targetVisualSTFdescriptorToOptimizeFor = 'x1.3 RsRcRatio';

% Generate the necessary mosaic params struct
pStruct = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
    coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptorToOptimizeFor);

% Whether to regenerate the mosaic at stage3A 
% (computation of input cone mosaic STF responses)
regenerateMosaicAtStage3A = true;

% Whether to visualize the input cone mosaic responses during their computation
visualizeInputConeMosaicSTFResponseSequences = ~true;

% Grid of (X,Y)-positions, (W,H)-sizes on which the surround will be optimized 
optimizationPositionsAndSizesGrids = RGCMosaicConstructor.compute.surroundOptimizationGrid(...
		pStruct.rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme, ...
		pStruct.rgcMosaicSurroundOptimization.minGridSize, ...
		pStruct.rgcMosaicSurroundOptimization.maxGridSize, ...
		pStruct.whichZernikeDataBase, pStruct.whichEye, pStruct.sourceLatticeSizeDegs, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
		'withExtremePositions', pStruct.rgcMosaicSurroundOptimization.addEightExtremePositions);

optimizationPositionIndicesToCompute = 1:size(optimizationPositionsAndSizesGrids,1);

% Generate the surroundRetinalConePoolingModel params struct
surroundRetinalConePoolingModelParamsStruct = ...
	RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateSurroundRetinalConePoolingStruct(...
		pStruct.rgcMosaicSurroundOptimization.optimizationStrategy)

% Generate targetVisualSTFmodifierStruct
targetVisualSTFmodifierStruct = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateTargetVisualSTFmodifiersStruct(...
	pStruct.rgcMosaicSurroundOptimization.targetVisualSTFdescriptor)

% Generate the surroundConnectivity simulation params struct 
surroundConnectivitySimulationParamsStruct = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateSurroundConnectivitySimulationParamsStruct(...
	pStruct.whichZernikeDataBase, pStruct.whichSubjectID, pStruct.rgcMosaic.employRFCenterOverlappingMosaic, ...
	optimizationPositionsAndSizesGrids, surroundRetinalConePoolingModelParamsStruct, ...
	targetVisualSTFmodifierStruct)


if (regenerateMosaicAtStage3A)
    visualizeSamplingPositionsForUncroppedMosaic = true;
	visualizeSpatialRelationshipToSourceMosaic = true;

	% Compute input cone mosaic STF responses
	RGCMosaicConstructor.compute.poolingFunctionsForSurroundOptimizationGrid(...
 		pStruct.whichEye, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
 		pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
 		pStruct.customLMSconeDensities, ...
 		surroundConnectivitySimulationParamsStruct, ...
 		'optimizationPositionIndicesToCompute', optimizationPositionIndicesToCompute, ...
 		'computeInputConeMosaicResponses', true, ...
		'optimizeSurroundConePooling', false, ...
		'visualizeInputConeMosaicSTFResponseSequences', visualizeInputConeMosaicSTFResponseSequences, ...
		'visualizeSamplingPositionsForUncroppedMosaic', visualizeSamplingPositionsForUncroppedMosaic, ...
    	'visualizeSpatialRelationshipToSourceMosaic', visualizeSpatialRelationshipToSourceMosaic);
    return;
end  % regenerateMosaicAtStage3A
