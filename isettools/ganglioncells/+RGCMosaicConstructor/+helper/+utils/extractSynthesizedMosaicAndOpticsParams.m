%
% RGCMosaicConstructor.helper.utils.extractSynthesizedMosaicAndOpticsParams
%
function [mosaicParams, opticsParams] = extractSynthesizedMosaicAndOpticsParams(...
    pStruct, targetVisualSTFdescriptor)

    surroundRetinalConePoolingModelParamsStruct = ...
	    RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateSurroundRetinalConePoolingStruct(...
		    pStruct.rgcMosaicSurroundOptimization.optimizationStrategy);
    
    targetVSTparams = ...
        RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateTargetVisualSTFmodifiersStruct(targetVisualSTFdescriptor);
    
    % Synthesize mosaicParams
    mosaicParams.sizeDegs = pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs;
    mosaicParams.eccDegs  = pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs;
    mosaicParams.spatialCompactnessSpectralPurityTradeoff = pStruct.rgcMosaic.spatialChromaticUniformityTradeoff;
    mosaicParams.surroundOptimizationSubString = sprintf('%s%s_vSTF_%2.1f_%2.1f', ...
        surroundRetinalConePoolingModelParamsStruct.name, ...
        surroundRetinalConePoolingModelParamsStruct.identifierString, ...
        targetVSTparams.surroundToCenterRcRatioMultiplier, ...
        targetVSTparams.surroundToCenterIntegratedSensitivityRatioMultiplier);

    % Synthesize opticsParams
    switch (pStruct.whichZernikeDataBase)
        case 'Polans2015'
            theSubjectRankings = PolansOptics.constants.subjectRanking;
        otherwise
            error('%s data based - Not yet coded.', pStruct.whichZernikeDataBase)
    end
    opticsParams.ZernikeDataBase = pStruct.whichZernikeDataBase;
    opticsParams.subjectRankOrder = theSubjectRankings(pStruct.whichSubjectID);
    opticsParams.type = 'nativeOptics';
    opticsParams.refractiveErrorDiopters = [];
    opticsParams.visualizePSFonTopOfConeMosaic = ~true;

end