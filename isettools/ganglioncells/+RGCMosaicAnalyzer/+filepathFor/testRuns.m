function [theSurroundConnectedMRGCMosaicFullFileName, theInputConeMosaicResponsesFullFileName, ...
    theMRGCMosaicResponsesFullFileName, surroundConnectedParamsStruct, targetVisualSTFmodifierStruct] = ...
        testRuns(whichEye, whichZernikeDataBase, whichSubjectID, ...
            mosaicEccDegs, mosaicSizeDegs, employRFCenterOverlappingMosaic, ...
            spatialChromaticUniformityTradeoff, customLMSconeDensities, ...
            targetVisualSTFdescriptor, surroundOptimizationStrategy, ...
            opticsModification, residualWithRespectToNativeOpticsDefocusDiopters, ...
            STFchromaticity, coneFundamentalsOptimizedForStimPosition, ...
            responsesSubDir)

    % Generate the surroundRetinalConePoolingModel params struct
    surroundRetinalConePoolingModelParamsStruct = ...
        RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateSurroundRetinalConePoolingStruct(surroundOptimizationStrategy);
    if (strcmp(surroundRetinalConePoolingModelParamsStruct.name, 'PackerDacey2002H1FixedCellIndex')) 
        theFixedH1CellIndex = RGCMosaicConstructor.helper.queryUserFor.fixedH1CellIndex(false, []);
        % Add to the surroundRetinalConePoolingModelParamsStruct  the fixed H1 cell index to be used
        surroundRetinalConePoolingModelParamsStruct.fixedH1CellIndex = theFixedH1CellIndex;
    end

    % Generate targetVisualSTFmodifierStruct
    targetVisualSTFmodifierStruct = ...
        RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateTargetVisualSTFmodifiersStruct(targetVisualSTFdescriptor);

    depictTargetVisualSTFinRelationshipToCronerKaplanData = false;
    if (depictTargetVisualSTFinRelationshipToCronerKaplanData)
        RGCMosaicConstructor.visualize.CronerKaplanRawData(...
            'surroundToCenterRcRatioMultiplier', targetVisualSTFmodifierStruct.surroundToCenterRcRatioMultiplier, ...
            'surroundToCenterIntegratedSensitivityRatioMultiplier', targetVisualSTFmodifierStruct.surroundToCenterIntegratedSensitivityRatioMultiplier);
    end

    % Generate the surroundConnectivity simulation params struct
    optimizationPositionsAndSizesGridsDummy = [0 0 0 0];
    surroundConnectivitySimulationParamsStruct = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateSurroundConnectivitySimulationParamsStruct(...
        whichZernikeDataBase, whichSubjectID, employRFCenterOverlappingMosaic, ...
        optimizationPositionsAndSizesGridsDummy, surroundRetinalConePoolingModelParamsStruct, ...
        targetVisualSTFmodifierStruct);

    % Generate filename for the surround-connected, compute-ready MRGCmosaic
    surroundConnectedParamsStruct.whichEye = whichEye;
    surroundConnectedParamsStruct.eccentricityDegs = mosaicEccDegs;
    surroundConnectedParamsStruct.sizeDegs = mosaicSizeDegs;
    surroundConnectedParamsStruct.customLMSconeDensities = customLMSconeDensities;
    surroundConnectedParamsStruct.spatialChromaticUniformityTradeoff = spatialChromaticUniformityTradeoff;
    surroundConnectedParamsStruct.surroundConnectivitySimulationParamsStruct = surroundConnectivitySimulationParamsStruct;
    surroundConnectedParamsStruct.surroundConnectivitySimulationParamsStruct.optimizationPosition = [];
    if (employRFCenterOverlappingMosaic)
        surroundConnectivityStage = 'surround connected with center overlap';
    else
        surroundConnectivityStage = 'surround connected';
    end

    % The surround-connected, compute-reday MRGCmosaic
    theSurroundConnectedMRGCMosaicFullFileName = RGCMosaicConstructor.filepathFor.exportedMosaicFileName(...
                    surroundConnectedParamsStruct, surroundConnectivityStage);

    % Assemble filenames for input cone mosaic and mRGCmosaic, adding optics and chromaticity info to the filename
    theInputConeMosaicResponsesFullFileName = strrep(theSurroundConnectedMRGCMosaicFullFileName, 'surroundConnectedOverlap', responsesSubDir);
    

    if (coneFundamentalsOptimizedForStimPosition)
        STFchromaticity = sprintf('%sPosOptim',STFchromaticity);
    end

    if (strcmp(opticsModification, 'customRefraction'))
        theInputConeMosaicResponsesFullFileName = strrep(theInputConeMosaicResponsesFullFileName, '.mat', sprintf('_inputConeMosaic_%s_%s.mat',sprintf('%s_%2.2fD', opticsModification, customRefractionDiopters), STFchromaticity));

    elseif (strcmp(opticsModification, 'refractionResidualWithRespectToNativeOptics'))
        if (residualWithRespectToNativeOpticsDefocusDiopters>0)
            theInputConeMosaicResponsesFullFileName = strrep(theInputConeMosaicResponsesFullFileName, '.mat', sprintf('_inputConeMosaic_%s_%s.mat',sprintf('defocused+%2.2fDdef', residualWithRespectToNativeOpticsDefocusDiopters), STFchromaticity));
        else
            theInputConeMosaicResponsesFullFileName = strrep(theInputConeMosaicResponsesFullFileName, '.mat', sprintf('_inputConeMosaic_%s_%s.mat',sprintf('defocused-%2.2fDdef', -residualWithRespectToNativeOpticsDefocusDiopters), STFchromaticity));
        end
    else
        theInputConeMosaicResponsesFullFileName = strrep(theInputConeMosaicResponsesFullFileName, '.mat', sprintf('_inputConeMosaic_%s_%s.mat', opticsModification, STFchromaticity));
    end
    theMRGCMosaicResponsesFullFileName = strrep(theInputConeMosaicResponsesFullFileName, 'inputConeMosaic', 'mRGCMosaic');
end