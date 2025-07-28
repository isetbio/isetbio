function [theInputConeMosaicSTFResponsesFullFileName, theMRGCMosaicSTFResponsesFullFileName] = stfResponsesFileName(...
        theSurroundConnectedMRGCMosaicFullFileName, surroundConnectedMRGCMosaicSubDir, chromaParamsStruct, paramsStruct)

    % Assemble STFresponses filename: add sub-directory- inputConeMosaicSTFresponses
    theInputConeMosaicSTFResponsesFullFileName = strrep(theSurroundConnectedMRGCMosaicFullFileName, 'SLIM', 'SLIM/inputConeMosaicSTFresponses');

    % Remove the 'surroundConnected' prefix inherited from theSurroundConnectedMRGCMosaicFullFileName
    theInputConeMosaicSTFResponsesFullFileName = strrep(theInputConeMosaicSTFResponsesFullFileName, surroundConnectedMRGCMosaicSubDir, '');

    % Update STFresponses filename to include the chromaticity info
    [~, ~, theInputConeMosaicSTFResponsesFullFileName] = ...
        RGCMosaicConstructor.helper.queryUserFor.stimulusChromaticity(theInputConeMosaicSTFResponsesFullFileName, ...
        'doNotQueryUserInsteadUseThisChromaParams', chromaParamsStruct);

    % Generate corresponding mRGCMosaic STF responses filaname
    theMRGCMosaicSTFResponsesFullFileName = strrep(theInputConeMosaicSTFResponsesFullFileName, 'inputConeMosaicSTFresponses', 'STFresponses');

	targetVisualSTFModifiersString = sprintf('_vSTF_%1.1f_%1.1f', ...
         paramsStruct.surroundConnectivitySimulationParamsStruct.targetVisualSTFmodifiersFromMeanValuesStruct.surroundToCenterRcRatioMultiplier, ...
         paramsStruct.surroundConnectivitySimulationParamsStruct.targetVisualSTFmodifiersFromMeanValuesStruct.surroundToCenterIntegratedSensitivityRatioMultiplier);

	theInputConeMosaicSTFResponsesFullFileName = strrep(theInputConeMosaicSTFResponsesFullFileName, ...
		targetVisualSTFModifiersString, '');
end