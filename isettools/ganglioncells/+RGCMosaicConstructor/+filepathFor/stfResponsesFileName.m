function [theInputConeMosaicSTFResponsesFullFileName, theMRGCMosaicSTFResponsesFullFileName] = stfResponsesFileName(...
        theSurroundConnectedMRGCMosaicFullFileName, surroundConnectedMRGCMosaicSubDir, chromaParamsStruct, paramsStruct, varargin)

    p = inputParser;
    p.addParameter('extraSubDirPath', '', @ischar);
    p.addParameter('generateMissingSubDirs', false, @islogical);
    % Execute the parser
    p.parse(varargin{:});
    extraSubDirPath = p.Results.extraSubDirPath;
    generateMissingSubDirs = p.Results.generateMissingSubDirs;

    % Root directory
    intermediateDataDir = RGCMosaicConstructor.filepathFor.intermediateDataDir();
    if (~isempty(extraSubDirPath))
        intermediateDataDir = fullfile(intermediateDataDir, extraSubDirPath);
    end

    % Assemble STFresponses filename from theSurroundConnectedMRGCMosaicFullFileName
    theInputConeMosaicSTFResponsesFullFileName = ...
        strrep(theSurroundConnectedMRGCMosaicFullFileName, surroundConnectedMRGCMosaicSubDir, 'inputConeMosaicSTFresponses/');

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

     
    theInputConeMosaicSTFResponsesFullFileName = strrep(theInputConeMosaicSTFResponsesFullFileName, intermediateDataDir , '');
    theInputConeMosaicSTFResponsesFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
        intermediateDataDir , theInputConeMosaicSTFResponsesFullFileName, ...
        'generateMissingSubDirs', generateMissingSubDirs);

    theMRGCMosaicSTFResponsesFullFileName = strrep(theMRGCMosaicSTFResponsesFullFileName, intermediateDataDir,  '');
    theMRGCMosaicSTFResponsesFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
        intermediateDataDir , theMRGCMosaicSTFResponsesFullFileName, ...
        'generateMissingSubDirs', generateMissingSubDirs);
end