function [theFullMRGCMosaicFileName, theIntermediateConnectivityStagesMetaDataFile, theMRGCMosaicFileName, mRGCMosaicSubDir] = ...
    exportedMosaicFileName(paramsStruct, connectivityStage, varargin)

    p = inputParser;
    p.addParameter('extraSubDirPath', '', @ischar);
    % Execute the parser
    p.parse(varargin{:});
    extraSubDirPath = p.Results.extraSubDirPath;


    % Root directory
    intermediateDataDir = RGCMosaicConstructor.filepathFor.intermediateDataDir();
    if (~isempty(extraSubDirPath))
        intermediateDataDir = fullfile(intermediateDataDir, extraSubDirPath);
    end

    % Brian wants to be able to use positionDegs insead of eccentricityDegs
    if ~isempty(paramsStruct.eccentricityDegs)
        % eccentricityDegs is always used if present, for
        % historical reasons 
        paramsStruct.eccentricityDegs = paramsStruct.eccentricityDegs;
    elseif (isfield(paramsStruct, 'positionDegs')) && (~isempty(paramsStruct.positionDegs))
        % Use positionDegs if available but there is not
        % eccentricityDegs
        paramsStruct.eccentricityDegs = paramsStruct.positionDegs;
    else                
        % Default
        paramsStruct.eccentricityDegs = [0,0];
    end

    switch (paramsStruct.whichEye)
        case 'right eye'
            theEyeString = 'RE';
        case 'left eye'
            theEyeString = 'LE';
        otherwise
            error('Unknown eye: ''%s''.', paramsStruct.whichEye);
    end

    % Spatial-chromatic purity tradeoff info
    spatialChromaticUniformityTradeoffString = sprintf('Phi_%2.2f', paramsStruct.spatialChromaticUniformityTradeoff);

    switch (connectivityStage)
        case {'center connected', 'center connected with overlap'}
             if (strcmp(connectivityStage, 'center connected'))
                mRGCMosaicSubDir = 'centerConnected/';
            else
                mRGCMosaicSubDir  = 'centerConnectedOverlap/';
            end
            theMRGCMosaicFileName = sprintf('%sMRGCMosaic_%s_Ecc%2.1f_%2.1f_Size%2.1fx%2.1f_%s.mat',...
                mRGCMosaicSubDir, theEyeString, paramsStruct.eccentricityDegs(1), paramsStruct.eccentricityDegs(2), ...
                paramsStruct.sizeDegs(1), paramsStruct.sizeDegs(2), spatialChromaticUniformityTradeoffString);

        case {'surround connected', 'surround connected with center overlap'}
            % Append optics params info
            if (isempty(paramsStruct.surroundConnectivitySimulationParamsStruct.optimizationPosition))

                surroundOptimizationString = RGCMosaicConstructor.filepathFor.surroundOptimizationString(...
                    paramsStruct.surroundConnectivitySimulationParamsStruct.poolingOptimizationParamsStruct.poolingModel);

                % Filename for compute-ready mosaic
                opticsSurroundOptimizationComboString = sprintf('Optics_%s-%d_%s_%s', ...
                    paramsStruct.surroundConnectivitySimulationParamsStruct.opticsParamsStruct.ZernikeDataBase, ...
                    paramsStruct.surroundConnectivitySimulationParamsStruct.opticsParamsStruct.subjectID, ...
                    paramsStruct.surroundConnectivitySimulationParamsStruct.opticsParamsStruct.modification, ...
                    surroundOptimizationString);
            else
                % Filename for optimization results 
                opticsSurroundOptimizationComboString = sprintf('OpticsAt_%2.2f_%2.2f_%s-%d_%s', ...
                    paramsStruct.surroundConnectivitySimulationParamsStruct.optimizationPosition(1), ...
                    paramsStruct.surroundConnectivitySimulationParamsStruct.optimizationPosition(2), ...
                    paramsStruct.surroundConnectivitySimulationParamsStruct.opticsParamsStruct.ZernikeDataBase, ...
                    paramsStruct.surroundConnectivitySimulationParamsStruct.opticsParamsStruct.subjectID, ...
                    paramsStruct.surroundConnectivitySimulationParamsStruct.opticsParamsStruct.modification);
            end

            if (strcmp(connectivityStage, 'surround connected'))
                mRGCMosaicSubDir = 'surroundConnected/';
            else
                mRGCMosaicSubDir  = 'surroundConnectedOverlap/';
            end

            targetVisualSTFModifiersString = sprintf('vSTF_%1.1f_%1.1f', ...
                paramsStruct.surroundConnectivitySimulationParamsStruct.targetVisualSTFmodifiersFromMeanValuesStruct.surroundToCenterRcRatioMultiplier, ...
                paramsStruct.surroundConnectivitySimulationParamsStruct.targetVisualSTFmodifiersFromMeanValuesStruct.surroundToCenterIntegratedSensitivityRatioMultiplier);

            
            theMRGCMosaicFileName = sprintf('%sMRGCMosaic_%s_Ecc%2.1f_%2.1f_Size%2.1fx%2.1f_%s_%s_%s.mat',...
                mRGCMosaicSubDir, theEyeString, paramsStruct.eccentricityDegs(1), paramsStruct.eccentricityDegs(2), ...
                paramsStruct.sizeDegs(1), paramsStruct.sizeDegs(2), ...
                spatialChromaticUniformityTradeoffString, ...
                opticsSurroundOptimizationComboString, ...
                targetVisualSTFModifiersString);

        otherwise
            error('Unknown connectivity stage: ''%s''.', connectivityStage);
    end


    if (~isempty(paramsStruct.customLMSconeDensities))
       customLMSconePercentages = round(paramsStruct.customLMSconeDensities*100);
       customLMSconeDensityString = sprintf('MRGCMosaic_%1.0f_%1.0f_%1.0f', ...
            customLMSconePercentages(1), customLMSconePercentages(2), customLMSconePercentages(3));
       theMRGCMosaicFileName = strrep(theMRGCMosaicFileName, 'MRGCMosaic', customLMSconeDensityString);
    end

    theFullMRGCMosaicFileName = fullfile(intermediateDataDir, theMRGCMosaicFileName);

    % Synthesize file with intermediate connectivity stages meta data so we can
    % visualize later connecticity at different stages
    theIntermediateConnectivityStagesMetaDataFile = strrep(theFullMRGCMosaicFileName, mRGCMosaicSubDir, 'centerConnectivityStages/');
end