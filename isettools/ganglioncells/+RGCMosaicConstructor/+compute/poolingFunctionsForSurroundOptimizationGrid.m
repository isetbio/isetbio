function [optimizationResultsFileNames, optimizationResultsTargetCenterConesNum, theSurroundConnectedMRGCMosaicFullFileNames] = ...
    poolingFunctionsForSurroundOptimizationGrid(whichEye, eccentricityDegs, sizeDegs, ...
    spatialChromaticUniformityTradeoff, ...
    customLMSconeDensities, ...
    surroundConnectivitySimulationParamsStruct, varargin)

% Parse input
p = inputParser;

% Required params
p.addRequired('whichEye', @(x)(ischar(x) && (ismember(x, {'left eye', 'right eye'}))));
p.addRequired('eccentricityDegs', @(x)(isnumeric(x) && (numel(x) == 2)));
p.addRequired('sizeDegs', @(x)(isnumeric(x) && (numel(x) == 2)));
p.addRequired('spatialChromaticUniformityTradeoff', @(x)(isscalar(x) && (x>=0) && (x<=1)));
p.addRequired('customLMSconeDensities', @(x)(isempty(x)||((isnumeric(x))&&(numel(x)==3))));
p.addRequired('surroundConnectivitySimulationParamsStruct', @(x)rfSurroundConnectivityValidationFunction(x));

% Actions to take
p.addParameter('computeInputConeMosaicResponses', false, @islogical);
p.addParameter('optimizeSurroundConePooling', false, @islogical);
p.addParameter('doNotWorryAboutMaximizingTargetRFcenterLMconeRatio', false, @islogical);
p.addParameter('optimizationPositionIndicesToCompute', [], @isnumeric);
p.addParameter('centerConeNumerositiesToOptimize', [], @isnumeric);
p.addParameter('centerConeDominanceToOptimize', cMosaic.LCONE_ID, @(x)(ismember(x, [cMosaic.LCONE_ID cMosaic.MCONE_ID])));

% Optional params
p.addParameter('userSuppliedInitialValuesForModelVariables', [], @(x)(isempty(x)||isstruct(x)));
p.addParameter('randomSeed', [], @(x)(isempty(x) || isscalar(x)));
p.addParameter('onlyReturnSurroundOptimizationResultFilenames', false, @islogical);

p.addParameter('visualizeInputConeMosaicSTFResponseSequences', false, @islogical);
p.addParameter('visualizeFullAndMaximalExcursionSTF', false, @islogical);
p.addParameter('visualizeGaussianFitToCenterSTF', false, @islogical);

p.addParameter('visualizeSamplingPositionsForUncroppedMosaic', false, @islogical);
p.addParameter('visualizeSpatialRelationshipToSourceMosaic', false, @islogical);

p.parse(whichEye, eccentricityDegs, sizeDegs, spatialChromaticUniformityTradeoff, ...
    customLMSconeDensities, surroundConnectivitySimulationParamsStruct, varargin{:});

computeInputConeMosaicResponses = p.Results.computeInputConeMosaicResponses;
optimizeSurroundConePooling = p.Results.optimizeSurroundConePooling;
doNotWorryAboutMaximizingTargetRFcenterLMconeRatio = p.Results.doNotWorryAboutMaximizingTargetRFcenterLMconeRatio;
optimizationPositionIndicesToCompute = p.Results.optimizationPositionIndicesToCompute;
centerConeNumerositiesToOptimize = p.Results.centerConeNumerositiesToOptimize;
centerConeDominanceToOptimize = p.Results.centerConeDominanceToOptimize;

userSuppliedInitialValuesForModelVariables = p.Results.userSuppliedInitialValuesForModelVariables;
onlyReturnSurroundOptimizationResultFilenames = p.Results.onlyReturnSurroundOptimizationResultFilenames;

visualizeInputConeMosaicSTFResponseSequences = p.Results.visualizeInputConeMosaicSTFResponseSequences;
visualizeFullAndMaximalExcursionSTF = p.Results.visualizeFullAndMaximalExcursionSTF;
visualizeGaussianFitToCenterSTF = p.Results.visualizeGaussianFitToCenterSTF;
visualizeSamplingPositionsForUncroppedMosaic = p.Results.visualizeSamplingPositionsForUncroppedMosaic;
visualizeSpatialRelationshipToSourceMosaic = p.Results.visualizeSpatialRelationshipToSourceMosaic;

randomSeed = p.Results.randomSeed;

if (surroundConnectivitySimulationParamsStruct.poolingOptimizationParamsStruct.employRFCenterOverlappingMosaic)
    centerConnectivityStage = 'center connected with overlap';
    surroundConnectivityStage = 'surround connected with center overlap';
else
    centerConnectivityStage = 'center connected';
    surroundConnectivityStage = 'surround connected';
end

% Chromaticity to use
chromaParamsStruct = struct(...
    'stimulusChromaticity', surroundConnectivitySimulationParamsStruct.STFparamsStruct.chromaticity, ...
    'coneFundamentalsOptimizedForStimPosition', surroundConnectivitySimulationParamsStruct.STFparamsStruct.coneFundamentalsOptimizedForStimPosition);


% Go through each optimization position
if (isempty(optimizationPositionIndicesToCompute))
    optimizationPositionIndicesToCompute = 1:size(surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid,1);
end

optimizationResultsFileNames = {};
optimizationResultsTargetCenterConesNum = [];
optimizationPositionsComputed = 0;

for iOptimizationPos = 1:size(surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid,1)

    if (~ismember(iOptimizationPos,optimizationPositionIndicesToCompute))
        fprintf('Skipping optimization position %d\n', iOptimizationPos);
        continue;
    end

    optimizationPositionsComputed = optimizationPositionsComputed + 1;

    % Generate the filename where the surround-connected RGC mosaic will be saved.
    % This encodes the optics under which the visual STFresponses were computed
    % as well as the optimization position within the mosaic
    pSurround = p.Results;

    pSurround.surroundConnectivitySimulationParamsStruct.optimizationPosition = surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(iOptimizationPos,:);
    [theSurroundConnectedMRGCMosaicFullFileName, ~, theSurroundConnectedMRGCMosaicFileName, surroundConnectedMRGCMosaicSubDir] = ...
        RGCMosaicConstructor.filepathFor.exportedMosaicFileName(...
        pSurround, surroundConnectivityStage, ...
        'generateMissingSubDirs', true);
    'clear pSurround';

    % Generate the STFresponses file names
    [theInputConeMosaicSTFResponsesFullFileName, theMRGCMosaicSTFResponsesFullFileName] = ...
        RGCMosaicConstructor.filepathFor.stfResponsesFileName(...
            theSurroundConnectedMRGCMosaicFullFileName, surroundConnectedMRGCMosaicSubDir, chromaParamsStruct, pSurround);

    if (computeInputConeMosaicResponses)
        startTime = cputime;
        fprintf('\n*******************\n Computing input cone mosaic responses for optimization position %d of %d\n************\n', ...
            iOptimizationPos, size(surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid,1));

        % Generate the filename where the center-connected RGC mosaic is loaded from
        paramsStructCenterOnly = p.Results;
        paramsStructCenterOnly.rfCenterConnectivityParams.spatialChromaticUniformityTradeoff  = p.Results.spatialChromaticUniformityTradeoff;

        % Generate the center-only connected mRGCMosaic filename
        [theCenterConnectedMRGCMosaicFullFileName, ~, theCenterConnectedMRGCMosaicFileName] = ...
            RGCMosaicConstructor.filepathFor.exportedMosaicFileName(paramsStructCenterOnly, centerConnectivityStage);

        try 
            % Load the center-connected mosaic
            load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');
            fprintf('Succesffuly loaded MRGC mosaic from\n%s\n', theCenterConnectedMRGCMosaicFullFileName);
        catch
            fprintf('The requested MRGCMosaic file\n\t %s\nwas not found.\n', theCenterConnectedMRGCMosaicFullFileName);
            if (strcmp(centerConnectivityStage, 'center connected with overlap'))
                fprintf('Will try loading mosaic without RF center overlap now.\n');
                [theCenterConnectedMRGCMosaicFullFileName, ~, theCenterConnectedMRGCMosaicFileName] = ...
                    RGCMosaicConstructor.filepathFor.exportedMosaicFileName(paramsStructCenterOnly, 'center connected');
                load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');
            end
            fprintf('Succesffuly loaded mosaic without RF center overlap from\n%s\n', theCenterConnectedMRGCMosaicFullFileName);
        end
        clear 'paramsStructCenterOnly';

        if (theMRGCMosaic.rgcRFcentersOverlap == false)
            error('The MRGCMosaic loaded does not have RF center overlap.');
        end

           
        if (visualizeSamplingPositionsForUncroppedMosaic) && (optimizationPositionsComputed == 1)
            superimposedSamplingRectsSmall.center = surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(:,1:2);
            superimposedSamplingRectsSmall.xRange = surroundConnectivitySimulationParamsStruct.optimizationSizeGrid(:,1)*0 + 0.1;
            superimposedSamplingRectsSmall.yRange = surroundConnectivitySimulationParamsStruct.optimizationSizeGrid(:,1)*0 + 0.1;

            fprintf('Please wait. Visualizing mosaic and ALL optimization position.');
            theMRGCMosaic.visualize(...
                    'superimposedRect', superimposedSamplingRectsSmall, ...
                    'superimposedRectColor', [1 1 0], ...
                    'superimposedRectLineWidth', 0.1, ...
                    'superimposedRectAlpha', 1.0, ...
                    'plotTitle','all positions');

            superimposedRectsActual.center = surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(:,1:2);
            superimposedRectsActual.xRange = surroundConnectivitySimulationParamsStruct.optimizationSizeGrid(:,1);
            superimposedRectsActual.yRange = surroundConnectivitySimulationParamsStruct.optimizationSizeGrid(:,2);

            fprintf('Please wait. Visualizing mosaic and ALL optimization position.');
            theMRGCMosaic.visualize(...
                    'superimposedRect', superimposedRectsActual, ...
                    'superimposedRectColor', [1 0.5 0.5], ...
                    'superimposedRectLineWidth', 0.5, ...
                    'superimposedRectAlpha', 0.2, ...
                    'plotTitle','all positions (actual sizes)');

            for iiOptimizationPos = 1:size(surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid,1)
                superimposedRect.center = surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(iiOptimizationPos,:);
                superimposedRect.xRange = surroundConnectivitySimulationParamsStruct.optimizationSizeGrid(iiOptimizationPos,1);
                superimposedRect.yRange = surroundConnectivitySimulationParamsStruct.optimizationSizeGrid(iiOptimizationPos,2);

                fprintf('Please wait. Visualizing mosaic and current optimization position.');
                theMRGCMosaic.visualize(...
                    'superimposedRect', superimposedRect, ...
                    'superimposedRectColor', [1 1 0], ...
                    'superimposedRectLineWidth', 1.0, ...
                    'superimposedRectAlpha', 0.4, ...
                    'plotTitle', sprintf('optimization position: %d of %d', iiOptimizationPos, size(surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid,1)));
            end % for iiOptimizationPos
        end % if (visualizeSamplingPositions)

        % Compute extra support for input cone mosaic at the optimization position
        extraSupportDegsForInputConeMosaic = mRGCMosaic.extraSupportDegsForMidgetRGCSurrounds(...
            surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(iOptimizationPos,:), ...
            surroundConnectivitySimulationParamsStruct.optimizationSizeGrid(iOptimizationPos,:));

        % Crop the mosaic at the current optimization position
        theMRGCMosaic.cropToSizeAtEccentricity(...
            surroundConnectivitySimulationParamsStruct.optimizationSizeGrid(iOptimizationPos,:), ...
            surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(iOptimizationPos,:), ...
            'extraSupportDegsForInputConeMosaic', extraSupportDegsForInputConeMosaic, ...
            'visualizeSpatialRelationshipToSourceMosaic', visualizeSpatialRelationshipToSourceMosaic);

        % STEP 1: Generate optics at optimizationPositionDegs
        [theOI, thePSF] = RGCMosaicConstructor.helper.optics.generate(...
            theMRGCMosaic.inputConeMosaic, ...
            surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(iOptimizationPos,:), ...
            surroundConnectivitySimulationParamsStruct.opticsParamsStruct, ...
            'visualizePSF', true);

        % STEP 2: Compute input cone mosaic responses to gratings of different orientations/SFs
        % Determine the stimulus pixel resolution to be a fraction of the minimum cone aperture or cone spacing in the mosaic
        % here, half of the cone spacing
        theMetric = 'cone aperture';  % choose from {'cone aperture' or cone spacing'}
        theFraction = 0.25;

        targetRGCindices =  1:theMRGCMosaic.rgcsNum;
        stimulusResolutionDegs = RGCMosaicConstructor.helper.simulateExperiment.stimulusResolutionFromConeApertureOrConeSpacing(...
            theMRGCMosaic, targetRGCindices, theFraction, theMetric);

        % Update the STF stimulus params struct
        surroundConnectivitySimulationParamsStruct.STFparamsStruct.resolutionDegs = stimulusResolutionDegs;
        surroundConnectivitySimulationParamsStruct.STFparamsStruct.positionDegs = theMRGCMosaic.eccentricityDegs;
        surroundConnectivitySimulationParamsStruct.STFparamsStruct.sizeDegs =  1.1*max(theMRGCMosaic.inputConeMosaic.sizeDegs);

        % Compute the STF responses and save them for later processing
        RGCMosaicConstructor.helper.simulateExperiment.spatialTransferFunction(...
            theMRGCMosaic, theOI, surroundConnectivitySimulationParamsStruct.STFparamsStruct, ...
            theInputConeMosaicSTFResponsesFullFileName, ...
            theMRGCMosaicSTFResponsesFullFileName, ...
            'computeInputConeMosaicResponses', computeInputConeMosaicResponses, ...
            'computeMRGCMosaicResponses', false, ...
            'visualizeResponse', visualizeInputConeMosaicSTFResponseSequences);

        endTime = cputime;
        fprintf('\n*******************\n Finished computing input cone mosaic responses for optimization position %d of %d in %2.1f minutes.\n************\n', ...
            iOptimizationPos, size(surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid,1), (endTime-startTime)/60);
    end %if (computeInputConeMosaicResponses)

    if (optimizeSurroundConePooling) || (onlyReturnSurroundOptimizationResultFilenames)
        % Size of optimization
        optimizationRegionString = sprintf('optimizationSizeDegs_%2.1f_%2.1f', ...
          surroundConnectivitySimulationParamsStruct.optimizationSizeGrid(iOptimizationPos,1), ...
          surroundConnectivitySimulationParamsStruct.optimizationSizeGrid(iOptimizationPos,2));

        % Assemble the surroundOptimizationString (to be added to theSurroundConnectedMRGCMosaicFullFileNames )
        surroundOptimizationString = RGCMosaicConstructor.filepathFor.surroundOptimizationString(...
            surroundConnectivitySimulationParamsStruct.poolingOptimizationParamsStruct.poolingModel);

        % Add the optimization region
        surroundOptimizationString = sprintf('_%s_%s.mat', surroundOptimizationString, optimizationRegionString);

        % Generate theSurroundConnectedMRGCMosaicFullFileNames
        theSurroundConnectedMRGCMosaicFullFileNames{iOptimizationPos} = ...
            strrep(theSurroundConnectedMRGCMosaicFullFileName, '.mat', surroundOptimizationString);
        theSurroundConnectedMRGCMosaicFileNames{iOptimizationPos} = ...
            strrep(theSurroundConnectedMRGCMosaicFileName, '.mat', surroundOptimizationString);

        % Compute components for the surround pooling optimization
        dataOut = RGCMosaicConstructor.compute.componentsForRFsurroundConnectivity(...
            theInputConeMosaicSTFResponsesFullFileName, ...
            surroundConnectivitySimulationParamsStruct.targetVisualSTFmodifiersFromMeanValuesStruct, ...
            centerConeNumerositiesToOptimize, centerConeDominanceToOptimize, ...
            surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(iOptimizationPos,:), ...
            doNotWorryAboutMaximizingTargetRFcenterLMconeRatio, ...
            ~onlyReturnSurroundOptimizationResultFilenames);

        if (onlyReturnSurroundOptimizationResultFilenames)
            for iCenterConesNum = 1:numel(dataOut.uniqueCenterConesNum)
                theOptimizedCenterConesNum = dataOut.uniqueCenterConesNum(iCenterConesNum);

                if (isempty(dataOut.centerConeNumerositiesToOptimize)) || (ismember(theOptimizedCenterConesNum, dataOut.centerConeNumerositiesToOptimize))
                    
                    optimizationResultsTargetCenterConesNum(numel(optimizationResultsTargetCenterConesNum)+1) = theOptimizedCenterConesNum;
                    optimizationResultsFileNames{numel(optimizationResultsFileNames)+1} = RGCMosaicConstructor.filepathFor.optimizationResults(...
                            theSurroundConnectedMRGCMosaicFullFileNames{iOptimizationPos}, ...
                            theOptimizedCenterConesNum, ...
                            centerConeDominanceToOptimize);
                end % if (ismember(theOptimizedCenterConesNum, dataOut.RFcenterConesNumToOptimizeFor))
            end % for iCenterConesNum
            fprintf('Only returning optimizationResults filename. \n')
            continue;
        end

        % STEP 3. Optimize weights of pooled cone mosaic responses to yield desired STF
        theMRGCMosaic = dataOut.theMRGCMosaic;
        targetVisualSTFparams = dataOut.targetVisualSTFparams;
        stimParams = dataOut.stimParams;
        theInputConeMosaicSTFresponses = dataOut.theInputConeMosaicSTFresponses;


        % Check whether we are using RF center overlap
        if ((surroundConnectivitySimulationParamsStruct.poolingOptimizationParamsStruct.employRFCenterOverlappingMosaic) && ...
            (theMRGCMosaic.rgcRFcentersOverlap == false))
            error('Specified to use mosaic with RF center overlap, but the RFcenter overlap has not been performed on the loaded mosaic.\n');
        end

        for iCenterConesNum = 1:numel(dataOut.uniqueCenterConesNum)
            theOptimizedCenterConesNum = dataOut.uniqueCenterConesNum(iCenterConesNum);

            if (isempty(dataOut.centerConeNumerositiesToOptimize)) || (ismember(theOptimizedCenterConesNum, dataOut.centerConeNumerositiesToOptimize))
                theOptimizedPositionDegs = stimParams.positionDegs;

                % Optimize surround pooling for theTargetRGCindex
                theTargetRGCindex = dataOut.theTargetRGCindex(iCenterConesNum);

                fprintf('----> Optimizing surround cone pooling for %d-center cone(s) RF center for mosaic centered at (%2.2f,%2.2f). Target RGC is located at (%2.2f,%2.2f)\n', ...
                    theOptimizedCenterConesNum, theOptimizedPositionDegs(1), theOptimizedPositionDegs(2), ...
                    theMRGCMosaic.rgcRFpositionsDegs(theTargetRGCindex,1), theMRGCMosaic.rgcRFpositionsDegs(theTargetRGCindex,2));

                minConeWeightVisualized = mRGCMosaic.sensitivityAtPointOfOverlap;
                connectivityVector = full(squeeze(theMRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,theTargetRGCindex)));
                identifiedInputConeIndices = find(connectivityVector>= minConeWeightVisualized);
                allInputConeIndices = find(connectivityVector >=  exp(-4));
                
                hFig = figure(667);
                ax = subplot(1,1,1);

                XLims(1) = min(squeeze(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(allInputConeIndices,1)))-0.1;
                XLims(2) = max(squeeze(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(allInputConeIndices,1)))+0.1;
                YLims(1) = min(squeeze(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(allInputConeIndices,2)))-0.1;
                YLims(2) = max(squeeze(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(allInputConeIndices,2)))+0.1;

                domainVisualizationLimits = [XLims(1) XLims(2) YLims(1) YLims(2) ];
                domainVisualizationTicks = struct('x', -40:0.5:0, 'y', -10:0.5:10);

                theMRGCMosaic.visualize(...
                    'figureHandle', hFig, ...
                    'axesHandle', ax, ...
                    'visualizedRGCindices', theTargetRGCindex, ...
                    'pooledConesLineWidth', 1.5, ...
                    'plottedRFoutlineLineWidth', 1.0, ...
                    'plottedRFoutlineFaceAlpha', 0.5, ...
                    'identifyPooledCones', true, ...
                    'identifyInputCones', true, ...
                    'inputConesAlpha', 1.0, ...
                    'identifiedInputConeIndices', identifiedInputConeIndices, ... % which cones to identity with their type
                    'minConeWeightVisualized', minConeWeightVisualized, ...       % where to draw the center profile
                    'centerSubregionContourSamples', 40, ...
                    'identifiedConeApertureThetaSamples', 40, ...
                    'identifiedConeAperture', 'lightCollectingArea6sigma', ...
                    'domainVisualizationLimits', domainVisualizationLimits, ...
                    'domainVisualizationTicks', domainVisualizationTicks, ...
                    'plotTitle', sprintf('RGCindex being fitted (%d)', theTargetRGCindex));

                if (isstruct(userSuppliedInitialValuesForModelVariables)) &&...
                        (isfield(userSuppliedInitialValuesForModelVariables, 'attemptToLoadExactOptimizationResultsFile')) && ...
                        (userSuppliedInitialValuesForModelVariables.attemptToLoadExactOptimizationResultsFile)

                        if (userSuppliedInitialValuesForModelVariables.skipOptimizationIfOptimizationFileExists)
                            optimizationResults = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.importOptimizationResults(...
                            theSurroundConnectedMRGCMosaicFullFileNames{iOptimizationPos}, ...
                            theOptimizedCenterConesNum, centerConeDominanceToOptimize, 'imported exact match', ...
                            'tryAlternateH1CellIndexAsLastResort', false);

                            if (~isempty(optimizationResults))
                                fprintf('>>>>> Found previous optimization file, and since the ''skipOptimizationIfOptimizationFileExists'' is set to true, we are skipping re-optimization !! <<<< \n');
                                continue;
                            end
                        end

                        if (strcmp(surroundConnectivitySimulationParamsStruct.poolingOptimizationParamsStruct.poolingModel.name, 'PackerDacey2002H1FixedCellIndex'))
                            tryAlternateH1CellIndexAsLastResort = true;
                        else
                            tryAlternateH1CellIndexAsLastResort = false;
                        end
                        optimizationResults = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.importOptimizationResults(...
                            theSurroundConnectedMRGCMosaicFullFileNames{iOptimizationPos}, ...
                            theOptimizedCenterConesNum, centerConeDominanceToOptimize, 'imported closest match', ...
                            'tryAlternateH1CellIndexAsLastResort', tryAlternateH1CellIndexAsLastResort);

                        if (~isempty(optimizationResults))
                            userSuppliedInitialValuesForModelVariables.optimizedValues = ...      
                                RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.extractPreviouslyOptimizedValues(optimizationResults);
                        end
                end % if (userSuppliedInitialValuesForModelVariables.attemptToLoadExactOptimizationResultsFile)
               
                % Optimize away !!
                optimizationResults = RGCMosaicConstructor.compute.optimizedSurroundPoolingForTargetRGC(...
                    theMRGCMosaic, theTargetRGCindex, ...
                    targetVisualSTFparams, surroundConnectivitySimulationParamsStruct, stimParams, ...
                    theInputConeMosaicSTFresponses, ...
                    'userSuppliedInitialValuesForModelVariables', userSuppliedInitialValuesForModelVariables, ...
                    'visualizeFullAndMaximalExcursionSTF', visualizeFullAndMaximalExcursionSTF, ...
                    'visualizeGaussianFitToCenterSTF', visualizeGaussianFitToCenterSTF);

                % Export optimization results
                RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.exportOptimizationResults(...
                    theSurroundConnectedMRGCMosaicFullFileNames{iOptimizationPos}, ...
                    theMRGCMosaic, theTargetRGCindex, targetVisualSTFparams, optimizationResults);
            else
                fprintf('Skipping optimization for %d RF centers\n', theOptimizedCenterConesNum);
            end % if (ismember(theOptimizedCenterConesNum, dataOut.RFcenterConesNumToOptimizeFor))
        end % iCenterConesNum

    end % if (optimizeSurroundConePooling)
end % iOptimizationPos
end

function s = rfSurroundConnectivityValidationFunction(x)
c(1) = isstruct(x);
c(numel(c)+1) = isfield(x, 'optimizationPositionsGrid');
c(numel(c)+1) = isfield(x, 'optimizationSizeGrid');
c(numel(c)+1) = isfield(x, 'opticsParamsStruct');
c(numel(c)+1) = isfield(x, 'STFparamsStruct');
c(numel(c)+1) = isfield(x, 'fitParams');
c(numel(c)+1) = isfield(x, 'poolingOptimizationParamsStruct');
c(numel(c)+1) = isfield(x, 'targetVisualSTFmodifiersFromMeanValuesStruct');

if (~all(c))
    x
    error('Invalid rfSurroundConnectivityStruct');
else
    s = true;
end

assert(size(x.optimizationPositionsGrid,2)==2, 'optimizationPositionsGrid must have 2 columns');
assert(size(x.optimizationSizeGrid,2)==2, 'optimizationSizeGrid must have 2 columns');
assert(size(x.optimizationPositionsGrid,1) == size(x.optimizationSizeGrid,1), 'optimizationPositionsGrid and optimizationSizeGrid must have the same number of rows');
end