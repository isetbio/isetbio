%% Generate mRGC mosaics at different substages of stage 2 synthesis (cone-RFcenter connectivity)
%
% Description:
%   Demonstrates how to generate an mRGC mosaic at stage 2A or 2C of connectivity
%   At stage 2A cones are connected to RGC RF centers in a
%   mutually-exclusive way (no RF center overlap). It is at this stage that
%   the user can specity the desired spatialChromaticUniformityTradeoff
%   At stage 2C cone connections diverge to nearby RGC RF centers
%   generating overlap between neighboring RF centers
%

% History:
%    08/28/25  NPC  Wrote it.


close all; clear all;

% Specify desired spatialChromaticUniformityTradeoff for the RGC RF centers
% A value of 1 corresponds to maximal spatial homogeneity
% A value of 0 corresponds to maximal spectral purity
spatialChromaticUniformityTradeoff = 1.0; 

% Name encoding rgcMosaic
rgcMosaicName = 'PLOSpaperNasal10DegsMosaic';
rgcMosaicName = 'PLOSpaperNasal14DegsMosaic';
rgcMosaicName = 'PLOSpaperNasal19DegsMosaic';
rgcMosaicName = 'PLOSpaperNasal25DegsMosaic';
rgcMosaicName = 'PLOSpaperNasal7DegsMosaic';

% Which optics to employ
opticsSubjectName = 'PLOSpaperDefaultSubject';

% Which species to employ 
% Choose between {'macaque', 'human'}. If 'macaque' is chosen, the input
% cone mosaic has a 1:1 L/M cone ratio.
coneMosaicSpecies = 'human';

% Cone types pooled by the RF center
coneTypesPooledByRFcenter = [cMosaic.LCONE_ID cMosaic.MCONE_ID];

% Generate the necessary mosaic params struct
pStruct = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
    coneMosaicSpecies, opticsSubjectName, rgcMosaicName, '');


% Whether to regenerate the mosaic at stage2A
regenerateMosaicAtStage2A = true;

% Whether to regenerate the mosaic at stage2C
regenerateMosaicAtStage2C = ~true;

% Whether to visualize the mosaic at stage2A
visualizeMosaicAtStage2A = true;

% Whether to visualize the mosaic at stage2C
visualizeMosaicAtStage2C = ~true;


if (regenerateMosaicAtStage2A)
    % Set the random seed so we get same cone assignments in the inputConeMosaic
    randomSeedForFreezingCmosaicConeAssignment = 1;
    
    % Optimization center: detemines the order in which RFs are optimized
    optimizationCenter = 'localSpacing'; % 'origin'; % localSpacing'; 'localConeToRGCdensityRatio'

    % Determine the max number of cones that can be
    % transferred or swapped between nearby RF centers, depending on
    % the mosaic eccentricity
    [maxConeInputsPerRGCToConsiderTransferToNearbyRGCs, ...
     maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs] = ...
        MosaicConnector.coneTransferAndSwapLimits(pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs);
    
    % Generate the RF center overlap - specific params
    rfCenterConnectivityParams = struct(...
            'optimizationCenter', optimizationCenter, ...
            'spatialChromaticUniformityTradeoff',  spatialChromaticUniformityTradeoff, ...
            'coneTypesToBeConnected', coneTypesPooledByRFcenter, ...
            'maxNeighborNormDistance', MosaicConnector.maxNeighborNormDistance, ...
            'maxNeighborsNum', MosaicConnector.maxNeighborsNum, ...
            'localSpacingFromCurrentCentroids', true, ...
            'maxConeInputsPerRGCToConsiderTransferToNearbyRGCs', maxConeInputsPerRGCToConsiderTransferToNearbyRGCs, ...
            'maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs', maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs, ...
            'maxPassesNum', 10, ...
            'visualizeIntermediateStagesOfCenterConnectivityOptimization', ~true, ...
            'saveIntermediateStagesOfCenterConnectivityOptimization', true);
    
    % Generate the center-connected (no RF center ovelap) mRGCMosaic
    [theCenterConnectedMRGCMosaicFullFileName, ...
     theIntermediateConnectivityStageMetaDataFile, ...
     theCenterConnectedMRGCMosaicFileName] = RGCMosaicConstructor.compute.centerConnectedMosaic(...
        pStruct.whichEye, ...
        pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
        pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
        rfCenterConnectivityParams,  ...
        pStruct.sourceLatticeSizeDegs, ...
        'randomSeed', randomSeedForFreezingCmosaicConeAssignment, ...
        'customLMSconeDensities', pStruct.customLMSconeDensities);
end


if (regenerateMosaicAtStage2C)
    % Load a previously-generated mosaic at stage 2A
    centerConnectedParamsStruct = struct(...
        'whichEye', pStruct.whichEye, ...
        'eccentricityDegs', pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
        'sizeDegs', pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
        'spatialChromaticUniformityTradeoff', spatialChromaticUniformityTradeoff, ...
        'customLMSconeDensities', pStruct.customLMSconeDensities);
    
    [theCenterConnectedMRGCMosaicFullFileName, ~, ...
     theCenterConnectedMRGCMosaicFileName] = ...
        RGCMosaicConstructor.filepathFor.exportedMosaicFileName(...
                centerConnectedParamsStruct, 'center connected');

    fprintf('Loading center-connected (stage 2A) MRGC mosaic from\n\t%s\n', theCenterConnectedMRGCMosaicFullFileName);
    % Load the center-connected mRGCMosaic
    load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');

    % RF center overlap params
    rfCenterOverlapParams = struct(...
        'sensitivityAtPointOfOverlap', mRGCMosaic.sensitivityAtPointOfOverlap, ...  
        'minSensitivityForInclusion', mRGCMosaic.minSensitivityForInclusionOfDivergentConeConnections, ...
        'maxNumberOfConesOutsideContour', 0, ...
        'overlappingWeightsDivergencePattern', 'inLineWithExclusiveConnections', ...  % choose from {'isotropic', 'inLineWithExclusiveConnections', 'orthogonalToExclusiveConnections'}
        'coneTypesIncluded', coneTypesPooledByRFcenter);

    % Generate RF center overlap
    theMRGCMosaic.generateRFcenterOverlap(rfCenterOverlapParams, ...
                'visualizeGenerationOfOverlappingRFcenterWeights', false, ...
                'minConeWeightVisualized', exp(-4), ...
                'visualizeMinSensitivityForInclusion', false, ...
                'visualizeSensitivityAtPointOfOverlap', false, ...
                'visualizeGauthierSensitivityAtPointOfOverlap', false);

    % Generate filename for exported mosaic
    [theCenterConnectedMRGCMosaicFullFileName, ~, ...
     theCenterConnectedMRGCMosaicFileName] = ...
        RGCMosaicConstructor.filepathFor.exportedMosaicFileName(...
            centerConnectedParamsStruct, 'center connected with overlap');
  
    % Export the RF center overlapping mosaic
    save(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic', '-v7.3');
    fprintf('The center-connected MRGC mosaic with RF center overlap (stage 2C) was saved to: \n\t%s\n', theCenterConnectedMRGCMosaicFullFileName);
end


if (visualizeMosaicAtStage2A)
    % Load a previously-generated mosaic at stage 2A
    centerConnectedParamsStruct = struct(...
        'whichEye', pStruct.whichEye, ...
        'eccentricityDegs', pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
        'sizeDegs', pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
        'spatialChromaticUniformityTradeoff', spatialChromaticUniformityTradeoff, ...
        'customLMSconeDensities', pStruct.customLMSconeDensities);
    
    [theCenterConnectedMRGCMosaicFullFileName, ~, ...
     theCenterConnectedMRGCMosaicFileName] = ...
        RGCMosaicConstructor.filepathFor.exportedMosaicFileName(...
                centerConnectedParamsStruct, 'center connected');

    fprintf('Loading center-connected (stage 2A) MRGC mosaic from\n\t%s\n', theCenterConnectedMRGCMosaicFullFileName);
    % Load the center-connected mRGCMosaic
    load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');
end

if (visualizeMosaicAtStage2C)
    % Load a previously-generated mosaic at stage 2C
    centerConnectedParamsStruct = struct(...
        'whichEye', pStruct.whichEye, ...
        'eccentricityDegs', pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
        'sizeDegs', pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
        'spatialChromaticUniformityTradeoff', spatialChromaticUniformityTradeoff, ...
        'customLMSconeDensities', pStruct.customLMSconeDensities);

    % Generate filename for exported mosaic
    [theCenterConnectedMRGCMosaicFullFileName, ~, ...
     theCenterConnectedMRGCMosaicFileName] = ...
        RGCMosaicConstructor.filepathFor.exportedMosaicFileName(...
            centerConnectedParamsStruct, 'center connected with overlap');
  
    fprintf('Loading center-connected (stage 2C) MRGC mosaic from\n\t%s\n', theCenterConnectedMRGCMosaicFullFileName);
    % Load the center-connected mRGCMosaic
    load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');
end


% Visualization params
% Whether to identify the input cones for each mRGC RF center
identifyInputConesInFullMosaic = ~true;

% Whether to identify the pooling of the input cones for each mRGC RF center
identifyPooledConesInFullMosaic = true;

% Visualization domain:
% If empty, we are visualizing the entire mosaic
domainVisualizationLimits = [];

% Or just visualize a portion by specifying the visualization domain limits
% minX
%domainVisualizationLimits(1) = 17;
% maxX
%domainVisualizationLimits(2) = 20;
% minY
%domainVisualizationLimits(3) =-0.3;   
% maxY
%domainVisualizationLimits(4) = 0.3;

% minConeWeightVisualized only has an effect for stage 2C mosaics
% (which contain divergent cone to RGC RF center connections, which had non-binary weights)

% Weight for a representation similar to that of Chichilnisky labs
minConeWeightVisualized = mRGCMosaic.sensitivityAtPointOfOverlap; 

% Weight = 10%, which is the noise floor in the measurements of cone weights 
% according to Greg Field and which is what they used in Fig 4 of their 2010 paper.
% minConeWeightVisualized = mRGCMosaic.minRFcenterConeWeightIncludedToMatchFigure4OfFieldEtAl2010

% Weight = 1%, the min sensitivity for inclusion of divergence cone connections
%minConeWeightVisualized = mRGCMosaic.minSensitivityForInclusionOfDivergentConeConnections

visualizationParamsStruct = struct(...
    'maxNumberOfConesOutsideContour', 0);

if (isempty(domainVisualizationLimits))
    figureFormat = PublicationReadyPlotLib.figureComponents('1x1 giant square mosaic');
else
    figureFormat = PublicationReadyPlotLib.figureComponents('1x1 double width figure');
end

RGCMosaicConstructor.visualize.fullMosaic(...
    theCenterConnectedMRGCMosaicFullFileName, ...
    theCenterConnectedMRGCMosaicFileName, ...
    minConeWeightVisualized, ...
    visualizationParamsStruct, ...
    'identifyInputCones', identifyInputConesInFullMosaic, ...
    'identifyPooledCones', identifyPooledConesInFullMosaic, ...
    'inputConesAlpha', 1.0, ...
    'withFigureFormat', figureFormat, ...
    'domainVisualizationLimits', domainVisualizationLimits, ...
    'domainVisualizationTicks', struct('x', -50:0.5:50, 'y', -50:0.5:50), ...
    'withoutPlotTitle', true);