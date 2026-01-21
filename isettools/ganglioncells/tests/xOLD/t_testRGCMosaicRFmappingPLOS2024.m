% Script to run m-sequence RF maps on mRGCMosaics
% This script is also used to generate materials for the validation
% figures for the PLOS2024 paper
%
% Usage:
%{
    t_testRGCMosaicRFmappingPLOS2024
%}

% Initialize session
close all; clear all;
    
% Actions to perform
computeInputConeMosaicResponses = true;
computeMRGCMosaicResponses = ~true;
computeRFmaps = ~true;

visualizeSingleChromaticityRFmaps  = ~true;

% Whether to visualize the RF map and the cone pooling map for each cell
generateVisualRFandConePoolingMapComboPlots = ~true;

% Whether to visualize the RF map of each cell together with each nearest neighbor RF map
% (Good for visualized RF overlap)
generateNearestNeighborOverlapPlots = ~true;

% Aggregation over multiple eccentricities
aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities = ~true;

% Aggregate over multiple chromaticities
visualizeAggregatedChromaticityRFmaps = false;

visualizeInputConeMosaicResponses = ~true;
visualizePSFonTopOfConeMosaic = ~true;

if (visualizeAggregatedChromaticityRFmaps)
    computeInputConeMosaicResponses = false;
    computeMRGCMosaicResponses = false;
    computeRFmaps = false;
end

if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)
    [surroundOptimizationStrategyToBeAggregated, ...
     targetNumberOfMappedCellsToBeAggregated, ...
     targetHorizontalEccentricitiesToBeAggregated, ...
     targetMosaicSizesToBeAggregated, ...
     targetMappedPositionDegsToBeAggregated ] = testsHelperFunctions.assembleEccentricityMatrices();
end

% Target RGCs with surround specificity around 0.5 (i.e., 50/50 L/M cone net weight in the surround)
targetedSurroundPurityRange = [0 1];     % Full range [0.4 0.6];
targetedRadialEccentricityRange = [];    % Empty means full range
targetedCenterConeNumerosityRange = [];  % Empty means full range
  

% Choose chromaticity to run
chromaticityForRFmapping = 'Achromatic';   % Choose between 'LconeIsolating', 'MconeIsolating', 'Achromatic'

if (visualizeAggregatedChromaticityRFmaps)
    chromaticityForRFmapping = 'Achromatic';
end

coneFundamentalsOptimizedForStimPosition = false;

% RF resolution (pixels)

% Did this for mosaics at eccentricities < 11 degrees

rfMapPixelsAcross = 45;
mSequenceBitLength = 12;

% for 25 degs (run imdex 52_
rfMapPixelsAcross = 63;
mSequenceBitLength = 13;

% Did this for mosaics at eccentricities > 10 degrees
%rfMapPixelsAcross = 89;
%mSequenceBitLength = 14;

%rfMapPixelsAcross = 127;
%mSequenceBitLength = 15;


% Choose optics to run
% Option 1. Native optics (what was used to optimize the surround)
opticsForRFmapping = 'nativeOptics';
residualWithRespectToNativeOpticsDefocusDiopters = [];

    % Option 2. residual defocus with respect to native optics
    %opticsForRFmapping = 'refractionResidualWithRespectToNativeOptics';
    %residualWithRespectToNativeOpticsDefocusDiopters = 0.75;
    %residualWithRespectToNativeOpticsDefocusDiopters = 0.50;
    %residualWithRespectToNativeOpticsDefocusDiopters = 0.25;
    %residualWithRespectToNativeOpticsDefocusDiopters =-0.25;

    % Option 3. adaptive optics, 6 mm pupil, no LCA
   opticsForRFmapping = 'adaptiveOptics6MM';

    % Option 4.  adaptive optics with LCA
    %opticsForRFmapping = 'adaptiveOptics6MMwithLCA';

    % Option 5. a custom refraction
    % opticsForRFmapping = 'customRefraction';
    % customRefractionDiopters = -4;

if strcmp(opticsForRFmapping, 'adaptiveOptics6MM')
    profileGain = 5;
else
    profileGain = 20;
end

profileGain = 1;

% RF isoresponse contour z-levels
zLevelsNegative = -1/profileGain*(0.002:0.05:0.5);
zLevelsPositive =  1/profileGain*(0.002:0.05:0.5);

% targetVisualSTFdescriptorToOptimizeFor: choose from
% 'default';              % Mean of Croner & Kaplan with default Rc/Rs
% 'x1.5 RsRcRatio';       % 50% higher than Croner & Kaplan
% 'higher RsRcRatio';
% 'higher intSCratio';
% 'very high intSCratio';
% 'lower intSCratio';
% 'ultra high intSCratio';
targetVisualSTFdescriptorToOptimizeFor = 'default';

% Human mosaic 
coneMosaicSpecies = 'human';

% Run with
% mRGCMosaic.amplificationInCenterOnlySensitivityCausedByInactiveSurrounds 0.6 
% (best agreement in RF center overlap with Gauthier)
rgcMosaicName = 'PLOSpaperTemporal32DegsMosaic';
rgcMosaicName = 'PLOSpaperTemporal25DegsMosaic';
%rgcMosaicName = 'PLOSpaperTemporal19DegsMosaic';
%rgcMosaicName = 'PLOSpaperTemporal7DegsMosaic';
%rgcMosaicName = 'PLOSpaperTemporal4DegsMosaic';
%rgcMosaicName = 'PLOSpaperTemporal2DegsMosaic';
%rgcMosaicName = 'PLOSpaperFovealMosaic';

% Optics subject
opticsSubjectName = 'PLOSpaperDefaultSubject';

% Generate params struct
pStruct = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
    coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptorToOptimizeFor);


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


if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)

    % Initialize vectors to be aggregated
    radialEccDegs = [];
    radialEccMMs = [];
    majorSigmaDegs = [];
    minorSigmaDegs = [];
    sigmaDegs = [];
    majorSigmaMicrons = [];
    minorSigmaMicrons = [];
    sigmaMicrons = [];

    % Agrregate data from all runs
    runIndicesToBeAggregated = 1:numel(targetHorizontalEccentricitiesToBeAggregated);
    runIndicesToBeAggregated = 1:57; % [34:37 38:56]
    
    %mminEccMMs = 6.0;
    %maxEccMMs = 6.6;
    minEccMMs = [];
    maxEccMMs = [];

    
    for iRun = runIndicesToBeAggregated
        mosaicEccDegs(1) = targetHorizontalEccentricitiesToBeAggregated(iRun);
        mosaicEccDegs(2) = pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs(2);
        aggregatedEccString = sprintf('Ecc%2.1f_%2.1f_Size%2.1fx%2.1f',mosaicEccDegs(1), mosaicEccDegs(2), ...
                targetMosaicSizesToBeAggregated(iRun,1),  targetMosaicSizesToBeAggregated(iRun,1));
        
        theSurroundConnectedMRGCMosaicFullFileName = strrep(...
            strrep(theOriginalSurroundConnectedMRGCMosaicFullFileName, targetString, aggregatedEccString), ...
            targetSurroundOptimizationStrategy, surroundOptimizationStrategyToBeAggregated{iRun});
        
    

        theMRGCMosaicRFmappingResponsesFullFileName = strrep(...
            strrep(theOriginalMRGCMosaicRFmappingResponsesFullFileName, targetString, aggregatedEccString), ...
            targetSurroundOptimizationStrategy, surroundOptimizationStrategyToBeAggregated{iRun});
    
        
        % Load theMRGCmosaic
        load(theSurroundConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');

        % Add the mapped position info
        mappedPositionDegs = targetMappedPositionDegsToBeAggregated(iRun,:);
        mappedSizeDegs = RGCMosaicAnalyzer.compute.mSequenceRFmapSizeDegsForMosaic(theMRGCMosaic, mappedPositionDegs, targetNumberOfMappedCellsToBeAggregated{iRun});
        
        mappedPositionString = sprintf('_@%2.1f_%2.1f.mat', mappedPositionDegs(1), mappedPositionDegs(2));
        theMRGCMosaicRFmappingResponsesFullFileName = strrep(theMRGCMosaicRFmappingResponsesFullFileName, ...
            '.mat', mappedPositionString);
        
        % Crop the mosaic
        fprintf('Trimming mRGCMosaic to %2.2f x %2.2f degs at (x,y) = (%2.1f,%2.1f)\n', ...
            mappedSizeDegs(1), mappedSizeDegs(2), mappedPositionDegs(1), mappedPositionDegs(2));

        % Extra support for input cone mosaic supporting the RF surrounds
        extraSupportDegsForInputConeMosaic = 0.3*mRGCMosaic.extraSupportDegsForMidgetRGCSurrounds(abs(mappedPositionDegs)+mappedSizeDegs, 0.0);
        visualizeSpatialRelationshipToSourceMosaic = ~true;

        theMRGCMosaic.cropToSizeAtEccentricity(mappedSizeDegs, mappedPositionDegs, ...
                'extraSupportDegsForInputConeMosaic', extraSupportDegsForInputConeMosaic, ...
                'visualizeSpatialRelationshipToSourceMosaic', visualizeSpatialRelationshipToSourceMosaic);


        if (~isempty(minEccMMs))&&(~isempty(maxEccMMs))
            eccMMs = theMRGCMosaic.angularEccDegsToLinearEccMMs(theMRGCMosaic.eccentricityDegs);
            mosaicRadialEccMMs = sqrt(sum(eccMMs.^2,2));
            if (mosaicRadialEccMMs<minEccMMs)||(mosaicRadialEccMMs>maxEccMMs)
                continue;
            end
            fprintf('Will accumulate run %d with mosaic from %s\n', iRun, theSurroundConnectedMRGCMosaicFullFileName);
        end

        % Generate theRFmapsFileName
        theRFmapsFileName = strrep(theMRGCMosaicRFmappingResponsesFullFileName, 'mSequenceResponses', 'demos/ReidShapleyAnalyses');
        theRFmapsFileName = strrep(theRFmapsFileName, '.mat', 'RF.mat');
        theRFmapsFileName = strrep(theRFmapsFileName, '_mRGCMosaic', '');

        % Load the RFmaps data file
        load(theRFmapsFileName, ...
            'theTemporalEquivalentEccentricitiesDegs', ...
            'theTemporalEquivalentEccentricitiesMMs', ...
            'theRFmaps', ...
            'spatialSupportDegs', ...
            'theFittedGaussianEllipsoids');

        theTargetVisualizedRGCindices = theMRGCMosaic.indicesOfRGCsWithinROI(mappedPositionDegs, [mappedSizeDegs(1) mappedSizeDegs(2)]);
        assert(numel(theTargetVisualizedRGCindices) == size(theTemporalEquivalentEccentricitiesDegs,1), ...
            sprintf('target RGCindices (%d) does not match analyzed data (%d)', ...
            numel(theTargetVisualizedRGCindices), size(theTemporalEquivalentEccentricitiesDegs,1) ...
            ));

        % Accumulate data
        numel(theTargetVisualizedRGCindices)
        for idx = 1:numel(theTargetVisualizedRGCindices)
            % The RGC position in degs
            theTargetVisualizedRGCindex = theTargetVisualizedRGCindices(idx);
            theRGCpositionDegs = theMRGCMosaic.rgcRFpositionsDegs(theTargetVisualizedRGCindex,:);

            % The fitted RF map
            theRFmap = theRFmaps{theTargetVisualizedRGCindex};
            theFittedEllipsoid  = theFittedGaussianEllipsoids{theTargetVisualizedRGCindex};

            % Radial eccentricities
            radialEccDegs(numel(radialEccDegs)+1) = sqrt(sum(theRGCpositionDegs.^2,2));
            radialEccMMs(numel(radialEccMMs)+1) = sqrt(sum(theRGCpositionDegs.^2,2)) * theMRGCMosaic.inputConeMosaic.micronsPerDegree * 1e-3;

            % Sigmas in degs
            majorSigmaDegs(numel(majorSigmaDegs)+1) = theFittedEllipsoid.majorSigmaDegs;
            minorSigmaDegs(numel(minorSigmaDegs)+1) = theFittedEllipsoid.minorSigmaDegs;
            sigmaDegs(numel(sigmaDegs)+1) = sqrt(theFittedEllipsoid.majorSigmaDegs*theFittedEllipsoid.minorSigmaDegs);

            % Sigmas in microns
            majorSigmaMicrons(numel(majorSigmaMicrons)+1) = theFittedEllipsoid.majorSigmaMicrons;
            minorSigmaMicrons(numel(minorSigmaMicrons)+1) = theFittedEllipsoid.minorSigmaMicrons;
            sigmaMicrons(numel(sigmaMicrons)+1) = sqrt(theFittedEllipsoid.majorSigmaMicrons*theFittedEllipsoid.minorSigmaMicrons);
        end
    end % for iRun

    % Finally, plot the fitted Gaussian radii as a function of eccentricity together with EJ's data
    % Filename for PDF graphic
    if (strcmp(opticsForRFmapping, 'refractionResidualWithRespectToNativeOptics'))
        opticsInfo = sprintf('%s_residualRefraction%2.2fD',opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters);
    else
        opticsInfo = opticsForRFmapping;
    end
    pdfFileName = sprintf('%s_%s_RFdiametersAsAFunctionOfEccentricity.pdf', opticsInfo, chromaticityForRFmapping);

    % Figure format
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    ff.box = 'off';
    ff.tickDir = 'in';
    ff.grid = 'off';
    ff.backgroundColor = 'none';

    % Render the plot
    displayMinorAndMajorAxisDiameters = false;
	employLogXaxis = true;
    employLogYaxis = ~true;


    RGCMosaicAnalyzer.visualize.mosaicRetinalDiametersAgainstMacaqueInVitroRFdiameters(201, ff, ...
        radialEccMMs, minorSigmaMicrons, majorSigmaMicrons, sigmaMicrons, ...
        displayMinorAndMajorAxisDiameters, employLogYaxis, employLogXaxis, ...
        pdfFileName);

    return;
end  % if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities)



%Configure a conservative parpool manager. This gives at least 8 GB RAM/core
%ASPPManager = AppleSiliconParPoolManager(16);
%
%ASPPManager = AppleSiliconParPoolManager('half max');
ASPPManager = AppleSiliconParPoolManager(8);

%ASPPManager = AppleSiliconParPoolManager('conservative');
%ASPPManager = AppleSiliconParPoolManager(10);


% This determines which mosaic data is analyzed
[surroundOptimizationStrategyToBeAggregated, ...
     targetNumberOfMappedCellsToBeAggregated, ...
     targetHorizontalEccentricitiesToBeAggregated, ...
     targetMosaicSizesToBeAggregated, ...
     targetMappedPositionDegsToBeAggregated ] = testsHelperFunctions.assembleEccentricityMatrices();


theRunIndexToCompute = 51 % 10: 19  33; 51 ..57
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

% Print info on how the surrounds where optimized
%RGCMosaicConstructor.helper.utils.printStruct(theMRGCMosaic.rfSurroundConnectivityParams, 'rfSurroundConnectivityParams')

% Crop mosaic to the mapped region
fprintf('Trimming mRGCMosaic to %2.1f x %2.1f degs at (x,y) = (%2.1f,%2.1f)\n', ...
        mappedSizeDegs(1), mappedSizeDegs(2), mappedPositionDegs(1), mappedPositionDegs(2));

% Extra support for input cone mosaic supporting the RF surrounds
extraSupportDegsForInputConeMosaic = 0.3*mRGCMosaic.extraSupportDegsForMidgetRGCSurrounds(abs(mappedPositionDegs)+mappedSizeDegs, 0.0);
visualizeSpatialRelationshipToSourceMosaic = true;

theMRGCMosaic.cropToSizeAtEccentricity(mappedSizeDegs, mappedPositionDegs, ...
        'extraSupportDegsForInputConeMosaic', extraSupportDegsForInputConeMosaic, ...
        'visualizeSpatialRelationshipToSourceMosaic', visualizeSpatialRelationshipToSourceMosaic);

if (visualizeSingleChromaticityRFmaps) || (visualizeAggregatedChromaticityRFmaps)
    theTargetVisualizedRGCindices = theMRGCMosaic.indicesOfRGCsWithinROI(mappedPositionDegs, [mappedSizeDegs(1) mappedSizeDegs(2)]);
end


if (computeInputConeMosaicResponses)
    switch (opticsForRFmapping) 
        case  'refractionResidualWithRespectToNativeOptics'
            computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters, ...
                chromaticityForRFmapping, rfMapPixelsAcross, mSequenceBitLength, ...
                theInputConeMosaicRFmappingResponsesFullFileName, theMRGCMosaicRFmappingResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic);

         case 'customRefraction'
            computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, customRefractionDiopters, ...
                chromaticityForRFmapping, rfMapPixelsAcross, mSequenceBitLength, ...
                theInputConeMosaicRFmappingResponsesFullFileName, theMRGCMosaicRFmappingResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic);

        otherwise
             computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, [],  ...
                chromaticityForRFmapping, rfMapPixelsAcross, mSequenceBitLength, ...
                theInputConeMosaicRFmappingResponsesFullFileName, theMRGCMosaicRFmappingResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic);
    end % switch
end % if (computeInputConeMosaicResponses)

if (computeMRGCMosaicResponses)
    RGCMosaicAnalyzer.simulateExperiment.mSequenceRFmapping(...
            theMRGCMosaic, [], [], ...
            theInputConeMosaicRFmappingResponsesFullFileName, ...
            theMRGCMosaicRFmappingResponsesFullFileName, ...
            'computeInputConeMosaicResponses', false, ...
            'computeMRGCMosaicResponses', true, ...
            'visualizeInputConeMosaicResponses', false);
end %if (computeMRGCMosaicResponses)

if (computeRFmaps)
    % How much smoothing to apply to the RF map before fitting the Gaussian
    smoothFittedRFmaps = ~true;

    % If visualizeSmoothingAndFitting is true, we do a serial fitting
    % If it is set to false, we use parfor (much faster)
    visualizeSmoothingAndFitting = ~true;

    RGCMosaicAnalyzer.compute.ReidShapleyMequenceRFmapSingleRun(...
        theMRGCMosaicRFmappingResponsesFullFileName, ...
        surroundConnectedParamsStruct, ...
        smoothFittedRFmaps, ...
        visualizeSmoothingAndFitting);
end


if (visualizeSingleChromaticityRFmaps)

    % Specify smoothing kernel params 
    [smoothinKernelSizePixels, smoothinKernelSigmaPixels] = testsHelperFunctions.smoothinKernelParams();

    coneFundamentalAndSpatialResolutionString = '';
    employTransparentBackground = false;

    [smoothinKernelSizePixels, smoothinKernelSigmaPixels] = RGCMosaicAnalyzer.visualize.mSequenceRFmaps(...
        theMRGCMosaicRFmappingResponsesFullFileName, ...
        surroundConnectedParamsStruct, ...
        theMRGCMosaic, theTargetVisualizedRGCindices, ...
        profileGain, zLevelsNegative, zLevelsPositive, ...
        smoothinKernelSizePixels, smoothinKernelSigmaPixels, ...
        coneFundamentalAndSpatialResolutionString, ...
        chromaticityForRFmapping, opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters, ...
        generateVisualRFandConePoolingMapComboPlots, ...
        generateNearestNeighborOverlapPlots, ...
        false, ...
        employTransparentBackground);
    fprintf('RF map smoothing kernel employed: size = %2.4f, sigma = %2.4f\n', smoothinKernelSizePixels, smoothinKernelSigmaPixels);
end

if (visualizeAggregatedChromaticityRFmaps)
    aggregatedChromaticities = {'Achromatic', 'LconeIsolating', 'MconeIsolating'} ;
    for iChromaticity = 1:numel(aggregatedChromaticities)
        aggregatedMRGCMosaicRFmappingResponsesFullFileNames{iChromaticity} = ...
            strrep(theMRGCMosaicRFmappingResponsesFullFileName, 'Achromatic', aggregatedChromaticities{iChromaticity});
    end

    employTransparentBackground = false;
    
    RGCMosaicAnalyzer.visualize.mSequenceRFmaps(...
        aggregatedMRGCMosaicRFmappingResponsesFullFileNames, ...
        surroundConnectedParamsStruct, ...
        theMRGCMosaic, theTargetVisualizedRGCindices, ...
        profileGain, zLevelsNegative, zLevelsPositive, ...
        aggregatedChromaticities, opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters, ...
        false, ...
        employTransparentBackground);
end


%
% SUPPORTING FUNCTIONS
%

function computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, customRefractionDiopters, ...
                chromaticityForRFmapping, rfMapPixelsAcross, mSequenceBitLength, ...
                theInputConeMosaicRFmappingResponsesFullFileName, ...
                theMRGCMosaicRFmappingResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic)

    % Generate the optics
    [theOI, thePSF] = RGCMosaicAnalyzer.compute.opticsForResponses(...
        theMRGCMosaic, opticsForRFmapping, customRefractionDiopters, visualizePSFonTopOfConeMosaic);

    % Determine the stimulus pixel resolution to be a fraction of the minimum cone aperture or cone spacing in the mosaic
    % here, half of the cone spacing
    theMetric = 'cone aperture';  % choose from {'cone aperture' or cone spacing'}
    mosaicRadialEcc = sqrt(sum(theMRGCMosaic.eccentricityDegs.^2,2));
    if ((strcmp(opticsForRFmapping, 'adaptiveOptics6MM')) || (strcmp(opticsForRFmapping, 'adaptiveOptics6MMwithLCA')))
        theFraction = 0.1;
    else
        theFraction = 0.25;
    end
    targetRGCindices =  1:theMRGCMosaic.rgcsNum;
    stimulusResolutionDegs = RGCMosaicConstructor.helper.simulateExperiment.stimulusResolutionFromConeApertureOrConeSpacing(...
                theMRGCMosaic, targetRGCindices, theFraction, theMetric);


    RFmappingParamsStruct = struct(...
        'backgroundChromaticity', [0.301 0.301], ...
        'backgroundLuminanceCdM2', 50.0, ...
        'chromaticity', chromaticityForRFmapping, ...
        'coneFundamentalsOptimizedForStimPosition', false, ...
        'ternaryInsteadOfBinaryMsequence', false, ...
        'rfPixelsAcross', rfMapPixelsAcross, ...
        'mSequenceBitLength', mSequenceBitLength,  ...
        'frameBatchSize', 128, ...
        'optimalResolutionDegs', stimulusResolutionDegs, ...               
        'positionDegs', theMRGCMosaic.inputConeMosaic.eccentricityDegs, ... 
        'sizeDegs', theMRGCMosaic.inputConeMosaic.sizeDegs);

    RGCMosaicAnalyzer.simulateExperiment.mSequenceRFmapping(...
            theMRGCMosaic, theOI, RFmappingParamsStruct, ...
            theInputConeMosaicRFmappingResponsesFullFileName, ...
            theMRGCMosaicRFmappingResponsesFullFileName, ...
            'computeInputConeMosaicResponses', true, ...
            'computeMRGCMosaicResponses', false, ...
            'visualizeInputConeMosaicResponses', visualizeInputConeMosaicResponses);

end