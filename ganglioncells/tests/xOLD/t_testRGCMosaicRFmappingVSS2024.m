% Script to examine how optics, chromatic aberration and incomplete silence substitution affects
% m-sequence RF maps on mRGCMosaics at different eccentricities. 
% Script used in the analyses for the ReidShapley experiment (presented first on the MURI grant meeting)
%
% Usage:
%{
    t_testRGCMosaicRFmappingVSS2024
%}

% Initialize session
close all; clear all;

% Actions to perform
computeInputConeMosaicResponses = ~true;
computeMRGCMosaicResponses = ~true;
computeRFmaps = ~true;

visualizeSingleChromaticityRFmaps = ~true;
generateVisualRFandConePoolingMapComboPlots = ~true;

% Aggregate over multiple chromaticities
visualizeAggregatedChromaticityRFmaps = true;
onlyPlotPreviouslyComputedRoRincRatios = true;

visualizeInputConeMosaicResponses = ~true;
visualizePSFonTopOfConeMosaic = true;

if (visualizeAggregatedChromaticityRFmaps)
    computeInputConeMosaicResponses = false;
    computeMRGCMosaicResponses = false;
    computeRFmaps = false;
end

% targetVisualSTFdescriptorToOptimizeFor:
targetVisualSTFdescriptorToOptimizeFor = 'default';

% Macaque mosaic 
coneMosaicSpecies = 'macaque';

% Target eccentricity
rgcMosaicName = 'VSS2024TalkTemporal3DegsMosaic';  % For Lee&Shapley analyses
rgcMosaicName = 'VSS2024TalkTemporal7DegsMosaic';  % For Reid&Shapley analyses (their RFs were in the range 3-13 degs)

% Optics subject
opticsSubjectName = 'VSS2024TalkFirstSubject';
opticsSubjectName = 'VSS2024TalkSecondSubject';


% Target RGCs with surround specificity around 0.5 (i.e., 50/50 L/M cone net weight in the surround)
targetedSurroundPurityRange = [0 1];       % [0.4 0.6];
targetedRadialEccentricityRange = [];       % Empty means full range
targetedCenterConeNumerosityRange = [1 1];  % single cone centers only
targetedCenterConeNumerosityRange = [1 2];  % one and two cone centers only
%targetedCenterConeNumerosityRange = [];     % Empty means any # of center cones
targetedCenterPurityRange = [1 1];         % single type of cone in the RF center

% RF mapping parameters
% 1. Chromaticity for RFmapping
chromaticityForRFmapping = 'Achromatic';   % Choose between 'LconeIsolating', 'MconeIsolating', 'Achromatic'

% 2. Cone fundamentals design: whether to design cone-isolating stimuli taking into account the spectral tuning of cones
% at the mapped position
coneFundamentalsOptimizedForStimPosition = ~true;


% 3. RF resolution (pixels): simulate the low-res experiment of Reid&Shapley

rfMapPixelsAcross = 9;
mSequenceBitLength = 11;

%rfMapPixelsAcross = 11;
%mSequenceBitLength = 11;

%rfMapPixelsAcross = 13;
%mSequenceBitLength = 11;

% High-resolution relative to the one employed by  Reid&Shapley 
%rfMapPixelsAcross = 17;
%mSequenceBitLength = 11;

% Resolution appropriate to the one employed by  Reid&Shapley 
%rfMapPixelsAcross = 15;
%SequenceBitLength = 11;

%rfMapPixelsAcross = 45;
%mSequenceBitLength = 12;

%rfMapPixelsAcross = 21;
%mSequenceBitLength = 11;

%rfMapPixelsAcross = 63;
%mSequenceBitLength = 13;

%rfMapPixelsAcross = 89;
%mSequenceBitLength = 14;

%rfMapPixelsAcross = 127;
%mSequenceBitLength = 15;

% 4. optics modification to run
% Option 1. Native optics (what was used to optimize the surround),
opticsForRFmapping = 'nativeOptics';
residualWithRespectToNativeOpticsDefocusDiopters = [];

% Option 2. residual defocus with respect to native optics
%opticsForRFmapping = 'refractionResidualWithRespectToNativeOptics';
%residualWithRespectToNativeOpticsDefocusDiopters = 0.75;
%residualWithRespectToNativeOpticsDefocusDiopters = -0.50;
%residualWithRespectToNativeOpticsDefocusDiopters = -0.5;
%residualWithRespectToNativeOpticsDefocusDiopters = 0.25;

% Option 3. adaptive optics, 6 mm pupil, no LCA
%opticsForRFmapping = 'adaptiveOptics6MM';

% Option 4.  adaptive optics with LCA
%opticsForRFmapping = 'adaptiveOptics6MMwithLCA';

% Option 5. a custom refraction
%opticsForRFmapping = 'customRefraction';
%customRefractionDiopters = -0.25;


if (visualizeAggregatedChromaticityRFmaps)
    chromaticityForRFmapping = 'Achromatic';
end

if strcmp(opticsForRFmapping, 'adaptiveOptics6MM')
    profileGain = 5;
else
    profileGain = 20;
end
profileGain = 1;

% RF isoresponse contour z-levels
zLevelsNegative = -1/profileGain*(0.002:0.05:0.5);
zLevelsPositive =  1/profileGain*(0.002:0.05:0.5);



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


%Configure a conservative parpool manager. This gives at least 8 GB RAM/core
%ASPPManager = AppleSiliconParPoolManager(16);
%
ASPPManager = AppleSiliconParPoolManager('half max');

%ASPPManager = AppleSiliconParPoolManager('conservative');
%ASPPManager = AppleSiliconParPoolManager(10);

% Load theMRGCmosaic
load(theSurroundConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');
% Print info on how the surrounds where optimized
%RGCMosaicConstructor.helper.utils.printStruct(theMRGCMosaic.rfSurroundConnectivityParams, 'rfSurroundConnectivityParams')

% Position within the mRGCmosaic that we will compute m-sequence RFmaps.
% - the center of the mRGCmosaic
mappedPositionDegs = pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs;

% - or an arbitrary position withing the mRGC mosaic
mappedPositionDegs  = [-5.5 0]; %with the  mosaic

% Mapped size within the mRGCmosaic (we only map mRGCs whose center is within this region)
targetNumberOfMappedCells = 450;
mappedSizeDegs = RGCMosaicAnalyzer.compute.mSequenceRFmapSizeDegsForMosaic(theMRGCMosaic, mappedPositionDegs, targetNumberOfMappedCells);
fprintf('Pixel size will be %f arc min\n', mappedSizeDegs(1)*60/rfMapPixelsAcross);

% Append the mapped position info in the filename
mappedPositionString = sprintf('_@%2.1f_%2.1f.mat', mappedPositionDegs(1), mappedPositionDegs(2));
theInputConeMosaicRFmappingResponsesFullFileName = strrep(theInputConeMosaicRFmappingResponsesFullFileName, ...
    '.mat', mappedPositionString);
theMRGCMosaicRFmappingResponsesFullFileName = strrep(theMRGCMosaicRFmappingResponsesFullFileName, ...
    '.mat', mappedPositionString);

% Append the RF map resolution in the filename
theInputConeMosaicRFmappingResponsesFullFileName = strrep(theInputConeMosaicRFmappingResponsesFullFileName, ...
    '.mat', sprintf('_%dx.mat',rfMapPixelsAcross));
theMRGCMosaicRFmappingResponsesFullFileName = strrep(theMRGCMosaicRFmappingResponsesFullFileName, ...
    '.mat', sprintf('_%dx.mat',rfMapPixelsAcross));


% Extra support for input cone mosaic supporting the RF surrounds
extraSupportDegsForInputConeMosaic = 0.3*mRGCMosaic.extraSupportDegsForMidgetRGCSurrounds(abs(mappedPositionDegs)+mappedSizeDegs, 0.0);
visualizeSpatialRelationshipToSourceMosaic = true;

% Crop mosaic to the mapped region
fprintf('Trimming mRGCMosaic to %2.1f x %2.1f degs at (x,y) = (%2.1f,%2.1f)\n', ...
        mappedSizeDegs(1), mappedSizeDegs(2), mappedPositionDegs(1), mappedPositionDegs(2));
theMRGCMosaic.cropToSizeAtEccentricity(mappedSizeDegs, mappedPositionDegs, ...
        'extraSupportDegsForInputConeMosaic', extraSupportDegsForInputConeMosaic, ...
        'visualizeSpatialRelationshipToSourceMosaic', visualizeSpatialRelationshipToSourceMosaic);

theMRGCMosaic.visualize()

if (visualizeSingleChromaticityRFmaps) || (visualizeAggregatedChromaticityRFmaps)
    % Inidices of RGCs within mapped region
    theTargetRGCindicesWithinMappedRegion = theMRGCMosaic.indicesOfRGCsWithinROI(...
        mappedPositionDegs, [mappedSizeDegs(1) mappedSizeDegs(2)]);

    % Indices of RGCs within target property ranges
    theTargetRGCindicesWithinTargetRanges = theMRGCMosaic.indicesOfRGCsWithinTargetedPropertyRanges( ...
                    targetedCenterConeNumerosityRange, ...
                    targetedSurroundPurityRange, ...
                    targetedRadialEccentricityRange, ...
                    targetedCenterPurityRange);

    % Intersection of these 2 indices
    theTargetVisualizedRGCindices = intersect(theTargetRGCindicesWithinMappedRegion,theTargetRGCindicesWithinTargetRanges);
end


if (computeInputConeMosaicResponses)
    switch (opticsForRFmapping) 
        case  'refractionResidualWithRespectToNativeOptics'
            computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters, ...
                chromaticityForRFmapping, rfMapPixelsAcross, mSequenceBitLength, ...
                coneFundamentalsOptimizedForStimPosition, ...
                theInputConeMosaicRFmappingResponsesFullFileName, theMRGCMosaicRFmappingResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic);

         case 'customRefraction'
            computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, customRefractionDiopters, ...
                chromaticityForRFmapping, rfMapPixelsAcross, mSequenceBitLength, ...
                coneFundamentalsOptimizedForStimPosition, ...
                theInputConeMosaicRFmappingResponsesFullFileName, theMRGCMosaicRFmappingResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic);

        otherwise
             computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, [],  ...
                chromaticityForRFmapping, rfMapPixelsAcross, mSequenceBitLength, ...
                coneFundamentalsOptimizedForStimPosition, ...
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
    % No smoothing
    smoothRFmaps = false;
    visualizeSmoothingAndFitting = false;

    RGCMosaicAnalyzer.compute.ReidShapleyMequenceRFmapSingleRun(...
        theMRGCMosaicRFmappingResponsesFullFileName, ...
        surroundConnectedParamsStruct, ...
        smoothRFmaps, ...
        visualizeSmoothingAndFitting);
end

% Specify smoothing kernel params
[smoothinKernelSizePixels, smoothinKernelSigmaPixels] = smoothinKernelParams(rfMapPixelsAcross);

if (visualizeSingleChromaticityRFmaps)
    [smoothinKernelSizePixels, smoothinKernelSigmaPixels] = RGCMosaicAnalyzer.visualize.mSequenceRFmaps(...
        theMRGCMosaicRFmappingResponsesFullFileName, ...
        surroundConnectedParamsStruct, ...
        theMRGCMosaic, theTargetVisualizedRGCindices, ...
        profileGain, zLevelsNegative, zLevelsPositive, ...
        smoothinKernelSizePixels, smoothinKernelSigmaPixels, ...
        coneFundamentalAndSpatialResolutionString, ...
        chromaticityForRFmapping, opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters, ...
        generateVisualRFandConePoolingMapComboPlots, false);
    fprintf('RF map smoothing kernel employed: size = %2.4f, sigma = %2.4f\n', smoothinKernelSizePixels, smoothinKernelSigmaPixels);
end

if (visualizeAggregatedChromaticityRFmaps)

    aggregatedChromaticities = {'Achromatic', 'LconeIsolating', 'MconeIsolating'} ;
    for iChromaticity = 1:numel(aggregatedChromaticities)
        aggregatedMRGCMosaicRFmappingResponsesFullFileNames{iChromaticity} = ...
            strrep(theMRGCMosaicRFmappingResponsesFullFileName, 'Achromatic', aggregatedChromaticities{iChromaticity});
    end

    RGCMosaicAnalyzer.visualize.mSequenceRFmaps(...
        aggregatedMRGCMosaicRFmappingResponsesFullFileNames, ...
        surroundConnectedParamsStruct, ...
        theMRGCMosaic, theTargetVisualizedRGCindices, ...
        profileGain, zLevelsNegative, zLevelsPositive, ...
        smoothinKernelSizePixels, smoothinKernelSigmaPixels, ...
        coneFundamentalAndSpatialResolutionString, ...
        aggregatedChromaticities, opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters, ...
        false, onlyPlotPreviouslyComputedRoRincRatios);
end


%
% SUPPORTING FUNCTIONS
%

function computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, customRefractionDiopters, ...
                chromaticityForRFmapping, rfMapPixelsAcross, mSequenceBitLength, ...
                coneFundamentalsOptimizedForStimPosition, ...
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
        'coneFundamentalsOptimizedForStimPosition', coneFundamentalsOptimizedForStimPosition, ...
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


function [smoothinKernelSizePixels, smoothinKernelSigmaPixels] = smoothinKernelParams(rfMapPixelsAcross)
     % Or specify a value for both to employ a user-specificed RF smoothing kernel
    smoothinKernelSizePixels = rfMapPixelsAcross/3;
    smoothinKernelSigmaPixels = smoothinKernelSizePixels / 6;
    smoothinKernelSizePixels = round(smoothinKernelSizePixels);
    if (mod(smoothinKernelSizePixels,2) == 0)
        % Make it odd
        smoothinKernelSizePixels = smoothinKernelSizePixels + 1;
    end
end