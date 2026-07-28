% Script to run m-sequence RF maps on mRGCMosaics
% This script is also used to generate materials for the validation
% figures for the VSS2024 presentation
%
% Usage:
%{
    t_testRGCMosaicRFmappingVSS2024
%}

% Initialize session
close all; clear all;

% Configure a conservative parpool manager. This gives at least 8 GB RAM/core
%ASPPManager = AppleSiliconParPoolManager(12);
%
ASPPManager = AppleSiliconParPoolManager('half max');
%ASPPManager = AppleSiliconParPoolManager('conservative');


computeInputConeMosaicResponses = ~true;
computeMRGCMosaicResponses = ~true;
computeRFmaps = ~true;
visualizeSingleChromaticityRFmaps = ~true;
visualizeAggregatedChromaticityRFmaps = true;

if (visualizeAggregatedChromaticityRFmaps)
    computeInputConeMosaicResponses = false;
    computeMRGCMosaicResponses = false;
    computeRFmaps = false;
end

visualizeInputConeMosaicResponses = true;
visualizePSFonTopOfConeMosaic = ~true;

% If fixedRGCindices is not empty, we use these indices.
% If fixedRGCindices is empty, we find the indices of RGCs whose (achromaticBPIs,centerIsolatingBPIs) 
% are closest to the specified exemplarRGCs.achromaticBPIs, exemplarRGCs.centerIsolatingBPIs

% Indices of exemplar cells visualized
McenterExemplars = [4204 4257 4054]
LcenterExemplars = [5051 4271 5020 4591]
exemplarRGCs.fixedRGCindices = [5051 4271 5020 4591 4204 4257 4054] ; %[4271 4054];
% The above cells were selected for the folling BPI targets 
% for the max Strehl ratio responses (0.00 defocus)
exemplarRGCs.achromaticBPIs      = [0.48  0.48];
exemplarRGCs.centerIsolatingBPIs = [0.9  0.95 1.0];

% Target RGCs with surround specificity around 0.5 (i.e., 50/50 L/M cone net weight in the surround)
targetedSurroundPurityRange = [0 1];    % [0.4 0.6];
targetedRadialEccentricityRange = [];       % Empty means full range
targetedCenterConeNumerosityRange = [1 1];  % single cone centers only
  

% Choose chromaticity to run
chromaticityForRFmapping = 'MconeIsolating';   % Choose between 'LconeIsolating', 'MconeIsolating', 'Achromatic'
if (visualizeAggregatedChromaticityRFmaps)
    chromaticityForRFmapping = 'Achromatic';
end

% RF resolution (pixels)
rfMapPixelsAcross = 63;

% Choose degree of bandpassiness
targetVisualSTFdescriptor = 'default';               % mean of Croner & Kaplan with default Rc/Rs
targetVisualSTFdescriptor = 'higher intSCratio';     % Higher than mean of Croner & Kaplan with default Rc/Rs
%targetVisualSTFdescriptor = 'very high intSCratio';  % Much higher than mean Croner & Kaplan with higher Rc/Rs
%targetVisualSTFdescriptor = 'lower intSCratio';      % Lower than mean of Croner & Kaplan with default Rc/Rs
%targetVisualSTFdescriptor = 'ultra high intSCratio';  % Twice the mean Croner & Kaplan with default Rc/Rs

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

if strcmp(opticsForRFmapping, 'adaptiveOptics6MM')
    profileGain = 1;
else
    profileGain = 4;
end

% Option 4.  adaptive optics with LCA
%opticsForRFmapping = 'adaptiveOptics6MMwithLCA';

% Option 5. a custom refraction
% opticsForRFmapping = 'customRefraction';
% customRefractionDiopters = -4;

% Which mRGC/cone mosaic lattice to use
sourceLatticeSizeDegs = 64; 
whichEye = 'right eye';

% Human retina
customLMSconeDensities = [];

% Macaque retina
customLMSconeDensities = [0.45 0.45 0.1];

% Examined spatialChromaticUniformityTradeoff [0: minimize chromatic variance 1: minimize spatial variance]
spatialChromaticUniformityTradeoff = 1.0;

% Which optics to use
whichZernikeDataBase = 'Polans2015';

% Second compute-ready mosaic (What was used for the VSS2024 talk)
mosaicEccDegs = [-3 0]; mosaicSizeDegs = 1*[1 1];
whichSubjectID = 3;


% Extra size (margin)
extraSizeDegs = 1;
mosaicSizeDegs = mosaicSizeDegs+extraSizeDegs*[1 1];

% Use mosaic with RFcenter overlap
employRFCenterOverlappingMosaic = true;

surroundOptimizationStrategy = 'LowerLeftQH1paramsMediumVisualSTFparamTolerance';

% Generate filename
[theSurroundConnectedMRGCMosaicFullFileName, theInputConeMosaicRFmappingResponsesFullFileName, ...
    theMRGCMosaicRFmappingResponsesFullFileName, surroundConnectedParamsStruct, targetVisualSTFmodifierStruct] = ...
        RGCMosaicAnalyzer.filepathFor.testRuns(whichEye, whichZernikeDataBase, whichSubjectID, ...
            mosaicEccDegs, mosaicSizeDegs, employRFCenterOverlappingMosaic, ...
            spatialChromaticUniformityTradeoff, customLMSconeDensities, ...
            targetVisualSTFdescriptor, surroundOptimizationStrategy, ...
            opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters, ...
            chromaticityForRFmapping, ...
            'mSequenceResponses');


% Load theMRGCmosaic
load(theSurroundConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');

% Position of the mRGCmosaic that we will map.
% This is the retion of the included RF centers.
% The actuall mapped area is larger to include the surrounds
theTargetVisualizedRGCindex = exemplarRGCs.fixedRGCindices(2);
mappedPositionDegs = theMRGCMosaic.rgcRFpositionsDegs(theTargetVisualizedRGCindex,:)
mappedSizeDegs = [0.5 0.5];
targetRGCstatsBeforeCropping = theMRGCMosaic.singleCellConnectivityStats(theTargetVisualizedRGCindex, 'center')

% Add the mapped position info
mappedPositionString = sprintf('_@%2.1f_%2.1f.mat', mappedPositionDegs(1), mappedPositionDegs(2));
theInputConeMosaicRFmappingResponsesFullFileName = strrep(theInputConeMosaicRFmappingResponsesFullFileName, ...
    '.mat', mappedPositionString);
theMRGCMosaicRFmappingResponsesFullFileName = strrep(theMRGCMosaicRFmappingResponsesFullFileName, ...
    '.mat', mappedPositionString);



% Crop mosaic
fprintf('Trimming mRGCMosaic to %2.1 x %2.1f degs at (x,y) = (%2.1f,%2.1f)\n', ...
        mappedPositionDegs(1), mappedPositionDegs(2), mappedSizeDegs(1), mappedSizeDegs(2));

% Extra support for input cone mosaic supporting the RF surrounds
extraSupportDegsForInputConeMosaic = 0.3*mRGCMosaic.extraSupportDegsForMidgetRGCSurrounds(abs(mappedPositionDegs)+mappedSizeDegs, 0.0);
visualizeSpatialRelationshipToSourceMosaic = true;

theMRGCMosaic.cropToSizeAtEccentricity(mappedSizeDegs, mappedPositionDegs, ...
        'extraSupportDegsForInputConeMosaic', extraSupportDegsForInputConeMosaic, ...
        'visualizeSpatialRelationshipToSourceMosaic', visualizeSpatialRelationshipToSourceMosaic);

idx = theMRGCMosaic.indicesOfRGCsWithinROI(mappedPositionDegs, [0.1 0.1])
theTargetVisualizedRGCindex = idx(1);
targetRGCstatsAfterCropping = theMRGCMosaic.singleCellConnectivityStats(theTargetVisualizedRGCindex, 'center');


if (computeInputConeMosaicResponses)
    switch (opticsForRFmapping) 
        case  'refractionResidualWithRespectToNativeOptics'
            computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters, ...
                chromaticityForRFmapping, rfMapPixelsAcross, ...
                theInputConeMosaicRFmappingResponsesFullFileName, theMRGCMosaicRFmappingResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic);

         case 'customRefraction'
            computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, customRefractionDiopters, ...
                chromaticityForRFmapping, rfMapPixelsAcross, ...
                theInputConeMosaicRFmappingResponsesFullFileName, theMRGCMosaicRFmappingResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic);

        otherwise
             computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, [],  ...
                chromaticityForRFmapping, rfMapPixelsAcross, ...
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
    RGCMosaicAnalyzer.compute.ReidShapleyMequenceRFmapSingleRun(...
        theMRGCMosaicRFmappingResponsesFullFileName, ...
        surroundConnectedParamsStruct ...
        );
end

if (visualizeSingleChromaticityRFmaps)
    RGCMosaicAnalyzer.visualize.mSequenceRFmaps(...
        theMRGCMosaicRFmappingResponsesFullFileName, ...
        surroundConnectedParamsStruct, ...
        theMRGCMosaic, ...
        theTargetVisualizedRGCindex, ...
        profileGain, ...
        chromaticityForRFmapping, ...
        opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters ...
        );
end

if (visualizeAggregatedChromaticityRFmaps)
    aggregatedChromaticities = {'Achromatic', 'LconeIsolating', 'MconeIsolating'} ;
    for iChromaticity = 1:numel(aggregatedChromaticities)
        aggregatedMRGCMosaicRFmappingResponsesFullFileNames{iChromaticity} = ...
            strrep(theMRGCMosaicRFmappingResponsesFullFileName, 'Achromatic', aggregatedChromaticities{iChromaticity});
    end

    theTargetVisualizedRGCindices = [6 42 28];

    RGCMosaicAnalyzer.visualize.mSequenceRFmaps(...
        aggregatedMRGCMosaicRFmappingResponsesFullFileNames, ...
        surroundConnectedParamsStruct, ...
        theMRGCMosaic, ...
        theTargetVisualizedRGCindices, ...
        profileGain, ...
        aggregatedChromaticities, ...
        opticsForRFmapping, residualWithRespectToNativeOpticsDefocusDiopters ...
        );
end


%
% SUPPORTING FUNCTIONS
%

function computeMSequenceResponses(theMRGCMosaic, opticsForRFmapping, customRefractionDiopters, ...
                chromaticityForRFmapping, rfMapPixelsAcross, ...
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
        'mSequenceBitLength', 13, ...
        'frameBatchSize', 256, ...
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