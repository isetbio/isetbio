% Script to compute retinal-space referred STFs on mRGCMosaics
% This script is also used to generate materials for the validation
% figures for the PLOS2024 paper
%
% Usage:
%{
    t_testRGCMosaicRetinalSTFsPLOS2024
%}

% Initialize session
close all; 
clear all;
    
% Whether to compute retinal space -- referred STFs (adaptive optics, in-vitro)
% of visual space referred STFS (physiological optics)
retinalSpaceReferredSTFs = true;

% Actions to perform: compute STF responses
computeInputConeMosaicResponses = ~true;
computeMRGCMosaicResponses = ~true;

% Set this to true to visualize responses of single cones
% Only happens if computeInputConeMosaicResponses is true
debugInputConeMosaicPcurrentResponse = ~true;


% Actions to perform: analyze computed STF responses
performCronerKaplanAnalysis = true;

% L and M cone isolating STFs if comparing to Wool
compareToWoolEtAl = ~true;

% Whether to re-analyze the STF data or load previously analyzed data
reAnalyzeSTFData = true;

% Visualization options
onlyVisualizeFittedDOGmodels = ~true;

% Whether to visualize the input cone mosaic responses
visualizeInputConeMosaicResponses = ~true;
visualizePSFonTopOfConeMosaic = ~true;

% Whether to visualize the full (2D) STF and the selected STF slice to be fitted
visualizeFullAndMaximalExcursionSTF = ~true;

% Whether to visualize the fitted STF
% Only has an effect when the model is being fit
visualizeSTFfits = ~true;

% Whether to visualize the constrained model parameter values and their range
% AND the sinusoidal fits to the time courses
% Only has an effect when the model is being fit
visualizeModelFitting = ~true;

% Whether to visuazize the fitted STF together
% with the retinal cone pooling map
visualizeSTFwithConeWeightsMap = true;


% Whether to visualize single cell cone-isolating STFs
visualizeSingleRGCConeIsolatingSTFsForOneMosaic = ~true;
singleRGCindexToAnalyze = [];


% Alter theMRGCMosaicSTFResponsesFullFileName is we simulate nonlinearities
mRGCNonLinearityParams = struct('type', 'photocurrent');
mRGCNonLinearityParams = [];

% Custom TF and contrast
customTemporalFrequencyAndContrast = [];  % the defaults, i.e., 1 Hz, 0.75 achromatic contrast

% non-default stimulus connfigurations. This will change the filenames

% Most-linear pCurrent configurationdrop
% 4Hz, 25% contrast, backgroundLuminance: 1x
customTemporalFrequencyAndContrast = struct(...
    'temporalFrequencyHz', 4.0, ...
    'totalContrast', 0.25, ...
    'backgroundLuminanceMultiplier', 1);


customTemporalFrequencyAndContrast = struct(...
    'temporalFrequencyHz', 1.0, ...
    'totalContrast', 0.25, ...
    'backgroundLuminanceMultiplier', 1);


customTemporalFrequencyAndContrast = struct(...
    'temporalFrequencyHz', 2.0, ...
    'totalContrast', 0.50, ...
    'backgroundLuminanceMultiplier', 1);

customTemporalFrequencyAndContrast = struct(...
    'temporalFrequencyHz', 4, ...
    'totalContrast', 0.25, ...
    'backgroundLuminanceMultiplier', 1);


customTemporalFrequencyAndContrast = [];

if (visualizeSingleRGCConeIsolatingSTFsForOneMosaic)

    singleRGCConeIsolatingSTFeccMMs = 3.0;

    % All cells in the selected mosaic
    singleRGCindexToAnalyze = [];

    % Specific cell in the selected mosaic
    % The below example cells are for singleRGCConeIsolatingSTFeccMMs = 3.0;
    %% High RF center cone purity (L) - L/M opponent response at 0 c/deg
    singleRGCindexToAnalyze = 108;

    % High RF center cone purity (M) - L/M opponent response at 0 c/deg
    singleRGCindexToAnalyze = 107;

    % Medium RF center cone purity (L) - L/M opponent response at 0 c/deg
    singleRGCindexToAnalyze = 84;

    % Medium RF center cone purity (M) - L/M opponent response at 0 c/deg
    singleRGCindexToAnalyze = 33;

    % Low RF center cone purity (L) - L/M non-opponent response at 0 c/deg
    singleRGCindexToAnalyze = 152;

    % Low  RF center cone purity (M) - L/M non-opponent response at 0 c/deg
    singleRGCindexToAnalyze = 67;

    singleRGCConeIsolatingSTFeccMMs = 3.9

    singleRGCindexToAnalyze = [];

    if (1==2)
        singleRGCConeIsolatingSTFeccMMs = 3.7;

        % Low  RF center cone purity (L) - L/M non-opponent response at 0 c/deg
        singleRGCindexToAnalyze = 148;

        % Medium  RF center cone purity (L) - notch

        % High RF center cone purity (L) - L/M opponent response at 0 c/deg
        singleRGCindexToAnalyze = 135;

        % All cells in the selected mosaic
        singleRGCindexToAnalyze = [];
    end

end


% Human mosaic 
coneMosaicSpecies = 'human';

% Optics subject for which the mRGC mosaics were optimized for
opticsSubjectName = 'PLOSpaperDefaultSubject';  % This is subject #2 (ranked #3), default
opticsSubjectName = 'VSS2024TalkFirstSubject';  % This is subject #3 (ranked #7), second subject

% Target eccentricity. Dummy, gets rewritten
RGCmosaicName = 'PLOSpaperTemporal32DegsMosaic';


% Optics for STF responses
if (retinalSpaceReferredSTFs)  
    opticsForSTFresponses = 'adaptiveOptics6MM';
    residualWithRespectToNativeOpticsDefocusDiopters = [];
else
    opticsForSTFresponses = 'nativeOptics';
    residualWithRespectToNativeOpticsDefocusDiopters = [];
end

% Which runs to compute (eccentricity range)
minEccToAnalyze = 0.01;
maxEccToAnalyze = 45.0;



% Did with pCurrent

%minEccToAnalyze = 16;
%maxEccToAnalyze = 16.0;

minEccToAnalyze = 8;
maxEccToAnalyze = 8.0;

%minEccToAnalyze = 4.5;
%maxEccToAnalyze = 4.6;

%minEccToAnalyze = 1.0;
%maxEccToAnalyze = 1.1;

%minEccToAnalyze = 0.1;
%maxEccToAnalyze = 0.12;


minEccToAnalyze = 3.0;
maxEccToAnalyze = 3.1


%minEccToAnalyze = 6.0;
%maxEccToAnalyze = 6.1;

%minEccToAnalyze = 12.5;
%maxEccToAnalyze = 12.6;

%minEccToAnalyze = 20.0;
%maxEccToAnalyze = 20.1;


%minEccToAnalyze = 2.9;
%maxEccToAnalyze = 3.0;

% Rs/Rc: spread a lot beween 8 - 17
%minEccToAnalyze = 4;
%maxEccToAnalyze = 5;

% Rs/Rc: spread a lot beween 8 - 17
%minEccToAnalyze = 4.5;
%maxEccToAnalyze = 4.5;

% Rs/Rc: spread a lot beween 5 - 19
%minEccToAnalyze = 5;
%maxEccToAnalyze = 5.0;

% Rs/Rc: focused high: 17
%minEccToAnalyze = 5.5;
%maxEccToAnalyze = 5.5;

% Rs/Rc: spread a lot beween 5 - 19
%minEccToAnalyze = 6;
%maxEccToAnalyze = 6.0;

% Rs/Rc: spread a lot beween 2 - 19
%minEccToAnalyze = 7;
%maxEccToAnalyze = 7.0;

% Rs/Rc: spread beween 5 - 17, but with some focus towards 7
%minEccToAnalyze = 8;
%maxEccToAnalyze = 8.0;

% Rs/Rc: spread beween 5 - 17, but with some focus towards 7
%minEccToAnalyze = 9;
%maxEccToAnalyze = 9.0;

% Rs/Rc: focused beween 8 - 9
%minEccToAnalyze = 10;
%maxEccToAnalyze = 10.0;

% Rs/Rc: focused beween 8 - 9
%minEccToAnalyze = 11;
%maxEccToAnalyze = 11.0;


%minEccToAnalyze = 2.9;
%maxEccToAnalyze = 5.1;

%minEccToAnalyze = 0.01;
%maxEccToAnalyze = 1;



% If comparing to Wool et al, run their specific ecc ranges
if (compareToWoolEtAl)
   
    if (visualizeSingleRGCConeIsolatingSTFsForOneMosaic)
        % Visualize STFs for a single mosaic
        eccRangeMMsInMacaqueRetina = singleRGCConeIsolatingSTFeccMMs*[1 1];
        eccRangeDegsInMacaqueRetina = RGCMosaicConstructor.helper.convert.eccentricityInMacaqueRetina('MMsToDegs', eccRangeMMsInMacaqueRetina );
        [eccDegsOut, eccMMsOut] = RGCMosaicConstructor.helper.convert.eccentricityDegsBetweenSpecies('MacaqueRetinaToHumanRetina', eccRangeDegsInMacaqueRetina);
        
        minEccToAnalyze = eccDegsOut(1)
        maxEccToAnalyze = eccDegsOut(2)
    else
        % Aggregated runs
        eccRangeMMsInMacaqueRetina = [0.05 8.5]; % Full RANGE
        eccRangeMMsInMacaqueRetina = [0.01 8];  % Selected RANGE
        eccRangeDegsInMacaqueRetina = RGCMosaicConstructor.helper.convert.eccentricityInMacaqueRetina('MMsToDegs', eccRangeMMsInMacaqueRetina );
        [eccDegsOut, eccMMsOut] = RGCMosaicConstructor.helper.convert.eccentricityDegsBetweenSpecies('MacaqueRetinaToHumanRetina', eccRangeDegsInMacaqueRetina)
        
        minEccToAnalyze = min(eccDegsOut);
        maxEccToAnalyze = max(eccDegsOut);
    end

end


% Split the computation between even and odd runs so we can run the 2 sets in different MATLAB instances
% Choose one of the options below
runOdd = true;
runEven = runOdd;

% targetVisualSTFdescriptorToOptimizeFor: choose from
%  - 'default';              % Mean of Croner & Kaplan with default Rc/Rs
%  - 'x1.5 RsRcRatio';       % 50% higher than Croner & Kaplan
%  - 'higher RsRcRatio';
%  - 'higher intSCratio';
%  - 'very high intSCratio';
%  - 'lower intSCratio';
%  - 'ultra high intSCratio';
targetVisualSTFdescriptorToOptimizeFor = 'default';

% Generate params struct
pStruct = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
            coneMosaicSpecies, opticsSubjectName, RGCmosaicName, targetVisualSTFdescriptorToOptimizeFor);


% Choose chromaticity for STF responses
chromaticityForSTFresponses = 'Achromatic';   % Choose between 'LconeIsolating', 'MconeIsolating', 'Achromatic'
coneFundamentalsOptimizedForStimPosition = false;

if (compareToWoolEtAl)
    % Select which cone isolating STF to run
    chromaticityForSTFresponses = 'LconeIsolating';
    %chromaticityForSTFresponses = 'MconeIsolating';

    coneFundamentalsOptimizedForStimPosition = true;
end


% Generate filenames
[theSurroundConnectedMRGCMosaicFullFileName, ...
 theInputConeMosaicSTFResponsesFullFileName, ...
 theMRGCMosaicSTFResponsesFullFileName, ...
 surroundConnectedParamsStruct, ...
 targetVisualSTFmodifierStruct] = RGCMosaicAnalyzer.filepathFor.testRuns(...
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
            opticsForSTFresponses, residualWithRespectToNativeOpticsDefocusDiopters, ...
            chromaticityForSTFresponses, ...
            coneFundamentalsOptimizedForStimPosition, ...
            'STFresponses');


theOriginalSurroundConnectedMRGCMosaicFullFileName = theSurroundConnectedMRGCMosaicFullFileName;
theOriginalMRGCMosaicSTFResponsesFullFileName = theMRGCMosaicSTFResponsesFullFileName;
theOriginalInputConeMosaicSTFResponsesFullFileName = theInputConeMosaicSTFResponsesFullFileName;

targetString = sprintf('Ecc%2.1f_%2.1f_Size%2.1fx%2.1f',...
    pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs(1), ...
    pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs(2), ...
    pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs(1), ...
    pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs(2));

if (~isempty(strfind(theOriginalMRGCMosaicSTFResponsesFullFileName, 'UpperH1paramsNarrowVisualSTFparamTolerance')))
    targetSurroundOptimizationStrategy = 'UpperH1paramsNarrowVisualSTFparamTolerance';
end

if (~isempty(strfind(theOriginalMRGCMosaicSTFResponsesFullFileName, 'MidH1paramsNarrowVisualSTFparamTolerance')))
    targetSurroundOptimizationStrategy  = 'MidH1paramsNarrowVisualSTFparamTolerance';
end

if (~isempty(strfind(theOriginalMRGCMosaicSTFResponsesFullFileName, 'LowH1paramsNarrowVisualSTFparamTolerance')))
    targetSurroundOptimizationStrategy  = 'LowH1paramsNarrowVisualSTFparamTolerance';
end

% Infos for all runs
[surroundOptimizationStrategyToBeAggregated, ...
     targetNumberOfMappedCellsToBeAggregated, ...
     targetHorizontalEccentricitiesToBeAggregated, ...
     targetMosaicSizesToBeAggregated, ...
     targetMappedPositionDegsToBeAggregated ] = testsHelperFunctions.assembleEccentricityMatrices();

if (performCronerKaplanAnalysis)
    % Pick the orientation for the STF analyses. 
    % This has an effect only if reAnalyzeSTFData == true
    
    if (retinalSpaceReferredSTFs)
        % No physiological optics, so pick the orientation for which the STF attains its max response
        % Doesnt make sense to choose the orientation for which the STF at half max extends to the highest
        % because we have perfect optics here, so no anisotropy as we could
        % have under physiological optics
        fixedOptimalOrientation = 'OptimalOrientation';
    else
        % Orientation for which the STF at half max extends to the highest
        % SF (presumably the orientation that corresponds to the highest
        % resolution phsyiological optics PSF)
        fixedOptimalOrientation = [];
    end

    % Target RGCs with surround specificity around 0.5 (i.e., 50/50 L/M cone net weight in the surround)
    targetedSurroundPurityRange = [];           % [0.4 0.6];
    targetedRadialEccentricityRange = [];       % Empty means full range
    targetedCenterConeNumerosityRange = [];     % single cone centers only
    targetedCenterPurityRange = [];
end % performCronerKaplanAnalysis

if (compareToWoolEtAl)
    % Orientation for which the STF attains its max response
    fixedOptimalOrientation = 'OptimalOrientation';

    % Target RGCs with surround specificity around 0.5 (i.e., 50/50 L/M cone net weight in the surround)
    targetedSurroundPurityRange = [];           % [0.4 0.6];
    targetedRadialEccentricityRange = [];       % Empty means full range
    targetedCenterConeNumerosityRange = [];     % single cone centers only
    targetedCenterPurityRange = [];
end

    
% Determine which runInidices to compute
runIndicesToCompute = find(...
    (abs(targetMappedPositionDegsToBeAggregated(:)) >= minEccToAnalyze) & ...
    (abs(targetMappedPositionDegsToBeAggregated(:)) <= maxEccToAnalyze));

if (isempty(runIndicesToCompute))
    % Find the closest ecc
    [~, runIndicesToCompute] = min(abs(abs(targetMappedPositionDegsToBeAggregated(:)) - minEccToAnalyze));
end

if (isempty(runIndicesToCompute))
    fprintf('Did not find any runs for mosaics with eccentricities in the range: [%2.2f - %2.2f]\n', ...
        minEccToAnalyze, maxEccToAnalyze);
    return;
else
    fprintf('Found %d runs for mosaics with eccentricities in the range: [%2.2f - %2.2f]\n', ...
        numel(runIndicesToCompute), minEccToAnalyze, maxEccToAnalyze);
end

%runIndicesToCompute = [2 10 25 32]


if (runOdd&&~runEven)
    idx = 1:2:numel(runIndicesToCompute);
    fprintf('Will run odd indices');
elseif (runEven&&~runOdd)
    idx = 2:2:numel(runIndicesToCompute);
    fprintf('Will run even indices');
else
    idx = 1:numel(runIndicesToCompute);
    fprintf('Will run all indices');
end


runIndicesToCompute = runIndicesToCompute(idx)


% Filenames
allMRGCMosaicSTFResponsesFullFileNames = cell(1, numel(runIndicesToCompute));
allCronerKaplanAnalysisFileNames = cell(1, numel(runIndicesToCompute));
allWoolAnalysisFileNames = cell(1, numel(runIndicesToCompute));

% Run only the specified run indiecs
for iRun = 1:numel(runIndicesToCompute)
    % Select the run index
    theRunIndexToCompute = runIndicesToCompute(iRun);

    mosaicEccDegs(1) = targetHorizontalEccentricitiesToBeAggregated(theRunIndexToCompute);
    mosaicEccDegs(2) = pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs(2);
    aggregatedEccString = sprintf('Ecc%2.1f_%2.1f_Size%2.1fx%2.1f',mosaicEccDegs(1), mosaicEccDegs(2), ...
                    targetMosaicSizesToBeAggregated(theRunIndexToCompute,1),  targetMosaicSizesToBeAggregated(theRunIndexToCompute,1));
    
    theSurroundConnectedMRGCMosaicFullFileName = strrep(...
                strrep(theOriginalSurroundConnectedMRGCMosaicFullFileName, targetString, aggregatedEccString), ...
                targetSurroundOptimizationStrategy, surroundOptimizationStrategyToBeAggregated{theRunIndexToCompute});
    
    theMRGCMosaicSTFResponsesFullFileName = strrep(...
                strrep(theOriginalMRGCMosaicSTFResponsesFullFileName, targetString, aggregatedEccString), ...
                targetSurroundOptimizationStrategy, surroundOptimizationStrategyToBeAggregated{theRunIndexToCompute});
        
    theInputConeMosaicSTFResponsesFullFileName = strrep(...
                strrep(theOriginalInputConeMosaicSTFResponsesFullFileName, targetString, aggregatedEccString), ...
                targetSurroundOptimizationStrategy, surroundOptimizationStrategyToBeAggregated{theRunIndexToCompute});
    
    % Load theMRGCmosaic
    load(theSurroundConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');
    
    % Add the mapped position info
    mappedPositionDegs = targetMappedPositionDegsToBeAggregated(theRunIndexToCompute,:);
    mappedSizeDegs = RGCMosaicAnalyzer.compute.mSequenceRFmapSizeDegsForMosaic(theMRGCMosaic, mappedPositionDegs, targetNumberOfMappedCellsToBeAggregated{theRunIndexToCompute});
            
    mappedPositionString = sprintf('_@%2.1f_%2.1f.mat', mappedPositionDegs(1), mappedPositionDegs(2));
    theMRGCMosaicSTFResponsesFullFileName = strrep(theMRGCMosaicSTFResponsesFullFileName, ...
                '.mat', mappedPositionString);

    theInputConeMosaicSTFResponsesFullFileName = strrep(theInputConeMosaicSTFResponsesFullFileName, ...
                '.mat', mappedPositionString);

    if (~isempty(mRGCNonLinearityParams))

        % Update inputConeMosaicSTFresponses filename with nonLinearity info
        nonLinearityPreFix = 'inputConeMosaic';
        switch (mRGCNonLinearityParams.type)
            case 'photocurrent'
                nonLinearityPostFix = 'inputConeMosaicPhotocurrent';
            otherwise
                error('Unknown mRGCNonLinearityParams.type')
        end
        theInputConeMosaicSTFResponsesFullFileName = strrep(theInputConeMosaicSTFResponsesFullFileName, ...
            nonLinearityPreFix, nonLinearityPostFix);


        % Update inputConeMosaicSTFresponses filename with stimulus info
        if (~isempty(customTemporalFrequencyAndContrast))
            stimulusPreFix = chromaticityForSTFresponses;
            stimulusPostFix = sprintf('%s%2.0f_x%1.0f_%02.0fHz', ...
                stimulusPreFix, ...
                customTemporalFrequencyAndContrast.totalContrast*100, ...
                customTemporalFrequencyAndContrast.backgroundLuminanceMultiplier, ...
                customTemporalFrequencyAndContrast.temporalFrequencyHz ...
                );


            theInputConeMosaicSTFResponsesFullFileName = strrep(theInputConeMosaicSTFResponsesFullFileName, ...
                stimulusPreFix, stimulusPostFix);
        end

        % Update mRGCMosaicSTFresponses filename with nonLinearity info
        nonLinearityPreFix = 'mRGCMosaic';
        switch (mRGCNonLinearityParams.type)
            case 'photocurrent'
                nonLinearityPostFix = 'mRGCMosaicPhotocurrent';
            otherwise
                error('Unknown mRGCNonLinearityParams.type')
        end
        theMRGCMosaicSTFResponsesFullFileName = strrep(theMRGCMosaicSTFResponsesFullFileName, ...
            nonLinearityPreFix, nonLinearityPostFix);

         % Update inputConeMosaicSTFresponses filename with stimulus info
        if (~isempty(customTemporalFrequencyAndContrast))
            theMRGCMosaicSTFResponsesFullFileName = strrep(theMRGCMosaicSTFResponsesFullFileName, ...
                stimulusPreFix, stimulusPostFix);
        end



    end


    if (computeInputConeMosaicResponses || computeMRGCMosaicResponses)
        % Crop mosaic to the mapped region
        fprintf('Trimming mRGCMosaic to %2.1f x %2.1f degs at (x,y) = (%2.1f,%2.1f)\n', ...
                mappedSizeDegs(1), mappedSizeDegs(2), mappedPositionDegs(1), mappedPositionDegs(2));
    
        % Extra support for input cone mosaic supporting the RF surrounds
        extraSupportDegsForInputConeMosaic = 0.3*mRGCMosaic.extraSupportDegsForMidgetRGCSurrounds(abs(mappedPositionDegs)+mappedSizeDegs, 0.0);
        visualizeSpatialRelationshipToSourceMosaic = true;
        
        theMRGCMosaic.cropToSizeAtEccentricity(mappedSizeDegs, mappedPositionDegs, ...
                'extraSupportDegsForInputConeMosaic', extraSupportDegsForInputConeMosaic, ...
                'visualizeSpatialRelationshipToSourceMosaic', visualizeSpatialRelationshipToSourceMosaic);



        if (computeInputConeMosaicResponses)
            RGCMosaicAnalyzer.compute.inputConeMosaicSTFsForStimulusChromaticityAndOptics(...
                theMRGCMosaic, opticsForSTFresponses, [], ...
                chromaticityForSTFresponses, ...
                coneFundamentalsOptimizedForStimPosition, ...
                theInputConeMosaicSTFResponsesFullFileName, ...
                theMRGCMosaicSTFResponsesFullFileName, ...
                visualizeInputConeMosaicResponses, ...
                visualizePSFonTopOfConeMosaic, ...
                'mRGCNonLinearityParams', mRGCNonLinearityParams, ...
                'customTemporalFrequencyAndContrast', customTemporalFrequencyAndContrast, ...
                'debugInputConeMosaicPcurrentResponse', debugInputConeMosaicPcurrentResponse);
        end % if (computeInputConeMosaicResponses)
    
        if (computeMRGCMosaicResponses)

            RGCMosaicConstructor.helper.simulateExperiment.spatialTransferFunction(...
                theMRGCMosaic, [], [],...
                theInputConeMosaicSTFResponsesFullFileName, ...
                theMRGCMosaicSTFResponsesFullFileName, ...
                'computeInputConeMosaicResponses', false, ...
                'computeMRGCMosaicResponses', true, ...
                'mRGCNonLinearityParams', mRGCNonLinearityParams, ...
                'visualizeResponse', ~true);
            continue;
        end % if (computeMRGCMosaicResponses)
    end % if (computeInputConeMosaicResponses || computeMRGCMosaicResponses)

    if (compareToWoolEtAl)
        % Just aggregate the filenames
        allMRGCMosaicSTFResponsesFullFileNames{iRun} = theMRGCMosaicSTFResponsesFullFileName;

        allWoolAnalysisFileNames{iRun} = ...
            strrep(theMRGCMosaicSTFResponsesFullFileName, 'STFresponses', 'demos/WoolAnalyses');
    end % if (compareToWoolEtAl)

    if (performCronerKaplanAnalysis)
        % Aggregate the filenames
        allMRGCMosaicSTFResponsesFullFileNames{iRun} = theMRGCMosaicSTFResponsesFullFileName;

        allCronerKaplanAnalysisFileNames{iRun} = ...
            strrep(theMRGCMosaicSTFResponsesFullFileName, 'STFresponses', 'demos/CronerKaplanAnalyses');

        if (reAnalyzeSTFData)
            % Analyze current file
            % Do not aggregate over anything
            aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs = false;
            aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities = false;

            RGCMosaicAnalyzer.compute.CronerAndKaplanSTFanalysis(...
                allMRGCMosaicSTFResponsesFullFileNames{iRun}, ...
                allCronerKaplanAnalysisFileNames{iRun}, ...
                aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
                aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities, ...
                targetedSurroundPurityRange, ...
                targetedRadialEccentricityRange, ...
                targetedCenterConeNumerosityRange, ...
                targetedCenterPurityRange, ...
                reAnalyzeSTFData, ...
                'onlyDisplayCronerKaplanData', ~true, ...
                'fixedOptimalOrientation', fixedOptimalOrientation, ...
                'deltaThresholdForLimitingFittedSTFtoPrimaryPeak', 0.01, ...
                'limitVisualizedCKdataToTheEccentricititesOfSyntheticCells', ~true, ...
                'visualizeFullAndMaximalExcursionSTF', visualizeFullAndMaximalExcursionSTF, ...
                'visualizeSTFfits', visualizeSTFfits, ...
                'visualizeModelFitting', visualizeModelFitting, ...
                'visualizeSTFwithConeWeightsMap', visualizeSTFwithConeWeightsMap, ...
                'retinalSpaceReferredSTFs', retinalSpaceReferredSTFs, ...
                'mRGCNonLinearityParams', mRGCNonLinearityParams, ...
                'showComponentLineWeightingFunctions', ~false, ...
                'customTemporalFrequencyAndContrast', customTemporalFrequencyAndContrast);

        else
            % Do not reanalyze data, just visualize 
            if (onlyVisualizeFittedDOGmodels) || (visualizeSTFwithConeWeightsMap)
                % Analyze current file
                % Do not aggregate over anything
                aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs = false;
                aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities = false;

                RGCMosaicAnalyzer.compute.CronerAndKaplanSTFanalysis(...
                    allMRGCMosaicSTFResponsesFullFileNames{iRun}, ...
                    allCronerKaplanAnalysisFileNames{iRun}, ...
                    aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
                    aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities, ...
                    targetedSurroundPurityRange, ...
                    targetedRadialEccentricityRange, ...
                    targetedCenterConeNumerosityRange, ...
                    targetedCenterPurityRange, ...
                    reAnalyzeSTFData, ...
                    'onlyVisualizeFittedDOGmodels', onlyVisualizeFittedDOGmodels, ...
                    'onlyDisplayCronerKaplanData', ~true, ...
                    'fixedOptimalOrientation', fixedOptimalOrientation, ...
                    'deltaThresholdForLimitingFittedSTFtoPrimaryPeak', 0.01, ...
                    'limitVisualizedCKdataToTheEccentricititesOfSyntheticCells', ~true, ...
                    'visualizeSTFwithConeWeightsMap', visualizeSTFwithConeWeightsMap, ...
                    'retinalSpaceReferredSTFs', retinalSpaceReferredSTFs, ...
                    'mRGCNonLinearityParams', mRGCNonLinearityParams, ...
                    'showComponentLineWeightingFunctions', ~false, ...
                    'customTemporalFrequencyAndContrast', customTemporalFrequencyAndContrast);

                fprintf('Finished visualizing iRun %d\n', iRun);
            end
        end
    end % if (performCronerKaplanAnalysis)

end % for iRun

if (computeInputConeMosaicResponses) || (computeMRGCMosaicResponses)
    return;
end


if (compareToWoolEtAl)

    % What to considered as 'low SFs' (< 0.5 cpd)
    theLowSF = [0 0.025];

    % Use full surround, not just surround-center
    computeSurroundConePurityExcludingCenterConnectedCones = false;

    % Classify L-center vs M-center based on L-- and M--cone isolating response at low SF (What Wool et al did.)
    classifyCenterDominanceBasedOnLowFrequencyResponse = true;
    
    % Classify L-center vs M-center based on RF center weights
    classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights = false;

    maxResponseStrengthRatio = 1.0;

    % Options for visualization of STF phase
    unwrapSTFphase = true;
    unwrapLMopponentSTFphase = true;
    employNativeSTFphaseUnwrapMethod = true;
    forceBandPassSTFphaseUnwrap = false;

    RGCMosaicAnalyzer.compute.WoolSTFanalysis(...
        allMRGCMosaicSTFResponsesFullFileNames, ...
        allWoolAnalysisFileNames, ...
        theLowSF, ...
        targetedSurroundPurityRange, ...
        targetedRadialEccentricityRange, ...
        targetedCenterConeNumerosityRange, ...
        targetedCenterPurityRange, ...
        'classifyCenterDominanceBasedOnLowFrequencyResponse', classifyCenterDominanceBasedOnLowFrequencyResponse, ...
        'classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights', classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights, ...
        'fixedOptimalOrientation', fixedOptimalOrientation, ...
        'singleRGCindexToAnalyze', singleRGCindexToAnalyze, ...
        'computeSurroundConePurityExcludingCenterConnectedCones', computeSurroundConePurityExcludingCenterConnectedCones, ...
        'visualizeSingleRGCConeIsolatingSTFs', visualizeSingleRGCConeIsolatingSTFsForOneMosaic, ...
        'unwrapSTFphase', unwrapSTFphase, ...
        'unwrapLMopponentSTFphase', unwrapLMopponentSTFphase, ...
        'forceBandPassSTFphaseUnwrap', forceBandPassSTFphaseUnwrap, ...
        'employNativeSTFphaseUnwrapMethod', employNativeSTFphaseUnwrapMethod, ...
        'maxResponseStrengthRatio', maxResponseStrengthRatio, ...
        'reAnalyzeData', reAnalyzeSTFData);
end

if (performCronerKaplanAnalysis) && (~reAnalyzeSTFData) && (~onlyVisualizeFittedDOGmodels)
    % Aggregate over eccentricities
    aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs = false;
    aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities = true;


    % Generate summary, over all eccentricities, filename
    theSummaryMatFile = allCronerKaplanAnalysisFileNames{1}
    idx1 = strfind(theSummaryMatFile, '@');
    idx2 = strfind(theSummaryMatFile, '.mat');
    substringToBeReplaced = theSummaryMatFile(idx1:idx2-1);
    theSummaryMatFile = strrep(theSummaryMatFile, substringToBeReplaced, 'AllEccentricities');
    targetMosaicString = 'MRGCMosaic_RE';
    idx1 = strfind(theSummaryMatFile,targetMosaicString);
    if (isempty(idx1))
        targetMosaicString = 'MRGCMosaic_LE';
        idx1 = strfind(theSummaryMatFile, targetMosaicString);
    end
    idx2 = strfind(theSummaryMatFile, '_Phi');
    substringToBeReplaced = theSummaryMatFile(idx1+numel(targetMosaicString):idx2-1);
    theSummaryMatFileName = strrep(theSummaryMatFile, substringToBeReplaced, '');


    [~, ~, ~, ~, radialTemporalEquivalentEccentricityDegsMosaic, ...
     RcDegsMosaic, RsToRcMosaic, intStoCsensMosaic, KcMosaic] = ...
        RGCMosaicAnalyzer.compute.CronerAndKaplanSTFanalysis(...
            allMRGCMosaicSTFResponsesFullFileNames, ...
            allCronerKaplanAnalysisFileNames, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
            aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities, ...
            targetedSurroundPurityRange, ...
            targetedRadialEccentricityRange, ...
            targetedCenterConeNumerosityRange, ...
            targetedCenterPurityRange, ...
            reAnalyzeSTFData, ...
            'onlyDisplayCronerKaplanData', ~true, ...
            'fixedOptimalOrientation', fixedOptimalOrientation, ...
            'showComponentLineWeightingFunctions', ~false, ...
            'limitVisualizedCKdataToTheEccentricititesOfSyntheticCells', ~true);

    save(theSummaryMatFileName, ...
        'radialTemporalEquivalentEccentricityDegsMosaic', ...
        'RcDegsMosaic', 'RsToRcMosaic', 'intStoCsensMosaic', 'KcMosaic');

    fprintf('\n\nSummary data across all analyzed eccentricities saved to\n%s\n\n', theSummaryMatFileName);

end % if (performCronerKaplanAnalysis)



