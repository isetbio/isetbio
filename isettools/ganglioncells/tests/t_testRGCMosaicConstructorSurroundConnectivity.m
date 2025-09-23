% Script to derive optimized surrounds for center-connected mRGC mosaics for either human or macaque retinas
% This script is also used to generate materials for the surround-connectivity
% figures of our PLOS2024 paper
%
% Usage:
%{
    t_testRGCMosaicConstructorSurroundConnectivity     
%}

% Initialize session
close all; clear all;

% Select major mosaic params:  
% - species, 
% - rgcMosaicName (eccentricity),
% - optics subject, and 
% - targetVisualSTFdescriptorToOptimizeFor

% Valid eccentricity strings are:
% rgcMosaicName: choose between {...
%   'VSS2024TalkTemporal3DegsMosaic' ...
%	'PLOSpaperFovealMosaic' ...
%	'PLOSpaperTemporal2DegsMosaic' ...
%	'PLOSpaperTemporal4DegsMosaic' ...
%	'PLOSpaperTemporal7DegsMosaic' ...
%	'PLOSpaperTemporal10DegsMosaic' ...
%	'PLOSpaperTemporal14DegsMosaic' ...
%	'PLOSpaperTemporal19DegsMosaic' ...
%	'PLOSpaperTemporal25DegsMosaic' ...
%	'PLOSpaperTemporal32DegsMosaic' ...
%   }

% Valid targetVisualSTFdescriptorToOptimizeFor strings are:
% 'default';   			  	% mean CK mean RsRcRatio, mean CK intSCratio (what was used in PLOSpaper)
% 'x1.5 RsRcRatio';       	% 1.5 x mean CK mean RsRcRatio, mean CK intSCratio
% 'x1.3 RsRcRatio';       	% 1.3 x mean CK mean RsRcRatio, mean CK intSCratio
% 'x0.75 RsRcRatio';       	% 0.75 x mean CK mean RsRcRatio, mean CK intSCratio
% 'higher RsRcRatio';     	% 1.3 x mean CK mean RsRcRatio, mean CK intSCratio
% 'higher intSCratio';	    % mean CK mean RsRcRatio, 1.3 x mean CK intSCratio
% 'very high intSCratio';   % mean CK mean RsRcRatio, 1.7 x mean CK intSCratio
% 'ultra high intSCratio';  % mean CK mean RsRcRatio, 2.0 x mean CK intSCratio
% 'lower intSCratio';       % mean CK mean RsRcRatio, 0.75 x mean CK intSCratio

targetVisualSTFdescriptorToOptimizeFor = 'default';
%targetVisualSTFdescriptorToOptimizeFor = 'x1.3 RsRcRatio';

% Valid optics subjects name are:
% 'PLOSpaperDefaultSubject'
% 'VSS2024TalkFirstSubject'
%v'VSS2024TalkSecondSubject'

if (1==1)
	% Human mosaic 
	coneMosaicSpecies = 'human';

	% Eccentricity
	rgcMosaicName = 'PLOSpaperFovealMosaic';

	% Run with mRGCMosaic.amplificationInCenterOnlySensitivityCausedByInactiveSurrounds set to 0.7 (from 0.85)
	rgcMosaicName = 'PLOSpaperTemporal32DegsMosaicLowerOverlap';
	rgcMosaicName = 'PLOSpaperTemporal25DegsMosaicLowerOverlap';
	rgcMosaicName = 'PLOSpaperTemporal19DegsMosaicLowerOverlap';
    rgcMosaicName = 'PLOSpaperTemporal14DegsMosaicLowerOverlap';
    rgcMosaicName = 'PLOSpaperTemporal10DegsMosaicLowerOverlap';
    rgcMosaicName = 'PLOSpaperTemporal7DegsMosaicLowerOverlap';
    rgcMosaicName = 'PLOSpaperTemporal4DegsMosaicLowerOverlap';
    rgcMosaicName = 'PLOSpaperTemporal2DegsMosaicLowerOverlap';
    rgcMosaicName = 'PLOSpaperFovealMosaicLowerOverlap';

    % Run with
    % mRGCMosaic.amplificationInCenterOnlySensitivityCausedByInactiveSurrounds 0.6 
    % (best agreement in RF center overlap with Gauthier)
    %rgcMosaicName = 'PLOSpaperTemporal32DegsMosaic';   % Fully Done. Also done with second  subject optics
    %rgcMosaicName = 'PLOSpaperTemporal25DegsMosaic';  % Fully Done. Also done with second  subject optics
    %rgcMosaicName = 'PLOSpaperTemporal19DegsMosaic';  % Fully done. Also done with second  subject optics
    %rgcMosaicName = 'PLOSpaperTemporal14DegsMosaic';  % Fully done. Also done with second  subject optics
    %rgcMosaicName = 'PLOSpaperTemporal10DegsMosaic';  % Fully done. Also done with second  subject optics
    %rgcMosaicName = 'PLOSpaperTemporal7DegsMosaic';   % Fully done. Also done with second  subject optics
    %rgcMosaicName = 'PLOSpaperTemporal4DegsMosaic';   % Fully done . Also done with second  subject optics
    rgcMosaicName = 'PLOSpaperTemporal2DegsMosaic';   % Fully done. DOing second on Crete now
    %rgcMosaicName = 'PLOSpaperFovealMosaic';          % Fully done. Also done with second  subject optics

	% Optics
	opticsSubjectName = 'PLOSpaperDefaultSubject';  % This is subject #2 (ranked #3)
    opticsSubjectName = 'VSS2024TalkFirstSubject';  % This is subject #3 (ranked #7)

	%opticsSubjectName = 'PLOSpaperStrehlRatio_0.21';  % Still very high RcDegs at fovea
	%opticsSubjectName = 'PLOSpaperStrehlRatio_0.19';  % Lets try this one
else
	% Macaque mosaic
	coneMosaicSpecies = 'macaque';

	% Eccentricity
	rgcMosaicName = 'VSS2024TalkTemporal3DegsMosaic';  % For Lee&Shapley analyses
	rgcMosaicName = 'VSS2024TalkTemporal7DegsMosaic';  % For Reid&Shapley analyses (their RFs were in the range 3-13 degs)

	% Optics
	opticsSubjectName = 'VSS2024TalkFirstSubject';
	opticsSubjectName = 'VSS2024TalkSecondSubject';
end



% Initialize RGCMosaic generation params
pStruct = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
    coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptorToOptimizeFor);


% Actions to perform
% STEP1: Compute input cone mosaic responses
performComputeInputConeMosaicSTFresponsesAction = ~true;
% Whether to visualize the input cone mosaic responses during their computation
visualizeInputConeMosaicSTFResponseSequences = ~true;

% STEP 2: Optimize surrounds so as to yield C&K '95 macaque Rs/Rc intS/C ratios
% appropriate for  the mosaic's eccentricity
performOptimizeSurroundConePoolingAction = ~true;

% STEP 3: Inspect aspects of the optimized center/surround RFs
performIspectOptimizationResultsAction = ~true;
if (performIspectOptimizationResultsAction)
	performOptimizeSurroundConePoolingAction = false;
end

% Whether to generate separate figures for each component (better for figures in a paper)
% or a single PDF with all components in a panel array. Summary PDFs have unique names
summaryFigureInsteadOfSeparateFigures = ~true;
% Whether to select a SINGLE optimization file via the finder GUI (true)or a collection of optimization files (false)
guiSelectedSingleOptimizationResult = ~true;

% STEP 4: Generate compute-ready mosaic by interpolating cone weights at the previously optimized positions
% at the vicinity of each and every RGC in the mosaic
performGenerateComputeReadyMosaic = true;
visualizeOptimizationGridOnTopOfMosaic = performGenerateComputeReadyMosaic;
visualizeSurroundConeWeightsInterpolationProcess = ~true;

optimizeAllPositions = true;
if  (performIspectOptimizationResultsAction)&&(guiSelectedSingleOptimizationResult)&&...
	(~performOptimizeSurroundConePoolingAction)&&(~performGenerateComputeReadyMosaic)
		optimizeAllPositions = ~true;
end

if (optimizeAllPositions)
	% Grid of (X,Y)-positions, (W,H)-sizes on which the surround will be optimized 
	optimizationPositionsAndSizesGrids = RGCMosaicConstructor.compute.surroundOptimizationGrid(...
		pStruct.rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme, ...
		pStruct.rgcMosaicSurroundOptimization.minGridSize, ...
		pStruct.rgcMosaicSurroundOptimization.maxGridSize, ...
		pStruct.whichZernikeDataBase, pStruct.whichEye, pStruct.sourceLatticeSizeDegs, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
		'withExtremePositions', pStruct.rgcMosaicSurroundOptimization.addEightExtremePositions);

	% How to divide the computation among different compute sessions that can be run in parallel
    optimizationPositionIndicesToCompute = 1:2:size(optimizationPositionsAndSizesGrids,1) 
    %optimizationPositionIndicesToCompute = 2:2:size(optimizationPositionsAndSizesGrids,1) 

else
	optimizationPositionsAndSizesGrids = [-14.0000   -2.0000    3.0000    3.0000];
	optimizationPositionIndicesToCompute = 1;
end

% If we generate the computeReadyMosaic use all optimization positions
if (performGenerateComputeReadyMosaic) || (performIspectOptimizationResultsAction) || (visualizeOptimizationGridOnTopOfMosaic)
    optimizationPositionIndicesToCompute = 1:size(optimizationPositionsAndSizesGrids,1);
end

% Print the optimization positions and sizes
optimizationPositionsAndSizesGrids

% Generate the surroundRetinalConePoolingModel params struct
surroundRetinalConePoolingModelParamsStruct = ...
	RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateSurroundRetinalConePoolingStruct(...
		pStruct.rgcMosaicSurroundOptimization.optimizationStrategy);


% ACTION !!
% Configure a conservative parpool manager. This gives at least 8 GB RAM/core
ASPPManager = AppleSiliconParPoolManager(10);
%
%ASPPManager = AppleSiliconParPoolManager('half max');
%ASPPManager = AppleSiliconParPoolManager('conservative');

if (performOptimizeSurroundConePoolingAction) || (performIspectOptimizationResultsAction)
	if (performOptimizeSurroundConePoolingAction) 
		% Optimize all numerosities
    	centerConeNumerositiesToOptimize = [];
    	% Query user as to which numerosities to optimize
    	% input('Numerosity of RF center cones for which to optimize surrounds. Enter [] for all numerosities found in the analyzed RGCMosaic patch: ');

    	% Query user as to which cone dominance to optimize, 1 or 2
    	centerConeDominanceToOptimize = input('Cone dominance of RF center for which to optimize surrounds. 1: L-cone dominance, 2:M-cone dominance. Your choice: ');

    	% Initial optimization params source options
		% Choose from {'none', 'default', 'imported exact match' or 'imported closest match'}
		%initialSurroundOptimizationValuesSource = 'imported exact match';	
		initialSurroundOptimizationValuesSource = 'imported closest match';
		%initialSurroundOptimizationValuesSource = 'none';	
		%initialSurroundOptimizationValuesSource = 'skip if previous file exists';
	end

	if (strcmp(surroundRetinalConePoolingModelParamsStruct.name, 'PackerDacey2002H1FixedCellIndex')) 
			fixedH1CellIndex = [];
		    % Query user regarding the H1 cell index to employ
		    [fixedH1CellIndex, inspectedOptimizationResultsTargetRFcenterConesNum, ...
		     inspectedOptimizationResultsTargetRFcenterDominantConeType, ...
		     theOptimizationResultsFileNameToBeInspected] = RGCMosaicConstructor.helper.queryUserFor.fixedH1CellIndex(...
			    (guiSelectedSingleOptimizationResult&&performIspectOptimizationResultsAction), pStruct.rgcMosaic.employRFCenterOverlappingMosaic);

		     % Add to the surroundRetinalConePoolingModelParamsStruct  the fixed H1 cell index to be used
		     surroundRetinalConePoolingModelParamsStruct.fixedH1CellIndex = fixedH1CellIndex;

	    if ((performIspectOptimizationResultsAction)&&(~guiSelectedSingleOptimizationResult))
	    	% Get the (inspectedOptimizationResultsTargetRFcenterConesNum, inspectedOptimizationResultsTargetRFcenterDominantConeType)
			inspectedOptimizationResultsTargetRFcenterConesNum = input('Numerosity of RF center cones for which surrounds were optimized. Enter [] for all numerosities found in the analyzed RGCMosaic patch: ');
	    	inspectedOptimizationResultsTargetRFcenterDominantConeType = input('Cone dominance of RF center for which surrounds were optimized. 1: L-cone dominance, 2:M-cone dominance. Your choice: ');
		end
    else
        if (~performOptimizeSurroundConePoolingAction)
		    if (guiSelectedSingleOptimizationResult)
			    [theFixedH1CellIndex, inspectedOptimizationResultsTargetRFcenterConesNum, inspectedOptimizationResultsTargetRFcenterDominantConeType, theOptimizationResultsFileNameToBeInspected] = ...
				    RGCMosaicConstructor.helper.queryUserFor.fixedH1CellIndex(true, pStruct.rgcMosaic.employRFCenterOverlappingMosaic);
		    end
    
		    if ((performIspectOptimizationResultsAction)&&(~guiSelectedSingleOptimizationResult))
			    % Get the (inspectedOptimizationResultsTargetRFcenterConesNum, inspectedOptimizationResultsTargetRFcenterDominantConeType)
			    inspectedOptimizationResultsTargetRFcenterConesNum = input('Numerosity of RF center cones for which surrounds were optimized. Enter [] for all numerosities found in the analyzed RGCMosaic patch: ');
	    	    inspectedOptimizationResultsTargetRFcenterDominantConeType = input('Cone dominance of RF center for which surrounds were optimized. 1: L-cone dominance, 2:M-cone dominance. Your choice: ');
            end
        end
	end
end

% Generate targetVisualSTFmodifierStruct
targetVisualSTFmodifierStruct = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateTargetVisualSTFmodifiersStruct(...
	pStruct.rgcMosaicSurroundOptimization.targetVisualSTFdescriptor);

% Generate the surroundConnectivity simulation params struct 
surroundConnectivitySimulationParamsStruct = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateSurroundConnectivitySimulationParamsStruct(...
	pStruct.whichZernikeDataBase, pStruct.whichSubjectID, pStruct.rgcMosaic.employRFCenterOverlappingMosaic, ...
	optimizationPositionsAndSizesGrids, surroundRetinalConePoolingModelParamsStruct, ...
	targetVisualSTFmodifierStruct);


if (performComputeInputConeMosaicSTFresponsesAction)
	visualizeSamplingPositionsForUncroppedMosaic = true;
	visualizeSpatialRelationshipToSourceMosaic = true;
	% Compute input cone mosaic STF responses
	RGCMosaicConstructor.compute.poolingFunctionsForSurroundOptimizationGrid(...
 		pStruct.whichEye, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
 		pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
 		pStruct.customLMSconeDensities, ...
 		surroundConnectivitySimulationParamsStruct, ...
 		'optimizationPositionIndicesToCompute', optimizationPositionIndicesToCompute, ...
 		'computeInputConeMosaicResponses', true, ...
		'optimizeSurroundConePooling', false, ...
		'visualizeInputConeMosaicSTFResponseSequences', visualizeInputConeMosaicSTFResponseSequences, ...
		'visualizeSamplingPositionsForUncroppedMosaic', visualizeSamplingPositionsForUncroppedMosaic, ...
    	'visualizeSpatialRelationshipToSourceMosaic', visualizeSpatialRelationshipToSourceMosaic);
    return;
end


if (performIspectOptimizationResultsAction) || (performGenerateComputeReadyMosaic) || (visualizeOptimizationGridOnTopOfMosaic&&(~performGenerateComputeReadyMosaic))

	% Assemble pdfExportsDir
	pdfExportSubDir = 'optResults';

	if ((guiSelectedSingleOptimizationResult)&&(contains(theOptimizationResultsFileNameToBeInspected, 'optResultsOverlap'))) || ...
	   ((~guiSelectedSingleOptimizationResult)&&(pStruct.rgcMosaic.employRFCenterOverlappingMosaic))
    	pdfExportSubDir = 'optResultsOverlap';
	end

	if (performIspectOptimizationResultsAction)
		if (~guiSelectedSingleOptimizationResult)
			% A collection of optimization 
			[theOptimizationResultsFileNameCollectionToBeInspected, theOptimizationResultsTargetCenterConesNumToBeInspected] = ...
				RGCMosaicConstructor.compute.poolingFunctionsForSurroundOptimizationGrid(...
			 		pStruct.whichEye,  ...
			 		pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
					pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
 					pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
			 		pStruct.customLMSconeDensities, ...
			 		surroundConnectivitySimulationParamsStruct, ...
			        'optimizationPositionIndicesToCompute', optimizationPositionIndicesToCompute, ...
			        'centerConeNumerositiesToOptimize', inspectedOptimizationResultsTargetRFcenterConesNum, ...
			        'centerConeDominanceToOptimize', inspectedOptimizationResultsTargetRFcenterDominantConeType, ...
			 		'computeInputConeMosaicResponses', ~true, ...
					'optimizeSurroundConePooling', ~true, ...
					'onlyReturnSurroundOptimizationResultFilenames', true);
		else
			% A single, user-selected optimization file
			theOptimizationResultsFileNameCollectionToBeInspected{1} = theOptimizationResultsFileNameToBeInspected;
			idx1 = strfind(theOptimizationResultsFileNameToBeInspected, 'Dominated_');
			idx2 = strfind(theOptimizationResultsFileNameToBeInspected, 'coneRF');
			theOptimizationResultsTargetCenterConesNumToBeInspected = ...
				str2num(theOptimizationResultsFileNameToBeInspected(idx1+numel('Dominated_'):idx2-1));
		end

        for iOptimizationResultsFileIndex = 1:numel(theOptimizationResultsFileNameCollectionToBeInspected)
			theOptimizationResultsFileNameToBeInspected = theOptimizationResultsFileNameCollectionToBeInspected{iOptimizationResultsFileIndex};
			dataAvailable = false;
			try
				load(theOptimizationResultsFileNameToBeInspected, ...
		   			'theMRGCMosaic', ...
			   		'targetVisualSTFparams', ...
			        'theTargetRGCindex', ...
			        'optimizationResults');
				dataAvailable = true;
			catch 
		   		fprintf('\nThe optimization results file \n\t %s \n was **NOT** found. Skipping inspection.\n', theOptimizationResultsFileNameToBeInspected);
		   		continue
			end

			% Assemble theSummaryPDFfileName
			idx = strfind(theOptimizationResultsFileNameToBeInspected,pdfExportSubDir);
			theSummaryPDFfileName = theOptimizationResultsFileNameToBeInspected(idx+numel(pdfExportSubDir):end);
			theSummaryPDFfileName = strrep(theSummaryPDFfileName, '.mat', '.pdf');

			if (isempty(inspectedOptimizationResultsTargetRFcenterConesNum))
				currentlyInspectedOptimizationResultsTargetRFcenterConesNum = theOptimizationResultsTargetCenterConesNumToBeInspected(iOptimizationResultsFileIndex);
			else
				currentlyInspectedOptimizationResultsTargetRFcenterConesNum = inspectedOptimizationResultsTargetRFcenterConesNum;
			end

			% Use fixed spatial support and scale bar when plotting the cone weights maps and line weighting functions
			% If [] is specified, these are computed automatically 
			if (sqrt(sum(pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs.^2,2)) > 20)
				fixedSpatialSupportTickSeparationArcMinForConePoolingMaps =  60;

				% Match the Field 2010 paper scale bar in Fig. 4 which was 25 microns
				Field2010Figure4CellsScaleBarMicrons = 25;

				% Convert to degs in macaque retina
				Field2010Figure4CellsScaleBarDegs = RGCMosaicConstructor.helper.convert.eccentricityInMacaqueRetina('MMsToDegs', Field2010Figure4CellsScaleBarMicrons * 1e-3);

				% Convert to degs in human retina
				fixedScaleBarDegsForConePoolingMaps = RGCMosaicConstructor.helper.convert.eccentricityDegsBetweenSpecies('MacaqueRetinaToHumanRetina', Field2010Figure4CellsScaleBarDegs);

		    elseif (sqrt(sum(pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs.^2,2)) < 1)
		    	% For foveal mosaics, 2.4 arc min which is 0.04 degs
		    	fixedSpatialSupportTickSeparationArcMinForConePoolingMaps =  2.4;
				fixedScaleBarDegsForConePoolingMaps =  0.05;

			else
				% 6 arc min which is 0.1 degs
				fixedSpatialSupportTickSeparationArcMinForConePoolingMaps =  6;
				fixedScaleBarDegsForConePoolingMaps =  0.1;
            end 

			% Visualize imported optimization results and export PDFs of the analyses
			RGCMosaicConstructor.visualize.optimizationResults(optimizationResults, ...
				theMRGCMosaic, theTargetRGCindex, targetVisualSTFparams, ...
				pdfExportSubDir, ...
				currentlyInspectedOptimizationResultsTargetRFcenterConesNum, ...
				inspectedOptimizationResultsTargetRFcenterDominantConeType, ...
				summaryFigureInsteadOfSeparateFigures, theSummaryPDFfileName, ...
				'fixedSpatialSupportTickSeparationArcMin', fixedSpatialSupportTickSeparationArcMinForConePoolingMaps, ...
				'fixedScalaBarDegs', fixedScaleBarDegsForConePoolingMaps, ...
				'contourGenerationMethod', 'ellipseFitToPooledConePositions', ...
				'maxNumberOfConesOutsideContour', 0);
		end % for iOptimizationResultsFileIndex

	else % if (performGenerateComputeReadyMosaic) || (visualizeOptimizationGridOnTopOfMosaic&&(~performGenerateComputeReadyMosaic))

		% Do not introduce additional variance to the intSens surround-to-center ratio
		surroundVarianceInComputeReadyMosaic = struct();

		% Or introduce additional variance and a bias
		if ((~isempty(pStruct.rgcMosaicSurroundOptimization.intSensRatioBias))&&(~isempty(pStruct.rgcMosaicSurroundOptimization.intSensRatioVariance)))
			surroundVarianceInComputeReadyMosaic = struct(...
				'intSensRatioBias', pStruct.rgcMosaicSurroundOptimization.intSensRatioBias, ...
				'intSensRatioSigma', sqrt(pStruct.rgcMosaicSurroundOptimization.intSensRatioVariance));
		end

		hFig = RGCMosaicConstructor.compute.computeReadyMosaic(...
			pStruct.whichEye, ...
			pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
			pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
 			pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
			pStruct.customLMSconeDensities, ...
			surroundConnectivitySimulationParamsStruct, ...
			optimizationPositionIndicesToCompute, ...
 			visualizeOptimizationGridOnTopOfMosaic, ...
            'employLconeDominanceOptimizationOnly', pStruct.rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly, ...
 			'surroundVarianceInComputeReadyMosaic', surroundVarianceInComputeReadyMosaic, ...
 			'onlyVisualizeOptimizationGrid', visualizeOptimizationGridOnTopOfMosaic&&(~performGenerateComputeReadyMosaic), ...
 			'visualizeInterpolationProcess', visualizeSurroundConeWeightsInterpolationProcess);

	end % if (performGenerateComputeReadyMosaic)
end % if (performIspectOptimizationResultsAction)||performGenerateComputeReadyMosaic)

if (performOptimizeSurroundConePoolingAction)
	% The default initial values
	defaultUserSuppliedInitialValuesForSurroundPoolingModel = struct(...
   		'poolingModel', struct(...
   			'name', surroundConnectivitySimulationParamsStruct.poolingOptimizationParamsStruct.poolingModel.name, ...
   			'weightsComputeHandle', surroundConnectivitySimulationParamsStruct.poolingOptimizationParamsStruct.poolingModel.weightsComputeHandle), ...
   		'optimizedValues', struct(...
   			'Kc', 1, ...
   			'RwDegs', 0.258, ...
   			'KsKcRatio', 0.723, ...
   			'VnVwRatio', 0.175, ...
   			'RnRwRatio', 0.227) ...
       		);

	switch (initialSurroundOptimizationValuesSource)
		case {'imported exact match', 'imported closest match'}
			userSuppliedInitialValuesForSurroundPoolingModel = defaultUserSuppliedInitialValuesForSurroundPoolingModel;

			% Try to load optimizationResults from a previously generated optimizationResults file
			optimizationResults = ...
				RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.findPreviouslyOptimizedSurroundPoolingModels(...
					pStruct.whichEye,  ...
			 		pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
					pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
 					pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
			 		pStruct.customLMSconeDensities, ...
			 		surroundConnectivitySimulationParamsStruct, ...
					'targetRFcenterConesNum', centerConeNumerositiesToOptimize, ...
					'targetRFcenterDominantConeType', centerConeDominanceToOptimize, ...
					'initialSurroundOptimizationValuesSource', initialSurroundOptimizationValuesSource, ...
					'optimizationPositionIndicesToCompute', optimizationPositionIndicesToCompute);
			if (~isempty(optimizationResults))
    			% Extract optimizedValues from the loaded optimizationResults
    			userSuppliedInitialValuesForSurroundPoolingModel.optimizedValues = ...		
    				RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.extractPreviouslyOptimizedValues(optimizationResults);
    		end
    		% Try to load an exact optimizationResults file at runtime
    		userSuppliedInitialValuesForSurroundPoolingModel.attemptToLoadExactOptimizationResultsFile = true;
    		userSuppliedInitialValuesForSurroundPoolingModel.skipOptimizationIfOptimizationFileExists = false;

    	case 'skip if previous file exists'
    		userSuppliedInitialValuesForSurroundPoolingModel = defaultUserSuppliedInitialValuesForSurroundPoolingModel;

			% Try to load optimizationResults from a previously generated optimizationResults file
			optimizationResults = ...
				RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.findPreviouslyOptimizedSurroundPoolingModels(...
					pStruct.whichEye,  ...
			 		pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
					pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
 					pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
			 		pStruct.customLMSconeDensities, ...
			 		surroundConnectivitySimulationParamsStruct, ...
					'targetRFcenterConesNum', centerConeNumerositiesToOptimize, ...
					'targetRFcenterDominantConeType', centerConeDominanceToOptimize, ...
					'initialSurroundOptimizationValuesSource', 'imported exact match', ...
					'optimizationPositionIndicesToCompute', optimizationPositionIndicesToCompute);
			if (~isempty(optimizationResults))
    			% Extract optimizedValues from the loaded optimizationResults
    			userSuppliedInitialValuesForSurroundPoolingModel.optimizedValues = ...		
    				RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.extractPreviouslyOptimizedValues(optimizationResults);
    		end

    		% Try to load an exact optimizationResults file at runtime
    		userSuppliedInitialValuesForSurroundPoolingModel.attemptToLoadExactOptimizationResultsFile = true;
    		userSuppliedInitialValuesForSurroundPoolingModel.skipOptimizationIfOptimizationFileExists = true;

		case 'default'
			% Generate userSuppliedInitialValuesForSurroundPoolingModel from scratch
	    	userSuppliedInitialValuesForSurroundPoolingModel = defaultUserSuppliedInitialValuesForSurroundPoolingModel;
	    	userSuppliedInitialValuesForSurroundPoolingModel.skipOptimizationIfOptimizationFileExists = false;

		case 'none'
			userSuppliedInitialValuesForSurroundPoolingModel = [];

		otherwise
			error('Unknown initialSurroundOptimizationValuesSource: ''%s''.\n', initialSurroundOptimizationValuesSource);
	end % switch (customInitialSurroundOptimizationValues)

	% Optimize surround pooling weights !
	RGCMosaicConstructor.compute.poolingFunctionsForSurroundOptimizationGrid(...
		pStruct.whichEye,  ...
		pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
 		pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
		pStruct.customLMSconeDensities, ...
		surroundConnectivitySimulationParamsStruct, ...
        'optimizationPositionIndicesToCompute', optimizationPositionIndicesToCompute, ...
        'centerConeNumerositiesToOptimize', centerConeNumerositiesToOptimize, ...
        'centerConeDominanceToOptimize',centerConeDominanceToOptimize, ...
 		'computeInputConeMosaicResponses', ~true, ...
		'optimizeSurroundConePooling', true, ...
		'doNotWorryAboutMaximizingTargetRFcenterLMconeRatio', pStruct.rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly, ...
		'userSuppliedInitialValuesForModelVariables', userSuppliedInitialValuesForSurroundPoolingModel, ...
		'visualizeFullAndMaximalExcursionSTF', true, ...
		'visualizeGaussianFitToCenterSTF', true, ...
		'onlyReturnSurroundOptimizationResultFilenames', false);
end % (performOptimizeSurroundConePoolingAction)