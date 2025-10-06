function t_generateMRGCmosaicAtStage3(options)
% Generate an mRGC mosaic at different stages of cone-to-mRGC RF surround connectivity
%
% Syntax:
%   t_generateMRGCmosaicAtStage3;
%
% Description:
%   Demonstrates how to generate an mRGC mosaic at stage 3A, 3B, or 3C of connectivity
%   At stage 3A cones are connected to RGC RF centers in a
%   mutually-exclusive way (no RF center overlap). It is at this stage that
%   the user can specity the desired spatialChromaticUniformityTradeoff
%   At stage 2C cone connections diverge to nearby RGC RF centers
%   generating overlap between neighboring RF centers
%
%  This is set up with key/value pairs that demonstate how to select different
%  options. Different choices are illustrated in the examples
%  in the source code.
%
% Optional key/value pairs
%    See source code arguments block for a list of key/value pairs.

% History:
%    08/28/25  NPC  Wrote it.

% Examples:
%{
    % Compute input cone mosaic STF responses (stage 3A)
    t_generateMRGCmosaicAtStage3(...
        'rgcMosaicName', 'PLOSpaperNasal2DegsTinyMosaic', ...
        'regenerateMosaicAtStage3A', true);

    % Derive optimized surround pooling functions for L-cone dominated mRGCs (stage 3B)
    t_generateMRGCmosaicAtStage3(...
        'rgcMosaicName', 'PLOSpaperNasal2DegsTinyMosaic', ...
        'regenerateMosaicAtStage3B', true, ...
        'centerConeDominanceToOptimize', cMosaic.LCONE_ID);

    % Derive optimized surround pooling functions for M-cone dominated mRGCs (stage 3B)
    t_generateMRGCmosaicAtStage3(...
        'rgcMosaicName', 'PLOSpaperNasal2DegsTinyMosaic', ...
        'regenerateMosaicAtStage3B', true, ...
        'centerConeDominanceToOptimize', cMosaic.MCONE_ID);

    % Inspect derived surround pooling functions for L-cone dominated mRGCs
    t_generateMRGCmosaicAtStage3(...
        'rgcMosaicName', 'PLOSpaperNasal2DegsTinyMosaic', ...
        'inspectMosaicAtStage3B', true, ...
        'centerConeDominanceToInspect', cMosaic.LCONE_ID);

    % Inspect derived surround pooling functions for M-cone dominated mRGCs
    t_generateMRGCmosaicAtStage3(...
        'rgcMosaicName', 'PLOSpaperNasal2DegsTinyMosaic', ...
        'inspectMosaicAtStage3B', true, ...
        'centerConeDominanceToInspect', cMosaic.MCONE_ID);

    % Inspect the surround pooling interpolation grid
    t_generateMRGCmosaicAtStage3(...
        'rgcMosaicName', 'PLOSpaperNasal2DegsTinyMosaic', ...
        'inspectSurroundPoolingInterpolationGrid', true);

    t_generateMRGCmosaicAtStage3(...
        'rgcMosaicName', 'PLOSpaperNasal2DegsTinyMosaic', ...
         regenerateMosaicAtStage3C', true);

%}


arguments
    % ---- Name encoding properties of the rgcMosaic, such as its eccentricity ---
    % See RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters
    % for what is available and to add new mosaics
    options.rgcMosaicName (1,:) char = 'PLOSpaperNasal7DegsMosaic';

    % ---- Which species to employ ----
    % Choose between {'macaque', 'human'}. If 'macaque' is chosen, the input
    % cone mosaic has a 1:1 L/M cone ratio.
    options.coneMosaicSpecies  (1,:) char {mustBeMember(options.coneMosaicSpecies,{'human','macaque'})} = 'human';


    % ----- Which subject optics to employ -----
    options.opticsSubjectName (1,:) char = 'PLOSpaperDefaultSubject';


    % ------ targetVisualSTF options ----
    % Options are : {'default', 'x1.3 RsRcRatio'}
    % These are with respect to the macaque data of the Croner & Kaplan '95 study
    % 'default': target the mean Rs/Rc, and the mean Ks/Kc (Rs/Rc)^2
    % 'x1.3 RsRcRatio': target Rs/Rc ratio that is 1.3 x mean, and target Ks/Kc (Rs/Rc)^2: mean
    options.targetVisualSTFdescriptorToOptimizeFor (1,:) char = 'default';

    % Which center cone dominance to derive surround cone pooling function
    % for. Choose from {cMosaic.LCONE_ID cMosaic.MCONE_ID}
    options.centerConeDominanceToOptimize (1,1) double = cMosaic.LCONE_ID;

    % Which center cone dominance to inspect surround cone pooling function
    % for. Choose from {cMosaic.LCONE_ID cMosaic.MCONE_ID}
    options.centerConeDominanceToInspect (1,1) double = cMosaic.LCONE_ID;

    % ---- Choices of actions to perform ----
    % Whether to regenerate the mosaic at stage3A 
    % (computation of input cone mosaic STF responses)
    options.regenerateMosaicAtStage3A (1,1) logical = false;

    % Whether to regenerate the mosaic at stage3B 
    % (determine optimized surround cone pooling functions so as to yield C&K '95 macaque Rs/Rc intS/C ratios
    % appropriate for  the mosaic's eccentricity)
    options.regenerateMosaicAtStage3B (1,1) logical = false;

    % Inspect the optimized surround pooling functions
    options.inspectMosaicAtStage3B (1,1) logical = false;

    % Inspect the interpolation grid for the surround pooling computation
    options.inspectSurroundPoolingInterpolationGrid (1,1) logical = false;

    % Whether to regenerate the mosaic at stage3C. This is the compute-ready mosaic
    options.regenerateMosaicAtStage3C (1,1) logical = false;


    % ---- Visualization options ----

    % Stage 3A visualizations (input cone mosaic STF responses)
    % Whether to visualize the optimization of the PSF Strehl ratio
    options.visualizeStrehlRatioOptimization (1,1) logical = false;

    % Whether to visualize the employed PSF
    options.visualizeEmployedPSF (1,1) logical = false;

    % Whether to visualize the input cone mosaic responses during their computation
    options.visualizeInputConeMosaicSTFResponseSequences (1,1) logical = false;

    options.visualizeSamplingPositionsForUncroppedMosaic (1,1) logical = false;
	options.visualizeSpatialRelationshipToSourceMosaic (1,1) logical = false;

    % Stage 3B visualizations (surround pooling function derivation)
    options.visualizeFullAndMaximalExcursionSTF (1,1) logical = false;
    options.visualizeGaussianFitToCenterSTF (1,1) logical = false;

    % Stage 3C visualizations (surround pooling function interpolation)
    options.visualizeOptimizationGridOnTopOfMosaic (1,1) logical = false;
    options.visualizeSurroundConeWeightsInterpolationProcess (1,1) logical = false;

    % Whether to generate separate figures for each component (better for figures in a paper)
    % or a single PDF with all components in a panel array. Summary PDFs have
    % unique names for each visualized location.
    % Non-summary PDFs get overriden for each visualized location (but the code
    % pauses at each location). 
    % This flag has an effect only for Stage 3B (inspection of optimized
    % surround cone pooling functions)
    options.summaryInsteadOfSeparateInspectionFigures (1,1) logical = true;

    % Whether to close previously open figures
    options.closeOpenFigures (1,1) logical = true;

end  % arguments

% Set flags from key/value pairs
rgcMosaicName = options.rgcMosaicName;
coneMosaicSpecies = options.coneMosaicSpecies;
opticsSubjectName = options.opticsSubjectName;


% Center cone dominance for which to derive surround pooling functions
centerConeDominanceToOptimize = options.centerConeDominanceToOptimize;

% Center cone dominance for which to inspect the derived surround pooling functions
centerConeDominanceToInspect = options.centerConeDominanceToInspect;

% Target STF descriptor
targetVisualSTFdescriptorToOptimizeFor = options.targetVisualSTFdescriptorToOptimizeFor;


% Actions to perform
regenerateMosaicAtStage3A = options.regenerateMosaicAtStage3A;
regenerateMosaicAtStage3B = options.regenerateMosaicAtStage3B;
inspectMosaicAtStage3B = options.inspectMosaicAtStage3B;
inspectSurroundPoolingInterpolationGrid = options.inspectSurroundPoolingInterpolationGrid;
regenerateMosaicAtStage3C = options.regenerateMosaicAtStage3C;

% Visualization options
% Stage 3A visualizations
visualizeEmployedPSF = options.visualizeEmployedPSF;
visualizeStrehlRatioOptimization = options.visualizeStrehlRatioOptimization;
visualizeInputConeMosaicSTFResponseSequences = options.visualizeInputConeMosaicSTFResponseSequences;
visualizeSamplingPositionsForUncroppedMosaic = options.visualizeSamplingPositionsForUncroppedMosaic;
visualizeSpatialRelationshipToSourceMosaic = options.visualizeSpatialRelationshipToSourceMosaic;

% Stage 3B visualizations (surround pooling function derivation)
visualizeFullAndMaximalExcursionSTF = options.visualizeFullAndMaximalExcursionSTF;
visualizeGaussianFitToCenterSTF = options.visualizeGaussianFitToCenterSTF;

% Stage 3B inspection visualization
summaryInsteadOfSeparateInspectionFigures = options.summaryInsteadOfSeparateInspectionFigures;

% Stage 3C visualizations
visualizeOptimizationGridOnTopOfMosaic = options.visualizeOptimizationGridOnTopOfMosaic;
visualizeSurroundConeWeightsInterpolationProcess = options.visualizeSurroundConeWeightsInterpolationProcess;

% Close previously open figures
closePreviouslyOpenFigures = options.closeOpenFigures;

if (closePreviouslyOpenFigures)
    % Close any stray figs
    close all;
end


% Generate the necessary mosaic params struct
pStruct = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
    coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptorToOptimizeFor);


% Generate spatial grid covering the extent of the synthesized mRGC. 
% Surround pooling functions will be derived at each node (position) of this grid
optimizationPositionsAndSizesGrids = RGCMosaicConstructor.compute.surroundOptimizationGrid(...
		pStruct.rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme, ...
		pStruct.rgcMosaicSurroundOptimization.minGridSize, ...
		pStruct.rgcMosaicSurroundOptimization.maxGridSize, ...
		pStruct.whichZernikeDataBase, pStruct.whichEye, pStruct.sourceLatticeSizeDegs, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
		'withExtremePositions', pStruct.rgcMosaicSurroundOptimization.addEightExtremePositions);



% Compute surround optimization functions for all grid positions 
optimizationPositionIndicesToCompute = 1:size(optimizationPositionsAndSizesGrids,1);

% Generate the surroundRetinalConePoolingModel params struct
surroundRetinalConePoolingModelParamsStruct = ...
	RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateSurroundRetinalConePoolingStruct(...
		pStruct.rgcMosaicSurroundOptimization.optimizationStrategy);


% Generate targetVisualSTFmodifierStruct
targetVisualSTFmodifierStruct = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateTargetVisualSTFmodifiersStruct(...
	    pStruct.rgcMosaicSurroundOptimization.targetVisualSTFdescriptor);


% Generate the surroundConnectivity simulation params struct 
surroundConnectivitySimulationParamsStruct = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateSurroundConnectivitySimulationParamsStruct(...
	    pStruct.whichZernikeDataBase, pStruct.whichSubjectID, pStruct.rgcMosaic.employRFCenterOverlappingMosaic, ...
	    optimizationPositionsAndSizesGrids, surroundRetinalConePoolingModelParamsStruct, ...
	    targetVisualSTFmodifierStruct);


if (regenerateMosaicAtStage3A)

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
        'visualizeEmployedPSF', visualizeEmployedPSF, ...
        'visualizeStrehlRatioOptimization', visualizeStrehlRatioOptimization, ...
		'visualizeInputConeMosaicSTFResponseSequences', visualizeInputConeMosaicSTFResponseSequences, ...
		'visualizeSamplingPositionsForUncroppedMosaic', visualizeSamplingPositionsForUncroppedMosaic, ...
    	'visualizeSpatialRelationshipToSourceMosaic', visualizeSpatialRelationshipToSourceMosaic);
    return;
end  % regenerateMosaicAtStage3A



if (regenerateMosaicAtStage3B) || (inspectMosaicAtStage3B)

    if (regenerateMosaicAtStage3B)
        % Optimize all numerosities
        centerConeNumerositiesToOptimize = [];
    	    
        % Options for loading the nitial optimization params
	    % Choose from {'none', 'default', 'imported exact match' or 'imported closest match'}
	    %initialSurroundOptimizationValuesSource = 'imported exact match';	
	    initialSurroundOptimizationValuesSource = 'imported closest match';
	    % initialSurroundOptimizationValuesSource = 'none';	
	    % initialSurroundOptimizationValuesSource = 'skip if previous file exists';
    end

    % User may select to use params from a fixed H1 cell index
    if (strcmp(surroundRetinalConePoolingModelParamsStruct.name, 'PackerDacey2002H1FixedCellIndex')) 
		% Query user regarding the H1 cell index to employ
		[fixedH1CellIndex, inspectedOptimizationResultsTargetRFcenterConesNum, ...
		     inspectedOptimizationResultsTargetRFcenterDominantConeType, ...
		     theOptimizationResultsFileNameToBeInspected] = RGCMosaicConstructor.helper.queryUserFor.fixedH1CellIndex(...
			    (guiSelectedSingleOptimizationResult&&performIspectOptimizationResultsAction), pStruct.rgcMosaic.employRFCenterOverlappingMosaic);

		% Add to the surroundRetinalConePoolingModelParamsStruct the fixed H1 cell index to be used
		surroundRetinalConePoolingModelParamsStruct.fixedH1CellIndex = fixedH1CellIndex;
    end

    if (inspectMosaicAtStage3B)
	    
        % Empty [] numerosity will inspect surround pooling functions for
        % all numerosities encountered in the synthesized RGCMosaic patch
		inspectedOptimizationResultsTargetRFcenterConesNum = []; 

        % Center cone dominance
	    inspectedOptimizationResultsTargetRFcenterDominantConeType = centerConeDominanceToInspect;

        % Get all the surround pooling optimization file names
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

        % Go through each optimization file and plot a summary of the optimization results
        for iOptimizationResultsFileIndex = 1:numel(theOptimizationResultsFileNameCollectionToBeInspected)

			theOptimizationResultsFileNameToBeInspected = theOptimizationResultsFileNameCollectionToBeInspected{iOptimizationResultsFileIndex};
			
            % See if the file exists
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

            % Specify visualization parameters
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

            [pdfExportSubDir, theSummaryPDFfileName] = RGCMosaicConstructor.filepathFor.optimizationResultsSummaryPDF(...
                theOptimizationResultsFileNameToBeInspected, ...
                'generateMissingSubDirs', true);

            % Visualize imported optimization results and export PDFs of the analyses
			RGCMosaicConstructor.visualize.optimizationResults(optimizationResults, ...
				theMRGCMosaic, theTargetRGCindex, targetVisualSTFparams, ...
				pdfExportSubDir, ...
				theOptimizationResultsTargetCenterConesNumToBeInspected(iOptimizationResultsFileIndex), ...
				inspectedOptimizationResultsTargetRFcenterDominantConeType, ...
				summaryInsteadOfSeparateInspectionFigures, theSummaryPDFfileName, ...
				'fixedSpatialSupportTickSeparationArcMin', fixedSpatialSupportTickSeparationArcMinForConePoolingMaps, ...
				'fixedScalaBarDegs', fixedScaleBarDegsForConePoolingMaps, ...
				'contourGenerationMethod', 'ellipseFitToPooledConePositions', ...
				'maxNumberOfConesOutsideContour', 0);

        end % for iOptimizationResultsFileIndex
        return;
    end % if (inspectMosaicAtStage3B)

    if (regenerateMosaicAtStage3B)
        % The default initial values
	    defaultInitialValuesForSurroundPoolingModel = struct(...
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
			    userSuppliedInitialValuesForSurroundPoolingModel = defaultInitialValuesForSurroundPoolingModel;
    
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
    		    userSuppliedInitialValuesForSurroundPoolingModel = defaultInitialValuesForSurroundPoolingModel;
    
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
	    	    userSuppliedInitialValuesForSurroundPoolingModel = defaultInitialValuesForSurroundPoolingModel;
	    	    userSuppliedInitialValuesForSurroundPoolingModel.skipOptimizationIfOptimizationFileExists = false;
    
		    case 'none'
			    userSuppliedInitialValuesForSurroundPoolingModel = [];
    
		    otherwise
			    error('Unknown initialSurroundOptimizationValuesSource: ''%s''.\n', initialSurroundOptimizationValuesSource);
	    end % switch (initialSurroundOptimizationValuesSource)

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
		    'visualizeFullAndMaximalExcursionSTF', visualizeFullAndMaximalExcursionSTF, ...
		    'visualizeGaussianFitToCenterSTF', visualizeGaussianFitToCenterSTF, ...
		    'onlyReturnSurroundOptimizationResultFilenames', false);

        return;
    end % if (regenerateMosaicAtStage3B)

end %if (regenerateMosaicAtStage3B) || (inspectMosaicAtStage3B)

if (inspectSurroundPoolingInterpolationGrid)
    hFig = RGCMosaicConstructor.compute.computeReadyMosaic(...
		pStruct.whichEye, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicEccDegs, ...
		pStruct.rgcMosaicSurroundOptimization.mosaicSizeDegs, ...
		pStruct.rgcMosaic.spatialChromaticUniformityTradeoff, ...
		pStruct.customLMSconeDensities, ...
		surroundConnectivitySimulationParamsStruct, ...
		optimizationPositionIndicesToCompute, ...
	    true, ...
        'onlyVisualizeOptimizationGrid', true, ...
        'employLconeDominanceOptimizationOnly', pStruct.rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly);
    return;
end

if (regenerateMosaicAtStage3C)
    % Default additional variance is zero
    surroundVarianceInComputeReadyMosaic = struct();

	% Or introduce additional variance and a bias
    if (...
            (~isempty(pStruct.rgcMosaicSurroundOptimization.intSensRatioBias)) && ...
            (~isempty(pStruct.rgcMosaicSurroundOptimization.intSensRatioVariance))...
        )
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
 			'visualizeInterpolationProcess', visualizeSurroundConeWeightsInterpolationProcess);

end % if (regenerateMosaicAtStage3C)

end