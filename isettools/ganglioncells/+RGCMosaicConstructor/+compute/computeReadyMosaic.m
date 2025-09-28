function hFig = computeReadyMosaic(whichEye,  mosaicEccDegs, mosaicSizeDegs, ...
	spatialChromaticUniformityTradeoff, customLMSconeDensities, surroundConnectivitySimulationParamsStruct, ...
	optimizationPositionIndicesToCompute, visualizeOptimizationGridOnTopOfMosaic, varargin)

	 % Parse input
    p = inputParser;
    p.addParameter('surroundVarianceInComputeReadyMosaic', @isstruct);
    p.addParameter('visualizeInterpolationProcess', false, @islogical);
    p.addParameter('employLconeDominanceOptimizationOnly', false, @islogical);
    p.addParameter('onlyVisualizeOptimizationGrid', false, @islogical);
    p.addParameter('visualizeNeighboringOptimizationGridNodesWithLines', false, @islogical);

    p.parse(varargin{:});
	onlyVisualizeOptimizationGrid = p.Results.onlyVisualizeOptimizationGrid;

	visualizeInterpolationProcess = p.Results.visualizeInterpolationProcess;
	surroundVarianceInComputeReadyMosaic = p.Results.surroundVarianceInComputeReadyMosaic;
    employLconeDominanceOptimizationOnly = p.Results.employLconeDominanceOptimizationOnly;
    visualizeNeighboringOptimizationGridNodesWithLines = p.Results.visualizeNeighboringOptimizationGridNodesWithLines;

	% All L-cone dominance optimization files
	[theLconeDominanceOptimizationResultsFileNameCollection, ...
	 theLconeDominanceOptimizationResultsTargetCenterConesNum] = ...
			RGCMosaicConstructor.compute.poolingFunctionsForSurroundOptimizationGrid(...
		 		whichEye,  mosaicEccDegs, mosaicSizeDegs, spatialChromaticUniformityTradeoff , ...
		 		customLMSconeDensities, surroundConnectivitySimulationParamsStruct, ...
		        'optimizationPositionIndicesToCompute', optimizationPositionIndicesToCompute, ...
		        'centerConeNumerositiesToOptimize', [], ...
		        'centerConeDominanceToOptimize', cMosaic.LCONE_ID, ...
		 		'computeInputConeMosaicResponses', ~true, ...
				'optimizeSurroundConePooling', ~true, ...
				'onlyReturnSurroundOptimizationResultFilenames', true);

    % All M-cone dominance optimization files
    if (employLconeDominanceOptimizationOnly)
        theMconeDominanceOptimizationResultsFileNameCollection = theLconeDominanceOptimizationResultsFileNameCollection;
	    theMconeDominanceOptimizationResultsTargetCenterConesNum = theLconeDominanceOptimizationResultsTargetCenterConesNum;
    else

	    [theMconeDominanceOptimizationResultsFileNameCollection, ...
	     theMconeDominanceOptimizationResultsTargetCenterConesNum] = ...
			    RGCMosaicConstructor.compute.poolingFunctionsForSurroundOptimizationGrid(...
		 		    whichEye,  mosaicEccDegs, mosaicSizeDegs, spatialChromaticUniformityTradeoff , ...
		 		    customLMSconeDensities, surroundConnectivitySimulationParamsStruct, ...
		            'optimizationPositionIndicesToCompute', optimizationPositionIndicesToCompute, ...
		            'centerConeNumerositiesToOptimize', [], ...
		            'centerConeDominanceToOptimize', cMosaic.MCONE_ID, ...
		 		    'computeInputConeMosaicResponses', ~true, ...
				    'optimizeSurroundConePooling', ~true, ...
				    'onlyReturnSurroundOptimizationResultFilenames', true);
    end

	% Generate the sampling grid for all optimized pooling functions
	nodesNum = numel(theLconeDominanceOptimizationResultsFileNameCollection) + numel(theMconeDominanceOptimizationResultsFileNameCollection);
	theOptimizationGrid.positionDegs = zeros(nodesNum,2);
	theOptimizationGrid.positionColor = zeros(nodesNum,3);
	theOptimizationGrid.centerConeNumerosity = zeros(nodesNum,1);
	theOptimizationGrid.centerConeDominance = zeros(nodesNum,1);
	theOptimizationGrid.computeStruct = cell(nodesNum,1);

	gridNodeIndex = 0;
	for coneDominance = [cMosaic.LCONE_ID cMosaic.MCONE_ID]
		switch (coneDominance)
			case cMosaic.LCONE_ID
				theOptimizationResultsFileNameCollection = theLconeDominanceOptimizationResultsFileNameCollection;
				theOptimizationResultsTargetCenterConesNum = theLconeDominanceOptimizationResultsTargetCenterConesNum;
				fittedPositionColor = [1 0 0];
			case cMosaic.MCONE_ID
				theOptimizationResultsFileNameCollection = theMconeDominanceOptimizationResultsFileNameCollection;
				theOptimizationResultsTargetCenterConesNum = theMconeDominanceOptimizationResultsTargetCenterConesNum;
				fittedPositionColor = [0 1 0];
		end

		for iOptimizationResultsFileIndex = 1:numel(theOptimizationResultsFileNameCollection)
			theOptimizationResultsFileName = theOptimizationResultsFileNameCollection{iOptimizationResultsFileIndex};

			% Load optimization results
			try
				load(theOptimizationResultsFileName, ...
		   			'theMRGCMosaic', ...
			        'theTargetRGCindex', ...
			        'optimizationResults');
			catch 
			   	error('\nThe optimization results file \n\t %s \n was **NOT** found.\n', theOptimizationResultsFileName);
			end % try

			% consistency check #1: RF center cone numerosity
            consistentRFcenterNumerosity = theOptimizationResultsTargetCenterConesNum(iOptimizationResultsFileIndex) == theMRGCMosaic.exclusivelyConnectedInputConeIndicesNum(theTargetRGCindex);
			assert(consistentRFcenterNumerosity, ...
				    'Inconsistency in RF center cone numerosity: %d vs %d\n', theOptimizationResultsTargetCenterConesNum(iOptimizationResultsFileIndex),  theMRGCMosaic.exclusivelyConnectedInputConeIndicesNum(theTargetRGCindex));

			% consistency check #2: RF center cone dominance
			s = theMRGCMosaic.singleCellConnectivityStats(theTargetRGCindex, 'center', ...
					'minConeWeightIncluded', mRGCMosaic.sensitivityAtPointOfOverlap);

            if (employLconeDominanceOptimizationOnly) && (s.dominantConeType == cMosaic.LCONE_ID) && (coneDominance == cMosaic.MCONE_ID)
                fprintf(2, '>>>>> WARNING: Reassigning L-cone dominance to M-cone dominance, as employLconeDominanceOptimizationOnly is set to true\n');
                s.dominantConeType = coneDominance;
            end

			assert(s.dominantConeType == coneDominance, ...
					'Inconsistency in RF center cones dominance: %d vs %d\n', s.dominantConeType, coneDominance);

			% Assemble compute struct
			theComputeStruct = struct();
			theComputeStruct.positionDegs = theMRGCMosaic.rgcRFpositionsDegs(theTargetRGCindex,:);
			theComputeStruct.coneDominance = coneDominance;
			theComputeStruct.coneNumerosity = theOptimizationResultsTargetCenterConesNum(iOptimizationResultsFileIndex);
			theComputeStruct.optimizationResults = optimizationResults;

			gridNodeIndex = gridNodeIndex + 1;
			theOptimizationGrid.positionDegs(gridNodeIndex,:) = theComputeStruct.positionDegs;
			theOptimizationGrid.positionColor(gridNodeIndex,:) = fittedPositionColor;
			theOptimizationGrid.centerConeNumerosity(gridNodeIndex) = theComputeStruct.coneNumerosity;
			theOptimizationGrid.centerConeDominance(gridNodeIndex) = coneDominance;
			theOptimizationGrid.computeStruct{gridNodeIndex} = theComputeStruct;

		end % for iOptimizationResultsFileIndex
	end % for coneDominance
	clear 'theMRGCMosaic';

	if (surroundConnectivitySimulationParamsStruct.poolingOptimizationParamsStruct.employRFCenterOverlappingMosaic)
    	centerConnectivityStage = 'center connected with overlap';
    	surroundConnectivityStage = 'surround connected with center overlap';
	else
    	centerConnectivityStage = 'center connected';
    	surroundConnectivityStage = 'surround connected';
	end

	% Generate filename for the center-connected mosaic and for the surround-connected (compute-ready) mosaic
	% that is to be generated 
	centerConnectedParamsStruct.whichEye = whichEye;
	centerConnectedParamsStruct.eccentricityDegs = mosaicEccDegs;
	centerConnectedParamsStruct.sizeDegs = mosaicSizeDegs;
	centerConnectedParamsStruct.spatialChromaticUniformityTradeoff = spatialChromaticUniformityTradeoff;
    centerConnectedParamsStruct.customLMSconeDensities = customLMSconeDensities;

	% The center-connected mRGC mosaic filename (input)
	[theCenterConnectedMRGCMosaicFullFileName, ~, theCenterConnectedMRGCMosaicFileName] = ...
			    RGCMosaicConstructor.filepathFor.exportedMosaicFileName(...
			    centerConnectedParamsStruct, centerConnectivityStage);

    % Generate the surround-connected (compute-ready) mRGC mosaic filename (output)
    surroundConnectedParamsStruct = centerConnectedParamsStruct;

    % Additional info specific to the surround process
    surroundConnectedParamsStruct.surroundConnectivitySimulationParamsStruct = surroundConnectivitySimulationParamsStruct;
    surroundConnectedParamsStruct.surroundConnectivitySimulationParamsStruct.optimizationPosition = [];
    theSurroundConnectedMRGCMosaicFullFileName = RGCMosaicConstructor.filepathFor.exportedMosaicFileName(...
			    surroundConnectedParamsStruct, surroundConnectivityStage);

    hFig = [];
	if (visualizeOptimizationGridOnTopOfMosaic)
		minConeWeightVisualized = mRGCMosaic.sensitivityAtPointOfOverlap;
		pStruct = struct(...
            'maxNumberOfConesOutsideContour', 0, ...
            'spatialChromaticUniformityTradeoff', spatialChromaticUniformityTradeoff ...
		);
		
		% Visualize it along with the optimized pooling function sampling position grid
		[hFig, theAxes] = RGCMosaicConstructor.visualize.fullMosaic(theCenterConnectedMRGCMosaicFullFileName, theCenterConnectedMRGCMosaicFileName, ...
			minConeWeightVisualized, pStruct, ...
			'visualizedSamplingPositionsGrid', theOptimizationGrid.positionDegs, ...
			'visualizedSamplingPositionsGridColor', theOptimizationGrid.positionColor, ...
			'withFigureFormat', PublicationReadyPlotLib.figureComponents('1x1 giant rectangular-wide mosaic'));

		% Superipose lines showing the relationship to the desired sampling position
		hold(theAxes, 'on')

		% Identify the intended optimization grid position with a star
		plot(theAxes, surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(:,1), surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(:,2), ...
			'yh', 'MarkerSize', 20, 'MarkerFaceColor', [1 0.5 0], 'LineWidth', 1.5);

		theTargetConeNumerosities = unique(theOptimizationGrid.centerConeNumerosity);
		theTargetConeDominances = [cMosaic.LCONE_ID cMosaic.MCONE_ID];

		for iNumerosity = 1:numel(theTargetConeNumerosities)
			theTargetConeNumerosity = theTargetConeNumerosities(iNumerosity);
			for iDominance = 1:numel(theTargetConeDominances)
				theTargetConeDominance = theTargetConeDominances(iDominance);
				if (employLconeDominanceOptimizationOnly) && (theTargetConeDominance == cMosaic.MCONE_ID)
					continue;
				end
				if (theTargetConeDominance == cMosaic.LCONE_ID)
					theDominantConeColor = [1 0 0];
				else
					theDominantConeColor = [0 1 0];
                end

                if (visualizeNeighboringOptimizationGridNodesWithLines)
				    for idx = 1:numel(optimizationPositionIndicesToCompute)
					    iPos = optimizationPositionIndicesToCompute(idx);
					    kdx = find((theOptimizationGrid.centerConeNumerosity == theTargetConeNumerosity) & (theOptimizationGrid.centerConeDominance == theTargetConeDominance));
					    ddd = sum((bsxfun(@minus, theOptimizationGrid.positionDegs(kdx,:),surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(iPos,:))).^2,2);
					    [~,iidx] = sort(ddd);
    
					    % Plot at most the 3 closest optimization positions with the same numerisity and cone dominance
					    for jjj = 1:min([3 numel(iidx)])
						    theGridNodeIndex = kdx(iidx(jjj));
    
						    xx = [surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(iPos,1) theOptimizationGrid.positionDegs(theGridNodeIndex,1)];
						    yy = [surroundConnectivitySimulationParamsStruct.optimizationPositionsGrid(iPos,2) theOptimizationGrid.positionDegs(theGridNodeIndex,2)];
						    plot(theAxes, xx, yy, 'r-', 'LineWidth', 2.0, 'Color', theDominantConeColor);
						    text(theAxes, xx(2), yy(2), sprintf('%d', theTargetConeNumerosity), 'FontSize', 14);
					    end
    
				    end % idx
                end % if (visualizeNeighboringOptimizationGridNodesWithLines)

			end % for theTargetConeDominance
		end % for theTargetConeNumerosity

		theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    	thePDFfileName = fullfile(theRawFiguresDir, strrep(theCenterConnectedMRGCMosaicFileName, '.mat', 'WithOptimizationGridSuprimposed.pdf' ));
    	NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
	end % visualizeOptimizationGridOnTopOfMosaic

	if (onlyVisualizeOptimizationGrid)
		fprintf('Skipping computation of computeReadyMosaic after visualizing the sampling grid.\n');
		return;
	end

    % Load the center-connected mosaic
    load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');

	

	%interpolationPolicy = 'match cone numerosity';
	interpolationPolicy = 'match cone numerosity and cone dominance';
	adjustPoolingWeightsToMatchCenterStrengthsWithOptimizedModel = ~true;
	adjustPoolingWeightsToMatchSurroundStrengthsWithOptimizedModel = ~true;
	maxNearbyCandidates = 3;

	% Indices to generate the sparse surround connectivity matrix
	surroundConeIndicesAllRGCs = cell(theMRGCMosaic.rgcsNum,1);
    surroundConeWeightsAllRGCs = cell(theMRGCMosaic.rgcsNum,1);
    surroundRGCindicesAllRGCs = cell(theMRGCMosaic.rgcsNum,1);

    % Find the indices of all compute structs with the same center cone dominance as this RGC
    gridNodeIndicesWithLconeDominance = find(theOptimizationGrid.centerConeDominance == cMosaic.LCONE_ID);
    gridNodeIndicesWithMconeDominance = find(theOptimizationGrid.centerConeDominance == cMosaic.MCONE_ID);



    % Determine unique # of center cones num in the mosaic
	allRGCCenterConesNum = theMRGCMosaic.allRFcenterConnectivityStats();
	allConeNumerosities = unique(allRGCCenterConesNum);

    % Downsample the population of centerConesNum
    downsampledCenterConesNum = RGCMosaicConstructor.compute.downsampledRFcenterConeNumerosityPopulation(allConeNumerosities);


    if (visualizeInterpolationProcess)
    	for theCurrentRGCindex = 1:theMRGCMosaic.rgcsNum
			[theCurrentRGCsurroundConeIndices, theCurrentRGCsurroundConeWeights] = computeSurroundConePoolingWeights(...
				theMRGCMosaic, theCurrentRGCindex, downsampledCenterConesNum, ...
				surroundConnectivitySimulationParamsStruct.poolingOptimizationParamsStruct, ...
				theOptimizationGrid, interpolationPolicy, ...
				gridNodeIndicesWithLconeDominance, gridNodeIndicesWithMconeDominance, ...
				adjustPoolingWeightsToMatchCenterStrengthsWithOptimizedModel, ...
				adjustPoolingWeightsToMatchSurroundStrengthsWithOptimizedModel, ...
				surroundVarianceInComputeReadyMosaic, ...
				maxNearbyCandidates, ...
				visualizeInterpolationProcess);

			surroundConeIndicesAllRGCs{theCurrentRGCindex} = theCurrentRGCsurroundConeIndices;
	        surroundConeWeightsAllRGCs{theCurrentRGCindex} = theCurrentRGCsurroundConeWeights;
	        surroundRGCindicesAllRGCs{theCurrentRGCindex} = repmat(theCurrentRGCindex, [numel(theCurrentRGCsurroundConeIndices) 1]);
		end % for iRGC
	else
		parfor theCurrentRGCindex = 1:theMRGCMosaic.rgcsNum
			[theCurrentRGCsurroundConeIndices, theCurrentRGCsurroundConeWeights] = computeSurroundConePoolingWeights(...
                theMRGCMosaic, theCurrentRGCindex, downsampledCenterConesNum, ...
					surroundConnectivitySimulationParamsStruct.poolingOptimizationParamsStruct, ...
					theOptimizationGrid, interpolationPolicy, ...
					gridNodeIndicesWithLconeDominance, gridNodeIndicesWithMconeDominance, ...
					adjustPoolingWeightsToMatchCenterStrengthsWithOptimizedModel, ...
					adjustPoolingWeightsToMatchSurroundStrengthsWithOptimizedModel, ...
					surroundVarianceInComputeReadyMosaic, ...
					maxNearbyCandidates, ...
					visualizeInterpolationProcess);

			surroundConeIndicesAllRGCs{theCurrentRGCindex} = theCurrentRGCsurroundConeIndices;
	        surroundConeWeightsAllRGCs{theCurrentRGCindex} = theCurrentRGCsurroundConeWeights;
	        surroundRGCindicesAllRGCs{theCurrentRGCindex} = repmat(theCurrentRGCindex, [numel(theCurrentRGCsurroundConeIndices) 1]);
		end % parfor iRGC
	end


	% Generate the sparse matrix encoding cone connections of cones  the RF surrounds of all cells in the MRGCMosaic
    surroundConeIndices = vertcat(surroundConeIndicesAllRGCs{:});
    surroundConeWeights = vertcat(surroundConeWeightsAllRGCs{:});
    rgcIndices = vertcat(surroundRGCindicesAllRGCs{:});
    rgcRFsurroundConeConnectivityMatrix = sparse(...
    	surroundConeIndices, rgcIndices, surroundConeWeights, ...
    	theMRGCMosaic.inputConeMosaic.conesNum, theMRGCMosaic.rgcsNum);

    % Set the surround connectivity matrix and also save the employed surroundConnectivitySimulationParamsStruct 
    % and theOptimizationGrid
    surroundConnectivitySimulationParamsStruct.employedOptimizationGrid = theOptimizationGrid;
    theMRGCMosaic.bakeSurroundConeConnectivityMatrixAndFreeze(rgcRFsurroundConeConnectivityMatrix, surroundConnectivitySimulationParamsStruct, surroundVarianceInComputeReadyMosaic);

    % Save the compute-ready mosaic
    fprintf('Saving compute-ready mosaic in %s\n', theSurroundConnectedMRGCMosaicFullFileName);
    save(theSurroundConnectedMRGCMosaicFullFileName, 'theMRGCMosaic', '-v7.3');
    fprintf('Done !!\n');
end



function [theSurroundConeIndices, theSurroundConeWeights] = computeSurroundConePoolingWeights(...
		theRGCMosaic, theCurrentRGCindex, downsampledCenterConesNum, ...
		poolingOptimizationParamsStruct, theOptimizationGrid, interpolationPolicy, ...
		gridNodeIndicesWithLconeDominance, gridNodeIndicesWithMconeDominance, ...
		adjustPoolingWeightsToMatchCenterStrengthsWithOptimizedModel, ...
		adjustPoolingWeightsToMatchSurroundStrengthsWithOptimizedModel, ...
		surroundVarianceInComputeReadyMosaic, ...
		maxNearbyCandidates, ...
		visualizeInterpolationProcess)

	% The cell's position
    theCurrentRGCposition = theRGCMosaic.rgcRFpositionsDegs(theCurrentRGCindex,:);

    % Spectral interpolation coefficients based on the cell's center cone numerosity at the default minConeWeight
    s = theRGCMosaic.singleCellConnectivityStats(theCurrentRGCindex, 'center');
    % Actual numerosity
    actualConeNumerosity = s.inputConesNum;

    % Downsampled center cone numerosity
    [~,idx] = min(abs(downsampledCenterConesNum-actualConeNumerosity));
    downSampledConeNumerosity = downsampledCenterConesNum(idx);
    fprintf('Downsampled numerosity for RGC with %d center cones: %d\n', actualConeNumerosity, downSampledConeNumerosity);

    % Spectral interpolation coefficients based on the cell's center cone composition
    % Including cones with weights down to 1%
    s = theRGCMosaic.singleCellConnectivityStats(theCurrentRGCindex, 'center', ...
    	'minConeWeightIncluded', 0.01, ...
    	'warnIfCenterConeInputNumerosityDiffersFromExclusiveOne', false);
    theLconeComputeStructWeight = s.netWeights(cMosaic.LCONE_ID);
    theMconeComputeStructWeight = s.netWeights(cMosaic.MCONE_ID);

    % Find the indices of all compute structs with the same center cone numerosity as this RGC
    gridNodeIndicesMatchingCenterConeNumerosity = find(theOptimizationGrid.centerConeNumerosity == downSampledConeNumerosity);
    

    % Form optimized cone numerosities string
    optimizedConeNumerosities = unique(theOptimizationGrid.centerConeNumerosity);
    optimizedConeNumerositiesString = '';
    for i = 1:numel(optimizedConeNumerosities)
    	optimizedConeNumerositiesString = sprintf('%s %d ', optimizedConeNumerositiesString, optimizedConeNumerosities(i));
    end

    if (visualizeInterpolationProcess)
    	visualizationData = struct();
    	visualizationData.allGridNodesPositionsDegs = theOptimizationGrid.positionDegs;
    	visualizationData.allGridNodesConeDominance = theOptimizationGrid.centerConeDominance;
    	visualizationData.currentRGCposition = theCurrentRGCposition;
    	visualizationData.spectralWeights(cMosaic.LCONE_ID) = theLconeComputeStructWeight;
    	visualizationData.spectralWeights(cMosaic.MCONE_ID) = theMconeComputeStructWeight;
    end

    
    switch (interpolationPolicy)
	    case  'match cone numerosity and cone dominance'

            
	    	candidateNodeIndicesWithMatchingConeNumerosityAndLconeDominance = ...
                intersect(gridNodeIndicesMatchingCenterConeNumerosity, gridNodeIndicesWithLconeDominance);

	    	[gridNodeIndicesWithLconeDominance, gridNodeWeightsWithLconeDominance] = ...
                triangulatingGridNodeIndicesAndWeights(maxNearbyCandidates,...
	    		candidateNodeIndicesWithMatchingConeNumerosityAndLconeDominance, theOptimizationGrid.positionDegs, ...
	    		theCurrentRGCposition);
	    	

	    	candidateNodeIndicesWithMatchingConeNumerosityAndMconeDominance = ...
                intersect(gridNodeIndicesMatchingCenterConeNumerosity, gridNodeIndicesWithMconeDominance);

	    	[gridNodeIndicesWithMconeDominance, gridNodeWeightsWithMconeDominance] = ...
                triangulatingGridNodeIndicesAndWeights(maxNearbyCandidates,...
	    		candidateNodeIndicesWithMatchingConeNumerosityAndMconeDominance, theOptimizationGrid.positionDegs, ...
	    		theCurrentRGCposition);

	    	% Adjust weights based on the cells L/M weights
	    	gridNodeWeightsWithLconeDominance = gridNodeWeightsWithLconeDominance * theLconeComputeStructWeight;
	    	gridNodeWeightsWithMconeDominance = gridNodeWeightsWithMconeDominance * theMconeComputeStructWeight;


            if isempty(gridNodeIndicesWithLconeDominance)
                gridNodeWeightsWithMconeDominance = gridNodeWeightsWithMconeDominance / sum(gridNodeWeightsWithMconeDominance(:));
            end

            if isempty(gridNodeIndicesWithMconeDominance)
                gridNodeWeightsWithLconeDominance = gridNodeWeightsWithLconeDominance / sum(gridNodeWeightsWithLconeDominance(:));
            end
            
	    	gridNodeIndices = cat(1, gridNodeIndicesWithLconeDominance(:), gridNodeIndicesWithMconeDominance(:));
	    	gridNodeWeights = cat(1, gridNodeWeightsWithLconeDominance(:), gridNodeWeightsWithMconeDominance(:));

	    	if (visualizeInterpolationProcess)
		    	visualizationData.triangulatingOptimizedModelSpatialPositions(cMosaic.LCONE_ID,:,:) = ...
		    		theOptimizationGrid.positionDegs(gridNodeIndicesWithLconeDominance,:);
		    	visualizationData.triangulatingOptimizedModelSpatialPositions(cMosaic.MCONE_ID,:,:) = ...
		    		theOptimizationGrid.positionDegs(gridNodeIndicesWithMconeDominance,:);
		    	visualizationData.triangulatingOptimizedModelSpatialWeights(cMosaic.LCONE_ID, :) = ...
		    		reshape(gridNodeWeightsWithLconeDominance, [1 numel(gridNodeWeightsWithLconeDominance)]);
		    	visualizationData.triangulatingOptimizedModelSpatialWeights(cMosaic.MCONE_ID, :) = ...
		    		reshape(gridNodeWeightsWithMconeDominance(:), [1 numel(gridNodeWeightsWithMconeDominance)]);
   		 	end

	    case 'match cone numerosity'
	    	[gridNodeIndices, gridNodeWeights] = triangulatingGridNodeIndicesAndWeights(maxNearbyCandidates,...
	    		gridNodeIndicesMatchingCenterConeNumerosity, theOptimizationGrid.positionDegs, ...
	    		theCurrentRGCposition);

	    	if (visualizeInterpolationProcess)
		    	visualizationData.triangulatingOptimizedModelSpatialPositions = reshape(gridNodeIndices, [1 numel(gridNodeIndices)]);
		    	visualizationData.triangulatingOptimizedModelSpatialWeights = reshape(gridNodeWeights, [1 numel(gridNodeWeights)]);
   		 	end

	    otherwise
	    	error('Unknown intepolationPolity: ''%s''.', interpolationPolicy);
    end

    if (numel(gridNodeIndices) == 0)
    	fprintf('----> No optimized models for RGC with %d center cone numerosity. Available numerosities: %s\n', s.inputConesNum, optimizedConeNumerositiesString);
    	
    else
    	fprintf('Evaluating %d optimized models at the vicinity of RGC %d to derive its surround cone weights\n', ...
			numel(gridNodeIndices), theCurrentRGCindex);
    end

    
    if (abs(sum(gridNodeWeights)-1) > 10*eps)
    	error('sum(gridNodeWeights (%d) does not sum to 1.0: %f', ...
            numel(gridNodeIndices), sum(gridNodeWeights));
    end

    
    % Evaluate each compute struct at the current RGC to compute cone weights
    pooledConeIndicesAndWeights = cell(1, numel(gridNodeIndices));
    randomVectorForSurroundVariance = randn(1,1);

    for idx = 1:numel(gridNodeIndices)

    	theGridNodeIndex = gridNodeIndices(idx);
    	theComputeStruct = theOptimizationGrid.computeStruct{theGridNodeIndex};
    	
    	optimizedModelTotalCenterStrength = sum(theComputeStruct.optimizationResults.theFinalPooledConeIndicesAndWeights.centerConeWeights);
    	optimizedModelTotalSurroundStrength = sum(theComputeStruct.optimizationResults.theFinalPooledConeIndicesAndWeights.surroundConeWeights);

    	% Update the modelConstants for theCurrentRGCindex
    	modelConstantsForCurrentRGC = updateModelConstantsForCurrentRGC(...
    		theComputeStruct.optimizationResults.modelConstants, ...
    		theRGCMosaic, theCurrentRGCindex, poolingOptimizationParamsStruct);

    	% Update the modelConstants with the passed surroundVarianceInComputeReadyMosaic
    	if (~isempty(fieldnames(surroundVarianceInComputeReadyMosaic)))
    		surroundVarianceInComputeReadyMosaic.randomVector = randomVectorForSurroundVariance;
    		modelConstantsForCurrentRGC.surroundVarianceInComputeReadyMosaic = surroundVarianceInComputeReadyMosaic;
    	end
    	
    	% Evaludate the optimized model{idx} for the modelConstantsForCurrentRGC
     	pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC = ...
     		theComputeStruct.optimizationResults.modelConstants.poolingModel.weightsComputeHandle(...
         		modelConstantsForCurrentRGC, theComputeStruct.optimizationResults.modelVariables.finalValues);
 
 		if (adjustPoolingWeightsToMatchCenterStrengthsWithOptimizedModel)
	 		% Find factor that makes the net center sensitivity at the current RGC equal to that of the optimized model
	 		% and multiply the surround weights by that factor
	 		factorToMatchCenterNetSensitivity = optimizedModelTotalCenterStrength/sum(pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC.centerConeWeights);
	 		pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC.surroundConeWeights = ...
	 			pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC.surroundConeWeights * factorToMatchCenterNetSensitivity;
		end

		if (adjustPoolingWeightsToMatchSurroundStrengthsWithOptimizedModel)
			factorToMatchSurroundNetSensitivity = optimizedModelTotalSurroundStrength/sum(pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC.surroundConeWeights);
			pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC.surroundConeWeights = ...
	 			pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC.surroundConeWeights * factorToMatchSurroundNetSensitivity;
		end

		% We do not need the center cone weights anymore
		pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC = ...
			rmfield(pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC, 'centerConeWeights');
		pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC = ...
			rmfield(pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC, 'centerConeIndices');

     	pooledConeIndicesAndWeights{idx} = pooledConeIndicesAndWeightsOfModelEvaluatedAtCurrentRGC;
    end % idx


    surroundConeIndices = [];
    surroundConeWeights = [];
    % Accumulate all the surround cones
    for idx = 1:numel(gridNodeIndices)
    	% Surround cones, multiplied by the weight of the grid node
        theModelSurroundConeWeights = pooledConeIndicesAndWeights{idx}.surroundConeWeights * gridNodeWeights(idx);
        theModelSurroundConeIndices = pooledConeIndicesAndWeights{idx}.surroundConeIndices;

        % Accumulate all surround pooling weights from the different grid nodes
        [surroundConeIndices, surroundConeWeights] = accumulateSurroundWeights(...
                surroundConeIndices, surroundConeWeights, ...
                theModelSurroundConeIndices, theModelSurroundConeWeights);

        if (visualizeInterpolationProcess)
        	switch (interpolationPolicy)
	    		case  'match cone numerosity and cone dominance'
	    			if (idx  <= numel(gridNodeWeightsWithLconeDominance))
	    				visualizationData.triangulatingOptimizedModelSurroundConeWeights{cMosaic.LCONE_ID, idx} = theModelSurroundConeWeights/gridNodeWeights(idx);
	    				visualizationData.triangulatingOptimizedModelSurroundConeIndices{cMosaic.LCONE_ID, idx} = theModelSurroundConeIndices;
	    			else
	    				visualizationData.triangulatingOptimizedModelSurroundConeWeights{cMosaic.MCONE_ID, idx-numel(gridNodeWeightsWithLconeDominance)} = theModelSurroundConeWeights;
	    				visualizationData.triangulatingOptimizedModelSurroundConeIndices{cMosaic.MCONE_ID, idx-numel(gridNodeWeightsWithLconeDominance)} = theModelSurroundConeIndices;
	    			end

	    		case 'match cone numerosity'
		   			visualizationData.triangulatingOptimizedModelSurroundConeWeights{idx} = theModelSurroundConeWeights/gridNodeWeights(idx);
		    		visualizationData.triangulatingOptimizedModelSurroundConeIndices{idx} = theModelSurroundConeIndices;
		    end
   		end
    end

    % Assemble struct with surround cone indices and weights
    theSurroundConeIndices = surroundConeIndices(:);
    theSurroundConeWeights = surroundConeWeights(:);

    if (visualizeInterpolationProcess)
    	thresholdSurroundConeWeightForInclusionInVisualization = 0.5/100;
    	surroundProfile = 'spatialMap2D';
    	surroundProfile = 'horizontalLineWeightingFunction';
    	surroundProfile = 'combo';
    	visualizeFullOptimizatioGrid = true;

    	visualizationData.currentRGCinterpolatedSurroundConeWeights = theSurroundConeWeights;
    	visualizationData.currentRGCinterpolatedSurroundConeIndices = theSurroundConeIndices;

    	hFig = visualizeInterpolationMethod(theRGCMosaic, visualizationData, ...
    		thresholdSurroundConeWeightForInclusionInVisualization, surroundProfile, visualizeFullOptimizatioGrid);
    end
end

function hFigSummary = visualizeInterpolationMethod(theRGCMosaic, visualizationData, ...
	thresholdSurroundConeWeightForInclusionInVisualization, surroundProfile, visualizeFullOptimizatioGrid)

	hFigSummary = figure(2000); clf;
    set(hFigSummary, 'Name', 'Summary', 'Position', [10 10 2000 1000], 'Color', [1 1 1]);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 4, ...
           'heightMargin',  0.16, ...
           'widthMargin',    0.06, ...
           'leftMargin',     0.05, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.10, ...
           'topMargin',      0.02);
    
    axInterpolatedSamples = subplot('Position', subplotPosVectors(1,1).v);
    axInterpolatedSurroundWeights = subplot('Position', subplotPosVectors(2,1).v);
    axTriangulatingOptimizedModelSurroundWeights{cMosaic.LCONE_ID,1} = subplot('Position', subplotPosVectors(1,2).v);
    axTriangulatingOptimizedModelSurroundWeights{cMosaic.LCONE_ID,2} = subplot('Position', subplotPosVectors(1,3).v);
    axTriangulatingOptimizedModelSurroundWeights{cMosaic.LCONE_ID,3} = subplot('Position', subplotPosVectors(1,4).v);
    axTriangulatingOptimizedModelSurroundWeights{cMosaic.MCONE_ID,1} = subplot('Position', subplotPosVectors(2,2).v);
    axTriangulatingOptimizedModelSurroundWeights{cMosaic.MCONE_ID,2} = subplot('Position', subplotPosVectors(2,3).v);
    axTriangulatingOptimizedModelSurroundWeights{cMosaic.MCONE_ID,3} = subplot('Position', subplotPosVectors(2,4).v);

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    hold(axInterpolatedSamples, 'on');

    % Lines depicting spatial interpolation weights
    for iTriangulatingSample = 1:size(visualizationData.triangulatingOptimizedModelSpatialPositions,2)
    	theSpatioSpectralWeight = visualizationData.triangulatingOptimizedModelSpatialWeights(cMosaic.LCONE_ID)*visualizationData.triangulatingOptimizedModelSpatialWeights(cMosaic.LCONE_ID,iTriangulatingSample);
    	if (theSpatioSpectralWeight>0)
	    	plot(axInterpolatedSamples, ...
	    		[visualizationData.currentRGCposition(1) visualizationData.triangulatingOptimizedModelSpatialPositions(cMosaic.LCONE_ID, iTriangulatingSample, 1)], ...
	    		[visualizationData.currentRGCposition(2) visualizationData.triangulatingOptimizedModelSpatialPositions(cMosaic.LCONE_ID, iTriangulatingSample, 2)], ...
	    		'r-', 'LineWidth', 10*ff.lineWidth*theSpatioSpectralWeight);
	    end
    end

    for iTriangulatingSample = 1:size(visualizationData.triangulatingOptimizedModelSpatialPositions,2)
    	theSpatioSpectralWeight = visualizationData.triangulatingOptimizedModelSpatialWeights(cMosaic.MCONE_ID) * ...
    	                          visualizationData.triangulatingOptimizedModelSpatialWeights(cMosaic.MCONE_ID,iTriangulatingSample);
    	if (theSpatioSpectralWeight>0)
	    	plot(axInterpolatedSamples, ...
	    		[visualizationData.currentRGCposition(1) visualizationData.triangulatingOptimizedModelSpatialPositions(cMosaic.MCONE_ID, iTriangulatingSample, 1)], ...
	    		[visualizationData.currentRGCposition(2) visualizationData.triangulatingOptimizedModelSpatialPositions(cMosaic.MCONE_ID, iTriangulatingSample, 2)], ...
	    		'g-', 'Color', [0 0.75 0.0], 'LineWidth', 10*ff.lineWidth*theSpatioSpectralWeight);
    	end
    end

    if (visualizeFullOptimizatioGrid)
	    idx = find(visualizationData.allGridNodesConeDominance == cMosaic.LCONE_ID);
	    lConeDominanceNodePositions = visualizationData.allGridNodesPositionsDegs(idx,:);
	    scatter(axInterpolatedSamples, lConeDominanceNodePositions(:,1), lConeDominanceNodePositions(:,2), ...
	    	(ff.markerSize-4)^2, 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.7 0 0], 'LineWidth', ff.lineWidth);

	    idx = find(visualizationData.allGridNodesConeDominance == cMosaic.MCONE_ID);
	    mConeDominanceNodePositions = visualizationData.allGridNodesPositionsDegs(idx,:);
	    scatter(axInterpolatedSamples, mConeDominanceNodePositions(:,1), mConeDominanceNodePositions(:,2), ...
	    	(ff.markerSize-4)^2, 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', [0.5 0.9 0.5], 'MarkerEdgeColor', [0 0.8 0], 'LineWidth', ff.lineWidth);
	else
	    scatter(axInterpolatedSamples, ...
	    	visualizationData.triangulatingOptimizedModelSpatialPositions(cMosaic.LCONE_ID, :, 1), ...
	    	visualizationData.triangulatingOptimizedModelSpatialPositions(cMosaic.LCONE_ID, :, 2), ...
	    	ff.markerSize^2, 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.7 0 0], 'LineWidth', ff.lineWidth);

	    scatter(axInterpolatedSamples, ...
	    	visualizationData.triangulatingOptimizedModelSpatialPositions(cMosaic.MCONE_ID, :, 1), ...
	    	visualizationData.triangulatingOptimizedModelSpatialPositions(cMosaic.MCONE_ID, :, 2), ...
	    	ff.markerSize^2, 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', [0.5 0.9 0.5], 'MarkerEdgeColor', [0 0.6 0],  'LineWidth', ff.lineWidth);
	end

    plot(axInterpolatedSamples, visualizationData.currentRGCposition(1), visualizationData.currentRGCposition(2), ...
    	'kh', 'MarkerFaceColor', [0.8 0.8 1], 'MarkerEdgeColor', [0.2 0.2 0.3], 'MarkerSize', ff.markerSize+4, 'LineWidth', ff.lineWidth);

    % Set the limits
    if (~visualizeFullOptimizatioGrid)
	    theRange = max(visualizationData.triangulatingOptimizedModelSpatialPositions(:)) - min(visualizationData.triangulatingOptimizedModelSpatialPositions(:));
	    XLims = visualizationData.currentRGCposition(1) + theRange * 0.51 *[-1 1];
	    YLims = visualizationData.currentRGCposition(2) + theRange * 0.51 *[-1 1];
	else
		XLims(1) = min(squeeze(visualizationData.allGridNodesPositionsDegs(:,1)));
		XLims(2) = max(squeeze(visualizationData.allGridNodesPositionsDegs(:,1)));
		deltaX = 0.02*(XLims(2)-XLims(1));
		XLims(1) = XLims(1)-deltaX;
		XLims(2) = XLims(2)+deltaX;
		YLims(1) = min(squeeze(visualizationData.allGridNodesPositionsDegs(:,2)));
		YLims(2) = max(squeeze(visualizationData.allGridNodesPositionsDegs(:,2)));
		deltaY = 0.02*(YLims(2)-YLims(1));
		YLims(1) = YLims(1)-deltaY;
		YLims(2) = YLims(2)+deltaY;
	end

    set(axInterpolatedSamples, 'XLim', XLims, 'YLim', YLims);
    
    set(axInterpolatedSamples, 'XTick', -(3):0.3:(3), 'YTick', -3:0.3:3);
    xlabel(axInterpolatedSamples, 'eccentricity, x (degs)');
	ylabel(axInterpolatedSamples, 'eccentricity, y (degs)');

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(axInterpolatedSamples, ff);


    maxSurroundConeWeight = 0;
    maxXYrange = 0;
    for coneDominance = [cMosaic.LCONE_ID cMosaic.MCONE_ID]
	    for iTriangulatingSample = 1:size(visualizationData.triangulatingOptimizedModelSpatialPositions,2)
	    	theWeights = visualizationData.triangulatingOptimizedModelSurroundConeWeights{coneDominance, iTriangulatingSample};
	    	theSurroundConeIndices = visualizationData.triangulatingOptimizedModelSurroundConeIndices{coneDominance, iTriangulatingSample};
			theSurroundConePositions = theRGCMosaic.inputConeMosaic.coneRFpositionsDegs(theSurroundConeIndices,:);
	    	if (max(theWeights)>maxSurroundConeWeight)
	    		maxSurroundConeWeight = max(theWeights);
	    	end
	    	xRange = max(theSurroundConePositions(:,1))-min(theSurroundConePositions(:,1));
			yRange = max(theSurroundConePositions(:,2))-min(theSurroundConePositions(:,2));
			xyRange = max([xRange yRange]);
			if (xyRange > maxXYrange)
				maxXYrange = xyRange;
			end
	    end
	end


	% The triangulated optimized model weights
    for coneDominance = [cMosaic.LCONE_ID cMosaic.MCONE_ID]
    	if (visualizationData.spectralWeights(coneDominance) > 0)
		    for iTriangulatingSample = 1:size(visualizationData.triangulatingOptimizedModelSpatialPositions,2)
				theSurroundConeWeights = visualizationData.triangulatingOptimizedModelSurroundConeWeights{coneDominance, iTriangulatingSample};
				theSurroundConeIndices = visualizationData.triangulatingOptimizedModelSurroundConeIndices{coneDominance, iTriangulatingSample};
				renderSurroundConePoolingWeightsMap(axTriangulatingOptimizedModelSurroundWeights{coneDominance, iTriangulatingSample}, ...
					coneDominance, visualizationData.triangulatingOptimizedModelSurroundConeWeights{coneDominance, iTriangulatingSample}, ...
					visualizationData.triangulatingOptimizedModelSurroundConeIndices{coneDominance, iTriangulatingSample}, ...
					visualizationData.spectralWeights(coneDominance), ...
					visualizationData.triangulatingOptimizedModelSpatialWeights(coneDominance,iTriangulatingSample), ...
					theRGCMosaic, surroundProfile, ...
					thresholdSurroundConeWeightForInclusionInVisualization, ...
					 maxXYrange);
				PublicationReadyPlotLib.applyFormat(axTriangulatingOptimizedModelSurroundWeights{coneDominance, iTriangulatingSample}, ff);
			end % for iTriangulatingSample
		end % if (visualizationData.spectralWeights(coneDominance) > 0)
	end % for coneDominance

	% The interpolated weights
	renderSurroundConePoolingWeightsMap(axInterpolatedSurroundWeights, ...
		[], visualizationData.currentRGCinterpolatedSurroundConeWeights, ...
		visualizationData.currentRGCinterpolatedSurroundConeIndices, ...
		[], [], theRGCMosaic, surroundProfile, ...
		thresholdSurroundConeWeightForInclusionInVisualization, ...
		maxXYrange);

	theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
	pdfExportSubDir = '';
    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'surroundInterpolationDemo.pdf');
    NicePlot.exportFigToPDF(thePDFfileName, hFigSummary,  300);
end


function renderSurroundConePoolingWeightsMap(ax, coneDominance, theSurroundConeWeights, theSurroundConeIndices, ...
	spectralWeight,  spatialWeight, theRGCMosaic, surroundProfile, ...
	thresholdSurroundConeWeightForInclusionInVisualization, ...
	maxXYrange)

	theSurroundConePositions = theRGCMosaic.inputConeMosaic.coneRFpositionsDegs(theSurroundConeIndices,:);
	theSurroundConeTypes = theRGCMosaic.inputConeMosaic.coneTypes(theSurroundConeIndices);

	c = median(theSurroundConePositions,1);
	XLims = c(1) + maxXYrange*0.51*[-1 1];
	YLims = c(2) + maxXYrange*0.51*[-1 1];

	nSpatialSupport = 256;
	xSupport = linspace(XLims(1), XLims(2), nSpatialSupport);
	ySupport = linspace(YLims(1), YLims(2), nSpatialSupport);
	spatialSupportXYDegs = [xSupport(:) ySupport(:)];
	flatTopSaturationLevel = 0.4;

	[subregionConeMap, subregionConeMapFlatTop] = mRGCMosaic.retinalSubregionConeMapFromPooledConeInputs(...
      theRGCMosaic.inputConeMosaic, theSurroundConeIndices, theSurroundConeWeights, spatialSupportXYDegs, ...
      flatTopSaturationLevel);

	subregionConeMap = subregionConeMap / max(subregionConeMap(:));
	subregionConeMapFlatTop = subregionConeMapFlatTop / max(subregionConeMapFlatTop(:));


	idx = find(theSurroundConeWeights >= thresholdSurroundConeWeightForInclusionInVisualization);
	theSurroundConeWeights = theSurroundConeWeights(idx);
	theSurroundConeIndices = theSurroundConeIndices(idx);
	theSurroundConePositions = theSurroundConePositions(idx,:);
	theSurroundConeTypes = theSurroundConeTypes(idx);
	
	
	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hold(ax, 'on');
	if (strcmp(surroundProfile, 'spatialMap2D')) || (strcmp(surroundProfile, 'combo'))
		imagesc(ax,xSupport, ySupport, subregionConeMapFlatTop);

		for iConeIdx = 1:numel(theSurroundConeIndices)
			switch (theSurroundConeTypes(iConeIdx))
	    		case cMosaic.LCONE_ID
	    			theConeColor = [1 0 0];
	    		case cMosaic.MCONE_ID
	    			theConeColor = [0 0.8 0];
	    		case cMosaic.SCONE_ID
	    			theConeColor = [0 0 1];
			end
			w = 0;
			plot(ax, theSurroundConePositions(iConeIdx,1), theSurroundConePositions(iConeIdx,2), '.', ...
					'MarkerSize', ff.markerSize, ...
					'MarkerEdgeColor', theConeColor, ...
					'MarkerFaceColor', theConeColor, ...
					'LineWidth', 2.0);
		end % for iConeIdx
	end
	if (strcmp(surroundProfile,'horizontalLineWeightingFunction')) || (strcmp(surroundProfile, 'combo'))
		lineWeightingProfile = sum(subregionConeMap, 1)/nSpatialSupport;
		boostFactor = 30;
		lineWeightingProfile = YLims(1) + boostFactor*(YLims(2)-YLims(1)) * lineWeightingProfile;
		plot(ax, xSupport, lineWeightingProfile, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 2*ff.lineWidth);
		plot(ax, xSupport, lineWeightingProfile, '-', 'Color', [0.5 0.5 0.8], 'LineWidth', ff.lineWidth);
	end

	% Set the limits
	axis(ax, 'image'); axis(ax, 'xy'); 
	set(ax, 'XLim', XLims, 'YLim', YLims, 'CLim', [0 1]);
	set(ax, 'XTick', -5:0.1:5, 'YTick', -5:0.1:5);
	xlabel(ax, 'eccentricity, x (degs)');
	ylabel(ax, 'eccentricity, y (degs)');
	colormap(ax, brewermap(512, 'greys'));

	if (~isempty(spectralWeight)) && (~isempty(spatialWeight)) && (~isempty(coneDominance))
		title(ax, sprintf('interpolation weight: %2.2f', spectralWeight * spatialWeight));
	end

	PublicationReadyPlotLib.applyFormat(ax, ff);
end

function [surroundConeIndices, surroundConeWeights] = accumulateSurroundWeights(...
                surroundConeIndices, surroundConeWeights, ...
                newSurroundConeIndices, newSurroundConeWeights)
    % Find which of the newSurroundConeIndices already exist in the surroundConeIndices
    [ia,ib] = ismember(newSurroundConeIndices, surroundConeIndices);

    for i = 1:numel(ia)
        if (ia(i) == 0)
            % Include weights for previoulsy NOT included cone indices
            surroundConeWeights = cat(1, surroundConeWeights, newSurroundConeWeights(i));
            surroundConeIndices = cat(1, surroundConeIndices, newSurroundConeIndices(i));
        else
            % Accumulate weights for previously included cone indices
            surroundConeWeights(ib(i)) = surroundConeWeights(ib(i)) + newSurroundConeWeights(i); 
        end
    end
end


function updatedModelConstants = updateModelConstantsForCurrentRGC(modelConstants, ...
	theMRGCMosaic, theTargetRGCindex, poolingOptimizationParamsStruct)

	cachedDataForCurrentRGC = RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.cache(...
		theMRGCMosaic, theTargetRGCindex, poolingOptimizationParamsStruct);

	% Update cache for current RGC
	updatedModelConstants = modelConstants;
	updatedModelConstants.cachedData = cachedDataForCurrentRGC;
end


function [gridNodeIndices, gridNodeWeights] = triangulatingGridNodeIndicesAndWeights(...
	maxCandidatesNum, candidateGridNodeIndices, gridPositionDegs, theCurrentRGCposition)

	% Compute distances to all candidate grid node indices
	distancesToAllCandidateOptimizationPositions = sqrt(sum((bsxfun(@minus, gridPositionDegs(candidateGridNodeIndices,:), theCurrentRGCposition)).^2,2));
	
	% Find the maxCandidates closest grid node indices
	[~, idx] = sort(distancesToAllCandidateOptimizationPositions, 'ascend');
	if (numel(idx) > maxCandidatesNum)
    	idx = idx(1:maxCandidatesNum);
    end

    gridNodeIndices = candidateGridNodeIndices(idx);
	distancesToNearbyOptimizationPositions = distancesToAllCandidateOptimizationPositions(idx);

	% Avoid division by 0
	dd = 0.10*max(distancesToNearbyOptimizationPositions(:));

	% Node weights proportional to inverse distances raised to the power of 4
	pwr = -4.0;
	gridNodeWeights  = (distancesToNearbyOptimizationPositions + dd).^pwr;

	% Sum to 1.0
	gridNodeWeights = gridNodeWeights /sum(gridNodeWeights(:));
end
