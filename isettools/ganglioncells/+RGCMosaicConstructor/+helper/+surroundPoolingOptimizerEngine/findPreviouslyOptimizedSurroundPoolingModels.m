function [optimizationResults, theMRGCMosaic, theTargetRGCindex, targetVisualSTFparams, optimizationResultsFileName] = ...
		findPreviouslyOptimizedSurroundPoolingModels(whichEye,  mosaicEccDegs, mosaicSizeDegs, ...
			chromaticSpatialVarianceTradeoff,  customLMSconeDensities, rfSurroundConnectivityParams, varargin)
	% Parse optional input
    p = inputParser;
    p.addParameter('targetRFcenterConesNum', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('targetRFcenterDominantConeType', cMosaic.LCONE_ID, @(x)(ismember(x, [cMosaic.LCONE_ID cMosaic.MCONE_ID])));
    p.addParameter('initialSurroundOptimizationValuesSource', '', @(x)(isempty(x)||ischar(x)));
    p.addParameter('optimizationPositionIndicesToCompute', [], @isnumeric);
    
    p.parse(varargin{:});
    theRFcenterConesNum = p.Results.targetRFcenterConesNum;
    theRFcenterDominantConeType = p.Results.targetRFcenterDominantConeType;
    initialSurroundOptimizationValuesSource = p.Results.initialSurroundOptimizationValuesSource;
    optimizationPositionIndicesToCompute = p.Results.optimizationPositionIndicesToCompute;

	% Retrieve mosaic filenames that may have been optimized
	[~,~,theSurroundConnectedMRGCMosaicFullFileNames] = ...
		RGCMosaicConstructor.compute.poolingFunctionsForSurroundOptimizationGrid(whichEye,  mosaicEccDegs, mosaicSizeDegs, ...
			chromaticSpatialVarianceTradeoff, customLMSconeDensities, rfSurroundConnectivityParams, ...
			'optimizationPositionIndicesToCompute', optimizationPositionIndicesToCompute(1), ...
			'computeInputConeMosaicResponses', ~true, ...
			'optimizeSurroundConePooling', ~true, ...
			'onlyReturnSurroundOptimizationResultFilenames', true);

    optimizationResults = [];
    for i = 1:numel(theSurroundConnectedMRGCMosaicFullFileNames)
        if (~isempty(theSurroundConnectedMRGCMosaicFullFileNames{i}))
            theChosenSurroundConnectedMRGCMosaicFullFileName = theSurroundConnectedMRGCMosaicFullFileNames{i};
            % Attempt to import the previously-optimized file
	        [optimizationResults, theMRGCMosaic, theTargetRGCindex, targetVisualSTFparams, optimizationResultsFileName, matchedAttributes] = ...
		        RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.importOptimizationResults(...
			        theChosenSurroundConnectedMRGCMosaicFullFileName, theRFcenterConesNum, ...
			        theRFcenterDominantConeType, initialSurroundOptimizationValuesSource);
        end
    end

	

	if (isempty(optimizationResults))
		switch (theRFcenterDominantConeType)
			case cMosaic.LCONE_ID
				fprintf('\nThere are no previously optimized surround pooling models for %d-cone RF centers with L-cone dominance.\n', ...
					theRFcenterConesNum);
			case cMosaic.MCONE_ID 
				fprintf('\nThere are no previously optimized surround pooling models for %d-cone RF centers with M-cone dominance.\n', ...
					theRFcenterConesNum);
			case cMosaic.SCONE_ID
				fprintf('\nThere are no previously optimized surround pooling models for %d-cone RF centers with S-cone dominance.\n', ...
					theRFcenterConesNum);
		end % switch
	else
		if ((matchedAttributes.theRFcenterDominantConeType==theRFcenterDominantConeType)) && ((matchedAttributes.theRFcenterConesNum==theRFcenterConesNum))
			fprintf('\nFound previously optimized surround pooling model for %d-cone RF centers and desired cone dominance (%d).\n\n', theRFcenterConesNum, theRFcenterDominantConeType);
		else
			if (matchedAttributes.theRFcenterDominantConeType==theRFcenterDominantConeType)
				fprintf('\nFound previously optimized surround pooling model for %d-cone RF centers (instead of %d) with desired cone dominance (%d).\n\n', matchedAttributes.theRFcenterConesNum, theRFcenterConesNum, theRFcenterDominantConeType);
			end
			if (matchedAttributes.theRFcenterConesNum==theRFcenterConesNum)
				fprintf('\nFound previously optimized surround pooling model for %d-cone RF centers with different cone dominance (%d instead of %d).\n\n',  theRFcenterConesNum, matchedAttributes.theRFcenterDominantConeType, theRFcenterDominantConeType);
			end
		end
	end

end

