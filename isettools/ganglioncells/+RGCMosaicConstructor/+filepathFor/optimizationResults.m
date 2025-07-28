function  fName = optimizationResults(theSurroundConnectedMRGCMosaicFullFileName, ...
    	theTargetRGCRFcenterConesNum, theTargetRGCRFcenterDominantConeType)

	if (contains(theSurroundConnectedMRGCMosaicFullFileName,'surroundConnectedOverlap'))
		theOptimizationResultsFullFileName = strrep(theSurroundConnectedMRGCMosaicFullFileName, 'surroundConnectedOverlap', '');
		theOptimizationResultsFullFileName = strrep(theOptimizationResultsFullFileName, 'SLIM', 'SLIM/optResultsOverlap');
	else
		theOptimizationResultsFullFileName = strrep(theSurroundConnectedMRGCMosaicFullFileName, 'surroundConnected', '');
		theOptimizationResultsFullFileName = strrep(theOptimizationResultsFullFileName, 'SLIM', 'SLIM/optResults');
	end

	switch (theTargetRGCRFcenterDominantConeType)
		case cMosaic.LCONE_ID
			targetInfoString = sprintf('_LconeDominated_%dconeRF.mat', theTargetRGCRFcenterConesNum);
		case cMosaic.MCONE_ID
			targetInfoString = sprintf('_MconeDominated_%dconeRF.mat', theTargetRGCRFcenterConesNum);
		case cMosaic.SCONE_ID
			targetInfoString = sprintf('_SconeDominated_%dconeRF.mat', theTargetRGCRFcenterConesNum);
	end
	fName = strrep(theOptimizationResultsFullFileName, '.mat', targetInfoString);
end