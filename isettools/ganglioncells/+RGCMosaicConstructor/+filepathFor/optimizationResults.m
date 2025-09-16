function  fName = optimizationResults(theSurroundConnectedMRGCMosaicFullFileName, ...
    	theTargetRGCRFcenterConesNum, theTargetRGCRFcenterDominantConeType, varargin)

    p = inputParser;
    p.addParameter('extraSubDirPath', '', @ischar);
    p.addParameter('generateMissingSubDirs', false, @islogical);
    % Execute the parser
    p.parse(varargin{:});
    extraSubDirPath = p.Results.extraSubDirPath;
    generateMissingSubDirs = p.Results.generateMissingSubDirs;

    % Root directory
    intermediateDataDir = RGCMosaicConstructor.filepathFor.intermediateDataDir();
    if (~isempty(extraSubDirPath))
        intermediateDataDir = fullfile(intermediateDataDir, extraSubDirPath);
    end

    theOptimizationResultsFullFileName = strrep(theSurroundConnectedMRGCMosaicFullFileName, 'surroundConnected', 'optResults');

    % OLD
	% if (contains(theSurroundConnectedMRGCMosaicFullFileName,'surroundConnectedOverlap'))
	% 	theOptimizationResultsFullFileName = strrep(theSurroundConnectedMRGCMosaicFullFileName, 'surroundConnectedOverlap', '');
	% 	theOptimizationResultsFullFileName = strrep(theOptimizationResultsFullFileName, 'SLIM', 'SLIM/optResultsOverlap');
	% else
	% 	theOptimizationResultsFullFileName = strrep(theSurroundConnectedMRGCMosaicFullFileName, 'surroundConnected', '');
	% 	theOptimizationResultsFullFileName = strrep(theOptimizationResultsFullFileName, 'SLIM', 'SLIM/optResults');
	% end

	switch (theTargetRGCRFcenterDominantConeType)
		case cMosaic.LCONE_ID
			targetInfoString = sprintf('_LconeDominated_%dconeRF.mat', theTargetRGCRFcenterConesNum);
		case cMosaic.MCONE_ID
			targetInfoString = sprintf('_MconeDominated_%dconeRF.mat', theTargetRGCRFcenterConesNum);
		case cMosaic.SCONE_ID
			targetInfoString = sprintf('_SconeDominated_%dconeRF.mat', theTargetRGCRFcenterConesNum);
	end
	fName = strrep(theOptimizationResultsFullFileName, '.mat', targetInfoString)

    fName = strrep(fName, intermediateDataDir , '')
    generateMissingSubDirs
    fName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
        intermediateDataDir , fName, ...
        'generateMissingSubDirs', generateMissingSubDirs)

end