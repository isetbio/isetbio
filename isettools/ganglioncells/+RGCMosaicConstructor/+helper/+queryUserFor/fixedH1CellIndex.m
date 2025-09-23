function [theFixedH1CellIndex, theRFcenterConesNum, theRFcenterDominantConeType, theSelectedOptimizationResultsFileName] = fixedH1CellIndex(performIspectOptimizationResultsUsingGUI, employRFCenterOverlappingMosaic)

	theFixedH1CellIndex = [];
	theRFcenterConesNum = [];
	theRFcenterDominantConeType = [];
	
	theSelectedOptimizationResultsFileName = '';

	if (performIspectOptimizationResultsUsingGUI)
		if (employRFCenterOverlappingMosaic)
			subDir = 'optResultsOverlap';
		else
			subDir = 'optResults';
		end

		% Ask the user to select a specific filename using the finder gui
		[~, theSelectedOptimizationResultsFileName] = RGCMosaicConstructor.gui.loadFileData(...
			'subDir',subDir, ...
			'onlyReturnSelectedFileName', true);

		% Attempt to extract fixedH1CellIndex, theRFcenterConesNum and theRFcenterDominantConeType from theSelectedOptimizationResultsFileName
		if (~isempty(theSelectedOptimizationResultsFileName))
			targetSubString = 'PackerDacey2002H1FixedCellIndex';
		 	idx = strfind(theSelectedOptimizationResultsFileName, targetSubString);
		 	if (~isempty(idx))
		 		theFixedH1CellIndex  = str2num(theSelectedOptimizationResultsFileName(idx+numel(targetSubString)))
		 		if (~ismember(theFixedH1CellIndex, RGCMosaicConstructor.constants.validPackerDacey2002H1CellIndices))
		 			theFixedH1CellIndex  = [];
		 		end
		 	end

		 	% Try for L-cone dominated substring
		 	targetSubString = 'LconeDominated_';
		 	idx = strfind(theSelectedOptimizationResultsFileName, targetSubString);
		 	p1 = idx+numel(targetSubString);
		 	if (~isempty(idx))
		 		theRFcenterDominantConeType = cMosaic.LCONE_ID;
		 		targetSubString = 'coneRFcenter';
		 		p2 = strfind(theSelectedOptimizationResultsFileName, targetSubString);
		 		if (~isempty(p2))
		 			theRFcenterConesNum = str2num(theSelectedOptimizationResultsFileName(p1:(p2-1)));
		 		end
		 	end

		 	% Try for M-cone dominated substring
		 	targetSubString = 'MconeDominated_';
		 	idx = strfind(theSelectedOptimizationResultsFileName, targetSubString);
		 	p1 = idx+numel(targetSubString);
		 	if (~isempty(idx))
		 		theRFcenterDominantConeType = cMosaic.MCONE_ID;
		 		targetSubString = 'coneRFcenter';
		 		p2 = strfind(theSelectedOptimizationResultsFileName, targetSubString);
		 		if (~isempty(p2))
		 			theRFcenterConesNum = str2num(theSelectedOptimizationResultsFileName(p1:(p2-1)));
		 		end
		 	end

		 	% Try for S-cone dominated substring
		 	targetSubString = 'SconeDominated_';
		 	idx = strfind(theSelectedOptimizationResultsFileName, targetSubString);
		 	p1 = idx+numel(targetSubString);
		 	if (~isempty(idx))
		 		theRFcenterDominantConeType = cMosaic.SCONE_ID;
		 		targetSubString = 'coneRFcenter';
		 		p2 = strfind(theSelectedOptimizationResultsFileName, targetSubString);
		 		if (~isempty(p2))
		 			theRFcenterConesNum = str2num(theSelectedOptimizationResultsFileName(p1:(p2-1)));
		 		end
		 	end

		end % if (~isempty(theSelectedOptimizationResultsFileName))
	else

		while (isempty(theFixedH1CellIndex)) || ((~isempty(theFixedH1CellIndex ))&&(~ismember(theFixedH1CellIndex, RGCMosaicConstructor.constants.validPackerDacey2002H1CellIndices)))
			validH1cellIndices = RGCMosaicConstructor.constants.validPackerDacey2002H1CellIndices
		    theFixedH1CellIndex  = input('Enter the target H1 horizontal cell index:');
		end
	end
	
end
