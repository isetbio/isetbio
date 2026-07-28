function [optimizationResults, theMRGCMosaic, theTargetRGCindex, targetVisualSTFparams, optimizationResultsFileName, matchedAttributes] = ...
		 importOptimizationResults(theSurroundConnectedMRGCMosaicFullFileName, theRFcenterConesNum, theRFcenterDominantConeType, initialSurroundOptimizationValuesSource, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('tryAlternateH1CellIndexAsLastResort', false, @islogical);
    p.parse(varargin{:});
    tryAlternateH1CellIndexAsLastResort = p.Results.tryAlternateH1CellIndexAsLastResort;

    % Try to import with desired H1CellIndex
    [optimizationResults, theMRGCMosaic, theTargetRGCindex, targetVisualSTFparams, optimizationResultsFileName, matchedAttributes] = ...
		 importOptimizationResultsForSpecificH1CellIndex(theSurroundConnectedMRGCMosaicFullFileName, theRFcenterConesNum, theRFcenterDominantConeType, ...
         initialSurroundOptimizationValuesSource);

    % Optimization results not found. Try with alternate H1CellIndex
    if (isempty(optimizationResults)) && (tryAlternateH1CellIndexAsLastResort)

        % Retrieve targetH1CellIndex from theSurroundConnectedMRGCMosaicFullFileName
        targetString = 'FixedCellIndex';
        idx = strfind(theSurroundConnectedMRGCMosaicFullFileName, targetString);
        targetH1CellIndex = str2num(strrep(theSurroundConnectedMRGCMosaicFullFileName(idx+(0:numel(targetString))), targetString, ''));
        targetH1CellString = sprintf('%s%d', targetString, targetH1CellIndex);

        % try with the alterateH1CellIndices
        alternateH1CellIndices = setdiff(RGCMosaicConstructor.constants.validPackerDacey2002H1CellIndices, targetH1CellIndex);
        for idx = 1:numel(alternateH1CellIndices)
            if (isempty(optimizationResults))
                
                % Generate optimization results filename for alternateH1CellIndex
                alternateH1CellIndex = alternateH1CellIndices(idx);
                alternateH1CellString = sprintf('%s%d', targetString, alternateH1CellIndex);
                fprintf('Trying to find optimization results for alternate H1cellIndex: %d instead of %d\n', alternateH1CellIndex, targetH1CellIndex)
                theAlternateH1CellIndexSurroundConnectedMRGCMosaicFullFileName = strrep(theSurroundConnectedMRGCMosaicFullFileName, targetH1CellString, alternateH1CellString);
                
                % Try loading this optimization file 
                [optimizationResults, theMRGCMosaic, theTargetRGCindex, targetVisualSTFparams, optimizationResultsFileName, matchedAttributes] = ...
		             importOptimizationResultsForSpecificH1CellIndex(theAlternateH1CellIndexSurroundConnectedMRGCMosaicFullFileName, theRFcenterConesNum, theRFcenterDominantConeType, ...
                     initialSurroundOptimizationValuesSource);

                if (~isempty(optimizationResults))
                    fprintf('Found optimization results for alternate H1cellIndex: %d instead of %d.\n', alternateH1CellIndex, targetH1CellIndex);
                end
            end
        end
    end
end

function [optimizationResults, theMRGCMosaic, theTargetRGCindex, targetVisualSTFparams, optimizationResultsFileName, matchedAttributes] = ...
		 importOptimizationResultsForSpecificH1CellIndex(theSurroundConnectedMRGCMosaicFullFileName, theRFcenterConesNum, theRFcenterDominantConeType, ...
         initialSurroundOptimizationValuesSource)

	% Assemble optimizationResultsFilename based on employed MRGC mosaic + target input cones num + target input cone dominance
    optimizationResultsFileName = RGCMosaicConstructor.filepathFor.optimizationResults(...
    	theSurroundConnectedMRGCMosaicFullFileName, ...
    	theRFcenterConesNum, ...
    	theRFcenterDominantConeType);

    switch (theRFcenterDominantConeType)
		case cMosaic.LCONE_ID
			fprintf('Trying to find a previously optimized pooling model for %d-cone RF center with L-cone dominance...\n', ...
			theRFcenterConesNum);
		case cMosaic.MCONE_ID
			fprintf('Trying to find a previously optimized pooling model for %d-cone RF center with M-cone dominance...\n', ...
			theRFcenterConesNum);
		case cMosaic.SCONE_ID
			fprintf('Trying to find a previously optimized pooling model for %d-cone RF center with S-cone dominance...\n', ...
			theRFcenterConesNum);
	end % switch

    try 
    	% Try for exact match
	   	load(optimizationResultsFileName, ...
	   			'theMRGCMosaic', ...
		   		'targetVisualSTFparams', ...
		        'theTargetRGCindex', ...
		        'theTargetRGCRFcenterConnectivityStats', ...
		        'optimizationResults');

	   	matchedAttributes.theRFcenterDominantConeType = theRFcenterDominantConeType;
	   	matchedAttributes.theRFcenterConesNum = theRFcenterConesNum;
        fprintf('\nFound optimization results file with matched cone dominance and numerosity.\n\t %s \nWill use it.\n', optimizationResultsFileName);
		   		
    catch % catch #1
	   	fprintf('\nThe requested optimization results file \n\t %s \n was **NOT** found.\n', optimizationResultsFileName);

	   	if (~strcmp(initialSurroundOptimizationValuesSource, 'imported closest match'))
            fprintf('\nNot importing any previous results\n');
	   		optimizationResults = [];
	   		theMRGCMosaic = [];
	   	    theTargetRGCindex = [];
	   	    targetVisualSTFparams = [];
	   	    matchedAttributes = [];
            return;
        end

   		% Try for alternate rf center cone dominance
   		if (theRFcenterDominantConeType == cMosaic.LCONE_ID)
   			theAlternateRFcenterDominantConeType = cMosaic.MCONE_ID;
   		else
   			theAlternateRFcenterDominantConeType = cMosaic.LCONE_ID;
   		end

   		optimizationResultsFileNameAlternateConeDominance = RGCMosaicConstructor.filepathFor.optimizationResults(...
			theSurroundConnectedMRGCMosaicFullFileName, ...
			theRFcenterConesNum, ...
			theAlternateRFcenterDominantConeType);

   		fprintf('--> Trying with alternate cone dominance\n');

   		try 
    		% Try for match of alternate cone dominance
	   		load(optimizationResultsFileNameAlternateConeDominance, ...
	   			'theMRGCMosaic', ...
		   		'targetVisualSTFparams', ...
		        'theTargetRGCindex', ...
		        'theTargetRGCRFcenterConnectivityStats', ...
		        'optimizationResults');

	   		matchedAttributes.theRFcenterDominantConeType = theAlternateRFcenterDominantConeType;
	   		matchedAttributes.theRFcenterConesNum = theRFcenterConesNum;
            fprintf('\nFound optimization results file with alternate cone dominance\n\t %s \nWill use it.\n', optimizationResultsFileNameAlternateConeDominance);
	   		
        catch % catch #2 didnt find one

	   		fprintf('\nThe alternate cone dominance optimization results file \n\t %s \n was **NOT** found.\n', optimizationResultsFileNameAlternateConeDominance);
	   		fprintf('Trying with alternate cone numerosity\n');

	   		% Try for alternate rf center numerosity
	   		if (theRFcenterConesNum == 1)
	   			theAlternateRFcenterConesNum = 2;
	   		else
	   			theAlternateRFcenterConesNum = theRFcenterConesNum-1;
	   		end
	   			
	   		optimizationResultsFileNameAlternateConeNumerosity = RGCMosaicConstructor.filepathFor.optimizationResults(...
				theSurroundConnectedMRGCMosaicFullFileName, ...
				theAlternateRFcenterConesNum, ...
				theRFcenterDominantConeType);

	   		try 
	    		% Try for match of alternate cone numerosity
		   		load(optimizationResultsFileNameAlternateConeNumerosity, ...
		   			'theMRGCMosaic', ...
			   		'targetVisualSTFparams', ...
			        'theTargetRGCindex', ...
			        'theTargetRGCRFcenterConnectivityStats', ...
			        'optimizationResults');

		   		matchedAttributes.theRFcenterDominantConeType = theRFcenterDominantConeType;
	   			matchedAttributes.theRFcenterConesNum = theAlternateRFcenterConesNum;

                fprintf('\nFound optimization results file with alternate cone numerosity \n\t %s \nWill use it.\n',optimizationResultsFileNameAlternateConeNumerosity);
	   		

            catch % catch #3
	   			fprintf('\n\tThe alternate cone numerosity optimization results file \n %s \n\t was **NOT** found.\n', optimizationResultsFileNameAlternateConeNumerosity);
	   	
	   			optimizationResults = [];
   				theMRGCMosaic = [];
   	    		theTargetRGCindex = [];
   	    		targetVisualSTFparams = [];
   	    		matchedAttributes = [];
	   		end % try-catch #3
	    end % try-catch #2
	end % try-catch #1
end
