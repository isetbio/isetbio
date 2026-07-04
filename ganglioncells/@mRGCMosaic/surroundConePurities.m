function [surroundConePurities, centerConeDominances, centerNumerosities, centerConePurities] = surroundConePurities(obj, ...
	theRGCindices, surroundConeSelection)

	assert(ismember(surroundConeSelection, {'all', 'surround pooling weights > center pooling weights'}), ...
		'surroundConeSelection must be set to either ''all'', or ''surround pooling weights > center pooling weights''. ''%s'' is not a valid option', surroundConeSelection);

	centerConeDominances = zeros(1, numel(theRGCindices));
	centerNumerosities = zeros(1, numel(theRGCindices));
	centerConePurities = zeros(1, numel(theRGCindices));
    surroundConePurities = nan(1, numel(theRGCindices));

	parfor iRGC = 1:numel(theRGCindices)
		theRGCindex = theRGCindices(iRGC);

		% Determine the RF center cone dominance
		rfCenterStats = obj.singleCellConnectivityStats(theRGCindex, 'center', ...
			'minConeWeightIncluded', 0.1);

		% Determine surround cones 
		if (strcmp(surroundConeSelection, 'all'))
			rfSurroundStats = obj.singleCellConnectivityStats(theRGCindex, 'surround');
		else
			rfSurroundStats = obj.singleCellConnectivityStats(theRGCindex, 'surround-center');
		end

        if (isempty(rfSurroundStats))
            totalLMsurroundSignal = nan;
        else
		    totalLMsurroundSignal = ...
			    rfSurroundStats.netWeights(cMosaic.LCONE_ID) + rfSurroundStats.netWeights(cMosaic.MCONE_ID);
        end

        totalLMcenterSignal = ...
			rfCenterStats.netWeights(cMosaic.LCONE_ID) + rfCenterStats.netWeights(cMosaic.MCONE_ID);

		if (rfCenterStats.dominantConeType == cMosaic.LCONE_ID)
			% surround purity score for L-center RFs: M/(L+M)
            if (~isnan(totalLMsurroundSignal))
			    surroundConePurities(iRGC) = rfSurroundStats.netWeights(cMosaic.MCONE_ID)/totalLMsurroundSignal;
            end
			centerConePurities(iRGC) = rfCenterStats.netWeights(cMosaic.LCONE_ID) / totalLMcenterSignal;
		else
			% surround purity score for M-center RFs: L/(L+M)
            if (~isnan(totalLMsurroundSignal))
			    surroundConePurities(iRGC) = rfSurroundStats.netWeights(cMosaic.LCONE_ID)/totalLMsurroundSignal;
            end
			centerConePurities(iRGC) = rfCenterStats.netWeights(cMosaic.MCONE_ID) / totalLMcenterSignal;
		end

		centerConeDominances(iRGC) = rfCenterStats.dominantConeType;
		centerNumerosities(iRGC) = rfCenterStats.inputConesNum;
	end % iRGC
end
