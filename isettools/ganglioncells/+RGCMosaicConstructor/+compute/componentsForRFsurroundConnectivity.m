function dataOut = componentsForRFsurroundConnectivity(theInputConeMosaicSTFResponsesFullFileName, ...
	targetVisualSTFmodificationsFromMeanValuesParamsStruct, ...
	centerConeNumerositiesToOptimize, centerConeDominanceToOptimize, positionDegsToOptimize, ...
	doNotWorryAboutMaximizingTargetRFcenterLMconeRatio, ...
	importInputConeMosaicSTFresponses)

	assert(ismember(centerConeDominanceToOptimize, [cMosaic.LCONE_ID cMosaic.MCONE_ID]), ...
		'centerConeDominanceToOptimize must be either cMosaic.LCONE_ID (1) or cMosaic.MCONE_ID (2)');

	dataOut = [];

    if (importInputConeMosaicSTFresponses)
    	fprintf('Loading mRGC mosaic, simulation params, and input cone mosaic responses from:\n%s\n', theInputConeMosaicSTFResponsesFullFileName);

		% Load the input cone mosaic STF responses and the corresponding STFparams struct 
		load(theInputConeMosaicSTFResponsesFullFileName, ...
			'theMRGCMosaic', ...
			'STFparamsStruct', ...
	    	'stimParams', ...
	    	'theInputConeMosaicSTFresponses');
		% Assemble dataOut
		dataOut.theMRGCMosaic = theMRGCMosaic;
		dataOut.stimParams = stimParams;
		dataOut.theInputConeMosaicSTFresponses = theInputConeMosaicSTFresponses;
		dataOut.theOptimizationPositionDegs = positionDegsToOptimize;
	else
		fprintf('Loading mRGC mosaic to determine the cone numerosities optimized at current position.\n%s\n', theInputConeMosaicSTFResponsesFullFileName);
		load(theInputConeMosaicSTFResponsesFullFileName, ...
			'theMRGCMosaic', ...
			'STFparamsStruct', ...
	    	'stimParams');
		% Assemble dataOut
		dataOut.theMRGCMosaic = theMRGCMosaic;
		dataOut.stimParams = stimParams;
		dataOut.theOptimizationPositionDegs = positionDegsToOptimize;
	end


	% Target visual STF params at the optimization position (STFparamsStruct.positionDegs)
    temporalEquivalentEccDegsInHumanRetina = theMRGCMosaic.temporalEquivalentEccentricityForEccXYDegs(dataOut.theOptimizationPositionDegs);
    radialTemporalEquivalentEccDegsInHumanRetina = sqrt(sum(temporalEquivalentEccDegsInHumanRetina.^2,2));

    % Convert eccentricity in human retina to eccentricity in macaque retina
    radialTemporalEquivalentEccDegsInMacaqueRetina = ...
     	RGCMosaicConstructor.helper.convert.eccentricityDegsBetweenSpecies(...
     		'HumanRetinaToMacaqueRetina', radialTemporalEquivalentEccDegsInHumanRetina);

    % Mean values from the Croner&Kaplan data
    meanSurroundToCenterRcRatio = RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio;
    meanSurroundToCenterIntegratedSensitivityRatio = RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(radialTemporalEquivalentEccDegsInMacaqueRetina);

    % Apply modifications from the mean visualSTF values
    fprintf('Modifiers from mean visual STF values. Rs/Rc: %f, intS/C: %f\n', ...
    	targetVisualSTFmodificationsFromMeanValuesParamsStruct.surroundToCenterRcRatioMultiplier, ...
    	targetVisualSTFmodificationsFromMeanValuesParamsStruct.surroundToCenterIntegratedSensitivityRatioMultiplier);

    targetSurroundToCenterRcRatio = meanSurroundToCenterRcRatio * ...
    	targetVisualSTFmodificationsFromMeanValuesParamsStruct.surroundToCenterRcRatioMultiplier;

    targetSurroundToCenterIntegratedSensitivityRatio =  meanSurroundToCenterIntegratedSensitivityRatio * ...
    	targetVisualSTFmodificationsFromMeanValuesParamsStruct.surroundToCenterIntegratedSensitivityRatioMultiplier;

    % Get C&K data for the equivalent eccentricity in the macaque retina
	dataOut.targetVisualSTFparams = struct(...
        'surroundToCenterRcRatio', targetSurroundToCenterRcRatio, ...
        'surroundToCenterIntegratedSensitivityRatio', targetSurroundToCenterIntegratedSensitivityRatio);

	% Determine unique # of center cones num in the mosaic
	[allRGCCenterConesNum, allRGCCenterDominantConeTypes, allRGCCenterRelativeConeWeights] = ...
		theMRGCMosaic.allRFcenterConnectivityStats();

  	% Find the indices of RGCs whose cone center cone dominance matches centerConeDominanceToOptimize
	idx = find(allRGCCenterDominantConeTypes == centerConeDominanceToOptimize);
	allConeNumerosities = unique(allRGCCenterConesNum(idx));

    % Downsample the population of centerConesNum 
    dataOut.uniqueCenterConesNum = RGCMosaicConstructor.compute.downsampledRFcenterConeNumerosityPopulation(allConeNumerosities);

    dataOut.rgcsNumWithRFcenterNumerosity = zeros(1, max(dataOut.uniqueCenterConesNum));

	for iCenterConesNum = 1:numel(dataOut.uniqueCenterConesNum)
		% Find the RGCs with this many center cones
		theTargetCenterConesNum = dataOut.uniqueCenterConesNum(iCenterConesNum);
		theRGCIndicesWithTargetCenterConesNum = theMRGCMosaic.indicesOfRGCsWithTargetCenterConesNumInRange(theTargetCenterConesNum);

		% Find indices of RGCs with desired # of center cones and desired cone dominance
		idx = find(allRGCCenterDominantConeTypes(theRGCIndicesWithTargetCenterConesNum) == centerConeDominanceToOptimize);
		theRGCIndicesWithTargetCenterConesNum = theRGCIndicesWithTargetCenterConesNum(idx);

		% Sort these rgcs based on their distance to theOptimizationPositionDegs
		d2 = sum((bsxfun(@minus, theMRGCMosaic.rgcRFpositionsDegs(theRGCIndicesWithTargetCenterConesNum,:), ...
			                     dataOut.theOptimizationPositionDegs)) .^ 2, 2);

		% OLD METHOD: Just pick the closest mRGC. But this may not have the maximal L/M cone ratio.
		[~,idx] = min(d2(:));
    	dataOut.theTargetRGCindex(iCenterConesNum) = theRGCIndicesWithTargetCenterConesNum(idx);

    	if (doNotWorryAboutMaximizingTargetRFcenterLMconeRatio)
    		% Do not worry about maximizing the L/M cone ratio
    	else
	    	if (theTargetCenterConesNum > 1)
	    		totalLMweight = ...
	    			allRGCCenterRelativeConeWeights(dataOut.theTargetRGCindex(iCenterConesNum),2) + ...
	    			allRGCCenterRelativeConeWeights(dataOut.theTargetRGCindex(iCenterConesNum),1);

				if (centerConeDominanceToOptimize == cMosaic.LCONE_ID)
		    		oldMethodLMconeRatio = allRGCCenterRelativeConeWeights(dataOut.theTargetRGCindex(iCenterConesNum),1) / totalLMweight;
		    	else
		    		oldMethodLMconeRatio = allRGCCenterRelativeConeWeights(dataOut.theTargetRGCindex(iCenterConesNum),2) / totalLMweight;
		    	end

		    	% NEW METHOD: SORT the mRGCs, and pick from the 20 closest ones, the one with the maximal L/M cone ratio.
		    	[~, distanceSortedIndices] = sort(d2(:), 'ascend');
		    	theRGCIndicesWithTargetCenterConesNum = theRGCIndicesWithTargetCenterConesNum(distanceSortedIndices);

		    	% Search the closest 64 RGCs to find the one with the largest L/M ratio (if centerConeDominanceToOptimize == cMosaic.LCONE_ID)
		    	% or the largest M/L ratio (if centerConeDominanceToOptimize == cMosaic.MCONE_ID)
		    	Nnearby = min([64 numel(distanceSortedIndices)]);
		    	candidateRGCindices = theRGCIndicesWithTargetCenterConesNum(1:Nnearby);
		    	totalCandidateLMweights = ...
		    		allRGCCenterRelativeConeWeights(candidateRGCindices,2) + ...
		    		allRGCCenterRelativeConeWeights(candidateRGCindices,1);
		    	if (centerConeDominanceToOptimize == cMosaic.LCONE_ID)
		    		theLMconeRatios = allRGCCenterRelativeConeWeights(candidateRGCindices,1) ./ totalCandidateLMweights;
		    	else
		    		theLMconeRatios = allRGCCenterRelativeConeWeights(candidateRGCindices,2) ./ totalCandidateLMweights;
		    	end
		    	[maxLMratio,idx] = max(theLMconeRatios);
		    	if (maxLMratio > oldMethodLMconeRatio)
		    		fprintf('>>>> Found nearby RGC with %d inputs (dominance: %d) with better LM ratio: %f vs %f\n', theTargetCenterConesNum, centerConeDominanceToOptimize, maxLMratio, oldMethodLMconeRatio);
		    		dataOut.theTargetRGCindex(iCenterConesNum) = candidateRGCindices(idx);
		    	end
		    end
		end

        dataOut.rgcsNumWithRFcenterNumerosity(theTargetCenterConesNum) = numel(theRGCIndicesWithTargetCenterConesNum);

        fprintf('There are %d mRGCs with cone dominance %d and %d center cones. Optimized RGC: #%d, located at (%2.2f,%2.2f) which is closest to the optimization position (%2.2f,%2.2f).\n', ...
			numel(theRGCIndicesWithTargetCenterConesNum), centerConeDominanceToOptimize, theTargetCenterConesNum, ...
            dataOut.theTargetRGCindex(iCenterConesNum), ...
            theMRGCMosaic.rgcRFpositionsDegs(dataOut.theTargetRGCindex(iCenterConesNum),1), ...
            theMRGCMosaic.rgcRFpositionsDegs(dataOut.theTargetRGCindex(iCenterConesNum),2), ...
            dataOut.theOptimizationPositionDegs(1), dataOut.theOptimizationPositionDegs(2));

	end % iCenterConesNum

    dataOut.centerConeNumerositiesToOptimize = centerConeNumerositiesToOptimize;
end


