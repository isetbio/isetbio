function [targetRGCindices, surroundConePurities, centerConeDominances, ...
	      centerConeNumerosities, centerConePurities] = indicesOfRGCsWithinTargetedPropertyRanges(obj, ...
			targetedCenterConeNumerosityRange, ...
			targetedSurroundPurityRange, ...
			targetedRadialEccentricityRange, ...
			targetedCenterPurityRange)


        if ((isempty(targetedCenterConeNumerosityRange))&& ...
        	(isempty(targetedSurroundPurityRange))&& ...
        	(isempty(targetedRadialEccentricityRange)) && ...
        	(isempty(targetedCenterPurityRange)) )
            	targetRGCindices = 1:obj.rgcsNum;
            return;
        end

		% Compute surround cone purities and center cone dominances
		% based on cones whose surround pooling weights > center pooling weights
		surroundConeSelection = 'surround pooling weights > center pooling weights';
		[surroundConePurities, centerConeDominances, centerConeNumerosities, centerConePurities] = ...
			obj.surroundConePurities(1:obj.rgcsNum, surroundConeSelection);

		% Compute radial eccentricities
		radialEccentricities = sqrt(sum(obj.rgcRFpositionsDegs(1:obj.rgcsNum,:).^2,2));
	
		if (~isempty(targetedCenterConeNumerosityRange)) 
	   		t1 = ( (centerConeNumerosities(:) >= targetedCenterConeNumerosityRange(1)) & ...
			       (centerConeNumerosities(:) <= targetedCenterConeNumerosityRange(2)) );
	   	else
	   		t1 = true(numel(centerConeNumerosities),1);
	   	end

	   	if (~isempty(targetedSurroundPurityRange)) 
	   		t2 = ( (surroundConePurities(:) >= targetedSurroundPurityRange(1)) & ...
			       (surroundConePurities(:) <= targetedSurroundPurityRange(2)) );
	   	else
	   		t2 = true(numel(surroundConePurities),1);
	   	end 

	   	if (~isempty(targetedRadialEccentricityRange))
	   		t3 = ((radialEccentricities(:) >= targetedRadialEccentricityRange(1)) & ...
		          (radialEccentricities(:) <= targetedRadialEccentricityRange(2)) );
	   	else
	   		t3 = true(numel(radialEccentricities),1);
	   	end

	   	if (~isempty(targetedCenterPurityRange)) 
	   		t4 = ( (centerConePurities(:) >= targetedCenterPurityRange(1)) & ...
			       (centerConePurities(:) <= targetedCenterPurityRange(2)) );
	   	else
	   		t4 = true(numel(centerConePurities),1);
	   	end 

		targetRGCindices = find(t1 & t2 &t3 & t4);

		surroundConePurities = surroundConePurities(targetRGCindices);
		centerConeDominances =  centerConeDominances(targetRGCindices);
		centerConeNumerosities = centerConeNumerosities(targetRGCindices);
		centerConePurities = centerConePurities(targetRGCindices);
end