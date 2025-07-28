function bakeSurroundConeConnectivityMatrixAndFreeze(obj, rgcRFsurroundConeConnectivityMatrix, ...
	rfSurroundConnectivityParams, surroundVarianceInComputeReadyMosaic)

	if (~isempty(obj.rgcRFsurroundConeConnectivityMatrix))
		error('rgcRFsurroundConeConnectivityMatrix is not empty.')
	end

	if (~isempty(obj.rfSurroundConnectivityParams))
		error('rfSurroundConnectivityParams is not empty.')
	end

	if (~isempty(obj.surroundVarianceInComputeReadyMosaic))
		error('surroundVarianceInComputeReadyMosaic is not empty.')
	end
	
	% Set the gains to 1/sum(center cone weights)
	responseGains = zeros(1,obj.rgcsNum);
	parfor iRGC = 1:obj.rgcsNum
		theInputConeIndices = obj.singleCellConnectivityStats(iRGC, 'center', ...
			'minConeWeightIncluded', 0.001, ...
			'inputConeIndicesOnly', true, ...
			'warnIfCenterConeInputNumerosityDiffersFromExclusiveOne', false);
		inputConeWeights = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(theInputConeIndices, iRGC)));
		responseGains(iRGC) = 1.0/sum(inputConeWeights);
	end

	obj.responseGains  = responseGains;
	obj.rgcRFsurroundConeConnectivityMatrix = rgcRFsurroundConeConnectivityMatrix;
	obj.rfSurroundConnectivityParams = rfSurroundConnectivityParams;
	obj.surroundVarianceInComputeReadyMosaic = surroundVarianceInComputeReadyMosaic;
end

