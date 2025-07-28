function [allRGCCenterConesNum, allRGCCenterDominantConeTypes, allRGCCenterRelativeConeWeights] = allRFcenterConnectivityStats(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('minConeWeightIncluded', mRGCMosaic.sensitivityAtPointOfOverlap, @isscalar);
    p.parse(varargin{:});
    minConeWeightIncluded = p.Results.minConeWeightIncluded;

	allRGCCenterConesNum = zeros(obj.rgcsNum,1);
	allRGCCenterDominantConeTypes = zeros(obj.rgcsNum,1);
	allRGCCenterRelativeConeWeights = zeros(obj.rgcsNum,3);

    parfor iRGC = 1:obj.rgcsNum
        s = obj.singleCellConnectivityStats(iRGC, 'center', ...
            'minConeWeightIncluded', minConeWeightIncluded);

        if (~isempty(s))
            allRGCCenterConesNum(iRGC) = s.inputConesNum;
            allRGCCenterDominantConeTypes(iRGC) = s.dominantConeType;
            allRGCCenterRelativeConeWeights(iRGC,:) = s.netWeights;
        end
    end

end