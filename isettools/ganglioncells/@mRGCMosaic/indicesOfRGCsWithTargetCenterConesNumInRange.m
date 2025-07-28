function theRGCIndices = indicesOfRGCsWithTargetCenterConesNumInRange(obj, targetCenterConesNumRange, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('minConeWeightIncluded', mRGCMosaic.sensitivityAtPointOfOverlap, @isscalar);
    p.parse(varargin{:});
    minConeWeightIncluded = p.Results.minConeWeightIncluded;

    centerConesNum = zeros(1, obj.rgcsNum);

    parfor iRGC = 1:obj.rgcsNum
        s = obj.singleCellConnectivityStats(...
            iRGC, 'center', ...
            'minConeWeightIncluded', minConeWeightIncluded);
        centerConesNum(iRGC) = s.inputConesNum;
    end

    if (numel(targetCenterConesNumRange) == 2)
        theRGCIndices = find(...
            (centerConesNum >= targetCenterConesNumRange(1)) & ...
            (centerConesNum <= targetCenterConesNumRange(2)));
    elseif (numel(targetCenterConesNumRange) == 1)
        theRGCIndices = find(centerConesNum == targetCenterConesNumRange);
    else
        error('targetCenterConesNumRange must be either a 1- or a 2-element vector');
    end
    
end