function [theCenterConeTypeWeights, theCenterConeTypeNum, theMajorityConeType] = centerConeTypeWeights(obj, theRGCindex)
    % Retrieve this cell's # of center cone indices
    connectivityVector = full(squeeze(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
    
    indicesOfCenterCones = find(connectivityVector > 0.0001);
    weightsOfCenterCones = connectivityVector(indicesOfCenterCones);
    typesOfCenterCones = obj.theRGCMosaic.inputConeMosaic.coneTypes(indicesOfCenterCones);

    coneTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
    theCenterConeTypeWeights = zeros(1, numel(coneTypes));
    theCenterConeTypeNum = zeros(1, numel(coneTypes));

    for iConeType = 1:numel(coneTypes)
        theConeType = coneTypes(iConeType);
        idx = find(typesOfCenterCones == theConeType);
        theCenterConeTypeWeights(theConeType) = sum(weightsOfCenterCones(idx));
        theCenterConeTypeNum(theConeType) = numel(idx);
    end
    
    % Determine majority cone type based on weights
    thresholdRatioForEqualWeights = 0.9;
    [~,idx] = sort(theCenterConeTypeWeights, 'descend');
    theMajorityConeType = coneTypes(idx(1));
    if (numel(find(theCenterConeTypeWeights>0)) > 1)
        ratio = theCenterConeTypeWeights(idx(2))/theCenterConeTypeWeights(idx(1));
        if (ratio  > thresholdRatioForEqualWeights)
            % nearly equal weights, majorityConeType = nan
            theMajorityConeType = nan;
        end
    end

end