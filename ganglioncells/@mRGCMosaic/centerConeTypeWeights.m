function [theCenterConeTypeWeights, theCenterConeTypeNum, theMajorityConeType, theCenterConeTypes, theCenterConeIndices] = centerConeTypeWeights(obj, theRGCindex)
    % Retrieve this cell's # of center cone indices
    if (isempty(obj.rgcRFcenterConePoolingMatrix))
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
    else
        connectivityVector = full(squeeze(obj.rgcRFcenterConePoolingMatrix(:, theRGCindex)));
    end

    theCenterConeIndices = find(connectivityVector > 0.0001);
    weightsOfCenterCones = connectivityVector(theCenterConeIndices);
    typesOfCenterCones = obj.inputConeMosaic.coneTypes(theCenterConeIndices);

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

    theCenterConeTypeWeights = theCenterConeTypeWeights(idx);
    theCenterConeTypeNum = theCenterConeTypeNum(idx);
    theMajorityConeType = coneTypes(idx(1));
    theCenterConeTypes = coneTypes(idx);

    if (numel(find(theCenterConeTypeWeights>0)) > 1)
        ratio = theCenterConeTypeWeights(2)/theCenterConeTypeWeights(1);
        if (ratio  > thresholdRatioForEqualWeights)
            % nearly equal weights, majorityConeType = nan
            theMajorityConeType = nan;
        end
    end

end