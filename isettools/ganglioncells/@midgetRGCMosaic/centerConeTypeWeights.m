function [theCenterConeTypeWeights, theCenterConeTypeNum] = centerConeTypeWeights(obj, theRGCindex)
    if (~isempty(obj.rgcRFcenterConeConnectivityMatrix))
        % Retrieve this cell's # of center cone indices
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
    else
        % Retrieve this cell's # of center cone indices
        connectivityVector = full(squeeze(obj.rgcRFcenterConePoolingMatrix(:, theRGCindex)));
    end
    indicesOfCenterCones = find(connectivityVector > 0.0001);

    weightsOfCenterCones = connectivityVector(indicesOfCenterCones);
    typesOfCenterCones = obj.inputConeMosaic.coneTypes(indicesOfCenterCones);

    coneTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
    theCenterConeTypeWeights = zeros(1, numel(coneTypes));
    theCenterConeTypeNum = zeros(1, numel(coneTypes));

    for iConeType = 1:numel(coneTypes)
        theConeType = coneTypes(iConeType);
        idx = find(typesOfCenterCones == theConeType);
        theCenterConeTypeWeights(theConeType) = sum(weightsOfCenterCones(idx));
        theCenterConeTypeNum(theConeType) = numel(idx);
    end
    
end
