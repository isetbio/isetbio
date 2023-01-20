function theMajorityCenterConeType = majorityCenterConeType(obj, theRGCindex)
    if (~isempty(obj.rgcRFcenterConeConnectivityMatrix))
        % Retrieve this cell's # of center cone indices
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
    else
        % Retrieve this cell's # of center cone indices
        connectivityVector = full(squeeze(obj.rgcRFcenterConePoolingMatrix(:, theRGCindex)));
    end

    indicesOfCenterCones = find(connectivityVector > 0.0001);

    coneTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
    lmsConesNum = zeros(1, numel(coneTypes));
    for iConeType = 1:numel(coneTypes)
        lmsConesNum(iConeType) = numel(find(obj.inputConeMosaic.coneTypes(indicesOfCenterCones) == coneTypes(iConeType)));
    end
    [~,idx] = max(lmsConesNum);
    theMajorityCenterConeType = coneTypes(idx);
end
