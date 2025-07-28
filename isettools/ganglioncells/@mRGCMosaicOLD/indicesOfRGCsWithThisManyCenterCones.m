function indicesOfRGCsIdentified = indicesOfRGCsWithThisManyCenterCones(obj, targetCenterConesNum)
    if (targetCenterConesNum > 0)

        if (isempty(obj.rgcRFcenterConePoolingMatrix))
            centerConesNum = full(sum(obj.rgcRFcenterConeConnectivityMatrix,1));
        else
            centerConesNum = full(sum(obj.rgcRFcenterConePoolingMatrix,1));
        end

        indicesOfRGCsIdentified = find(centerConesNum == targetCenterConesNum);

    else
        indicesOfRGCsIdentified  = [];
    end
end
