function centerSubregionMajorityConeTypes  = majorityConeTypes(obj, theRGCindices)
    
    theMajorityCenterConeTypes = zeros(1, numel(theRGCindices));

    parfor i = 1:numel(theRGCindices)
        theRGCindex = theRGCindices(i);
        connectivityVector = full(squeeze(obj.centerConePoolingMatrix(:, theRGCindex)));
        indicesOfCenterCones = find(connectivityVector > 0.0001);
    
        coneTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
        lmsConesNum = zeros(1, numel(coneTypes));
        for iConeType = 1:numel(coneTypes)
            lmsConesNum(iConeType) = numel(find(obj.inputConeMosaic.coneTypes(indicesOfCenterCones) == coneTypes(iConeType)));
        end
        [~,idx] = max(lmsConesNum);
        centerSubregionMajorityConeTypes(i) = coneTypes(idx);
    end

end