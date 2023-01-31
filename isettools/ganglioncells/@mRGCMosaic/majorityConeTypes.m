function [centerSubregionMajorityConeTypes, theCenterConeTypesNums] = majorityConeTypes(obj, theRGCindices)
    
    coneTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];

    theCenterConeTypesNums = zeros(numel(theRGCindices), numel(coneTypes));

    parfor i = 1:numel(theRGCindices)
        theRGCindex = theRGCindices(i);
        connectivityVector = full(squeeze(obj.centerConePoolingMatrix(:, theRGCindex)));
        indicesOfCenterCones = find(connectivityVector > 0.0001);
    
        lmsConesNum = zeros(1, numel(coneTypes));
        theCenterConeTypesNumsForThisRGC = zeros(1, numel(coneTypes));

        for iConeType = 1:numel(coneTypes)
            theConeType = coneTypes(iConeType);
            lmsConesNum(iConeType) = numel(find(obj.inputConeMosaic.coneTypes(indicesOfCenterCones) == theConeType));
            theCenterConeTypesNumsForThisRGC(theConeType) = lmsConesNum(iConeType);
        end

        theCenterConeTypesNums(i,:) = theCenterConeTypesNumsForThisRGC;
        [~,idx] = max(lmsConesNum);
        centerSubregionMajorityConeTypes(i) = coneTypes(idx);
    end

end