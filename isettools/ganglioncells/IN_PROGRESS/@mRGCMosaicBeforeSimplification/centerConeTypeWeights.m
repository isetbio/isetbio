function [allCellsCenterConeTypeWeights, allCellsCenterConeTypeNum, allCellsMajorityConeType, allCellsCenterConesNum] = ...
            centerConeTypeWeights(obj, theRGCindices)

    coneTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];

    allCellsCenterConeTypeWeights = zeros(numel(theRGCindices), numel(coneTypes));
    allCellsCenterConeTypeNum = zeros(numel(theRGCindices), numel(coneTypes));
    allCellsMajorityConeType = zeros(numel(theRGCindices),1);
    allCellsCenterConesNum = zeros(numel(theRGCindices),1);
    
    parfor iRGC = 1:numel(theRGCindices)
        theRGCindex = theRGCindices(iRGC);

        % Retrieve this cell's # of center cone indices
        connectivityVector = full(squeeze(obj.centerConePoolingMatrix(:, theRGCindex)));

        indicesOfCenterCones = find(connectivityVector > 0.0001);
        weightsOfCenterCones = connectivityVector(indicesOfCenterCones);
        typesOfCenterCones = obj.inputConeMosaic.coneTypes(indicesOfCenterCones);

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

        allCellsCenterConeTypeNum(iRGC,:) = theCenterConeTypeNum;
        allCellsCenterConeTypeWeights(iRGC,:) = theCenterConeTypeWeights;
        allCellsMajorityConeType(iRGC) = theMajorityConeType;
        allCellsCenterConesNum(iRGC) = numel(indicesOfCenterCones);
    end
     
end