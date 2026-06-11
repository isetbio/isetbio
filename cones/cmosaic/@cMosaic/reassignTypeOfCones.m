function reassignTypeOfCones(obj, coneIndices, newConeType)

    
    oldConesNum = numel(obj.coneTypes);
    if (isempty(coneIndices))
        coneIndices = 1:oldConesNum;
    end
    coneIndices = coneIndices(:)';

    % Assert the the cone indices are valid
    assert(min(coneIndices)>0, 'coneIndices must be > 0');
    assert(min(coneIndices)<= oldConesNum, 'coneIndices must be <= %d', oldConesNum);
    assert(max(coneIndices)>0, 'coneIndices must be > 0');
    assert(max(coneIndices)<= oldConesNum, 'coneIndices must be <= %d', oldConesNum);
    
    % Assert that the new cone type is valid
    assert(ismember(newConeType, [cMosaic.LCONE_ID, cMosaic.MCONE_ID, cMosaic.SCONE_ID, cMosaic.KCONE_ID]), ...
        'newConeType (%d) is invalid');

    % Update the lmskConeIndices
    oldConeTypes = obj.coneTypes(coneIndices);
    idx = find(oldConeTypes == cMosaic.LCONE_ID);
    if (~isempty(idx))
        obj.lConeIndices = setdiff(obj.lConeIndices, coneIndices(idx));
    end

    idx = find(oldConeTypes == cMosaic.MCONE_ID);
    if (~isempty(idx))
        obj.mConeIndices = setdiff(obj.mConeIndices, coneIndices(idx));
    end

    idx = find(oldConeTypes == cMosaic.SCONE_ID);
    if (~isempty(idx))
        obj.sConeIndices = setdiff(obj.sConeIndices, coneIndices(idx));
    end

    idx = find(oldConeTypes == cMosaic.KCONE_ID);
    if (~isempty(idx))
        obj.kConeIndices = setdiff(obj.kConeIndices, coneIndices(idx));
    end


    switch (newConeType)
        case cMosaic.LCONE_ID
            obj.lConeIndices = cat(1, obj.lConeIndices(:), coneIndices(:));
        case cMosaic.MCONE_ID(:)
            obj.mConeIndices = cat(1, obj.mConeIndices(:), coneIndices(:));
        case cMosaic.SCONE_ID
            obj.sConeIndices = cat(1, obj.sConeIndices(:), coneIndices(:));
        case cMosaic.KCONE_ID
            obj.kConeIndices = cat(1, obj.kConeIndices(:), coneIndices(:));
    end

    % Change the cone types
    obj.coneTypes(coneIndices) = newConeType;

    conesNum = numel(obj.coneTypes);
    assert(conesNum == oldConesNum, 'fatal error in cone type reassignement.')


    % Compute achieved cone densities
    achievedConeDensities = [...
        numel(obj.lConeIndices)/conesNum ...
        numel(obj.mConeIndices)/conesNum ...
        numel(obj.sConeIndices)/conesNum ...
        numel(obj.kConeIndices)/conesNum];
    
    % Update coneDensities
    obj.achievedConeDensities = achievedConeDensities;

end
