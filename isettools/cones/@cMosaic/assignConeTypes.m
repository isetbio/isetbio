function assignConeTypes(obj)

    % Ensure densities sum up to 1.0
    if (numel(obj.coneDensities) == 3)
        obj.coneDensities = [obj.coneDensities 0];
    end
    remaining = 1-obj.coneDensities(obj.SCONE_ID);
    f = 1/remaining;
    obj.coneDensities(obj.LCONE_ID) = obj.coneDensities(obj.LCONE_ID) * f;
    obj.coneDensities(obj.MCONE_ID) = obj.coneDensities(obj.MCONE_ID) * f;
    obj.coneDensities(obj.KCONE_ID) = obj.coneDensities(obj.KCONE_ID) * f;
    
    % Reserve all cones within the tritanopic area to be either L or M
    ecc = sqrt(sum(obj.coneRFpositionsDegs.^2,2));
    fovealLMconeIndices = find(ecc <= 0.5*obj.tritanopicRadiusDegs);
    
    % Determine non-foveal LM cone indices, leaving room for regularly
    % spaced S-cones with density = obj.coneDensities(obj.SCONE_ID)
    conesNum = size(obj.coneRFpositionsDegs,1);
    peripheralConeIndices = setdiff(1:conesNum, fovealLMconeIndices);
    idx = determineLMconeIndices(...
        obj.coneRFpositionsDegs(peripheralConeIndices,:), ...
        obj.coneRFspacingsDegs (peripheralConeIndices), ...
        obj.coneDensities(obj.SCONE_ID));
    peripheralLMConeIndices = peripheralConeIndices(idx);
    
    % All LM cone indices
    allLMconeIndices = [fovealLMconeIndices(:); peripheralLMConeIndices(:)];
    
    % Assing L and M-cone indices based on relative ratio of L:M cone density
    LMratio = obj.coneDensities(obj.LCONE_ID) / obj.coneDensities(obj.MCONE_ID);
    
    p = rand(1,numel(allLMconeIndices));
    if (isinf(LMratio))
        idx = 1:numel(allLMconeIndices);
    else
        idx = find(p<LMratio/(1+LMratio));
    end
    
    obj.lConeIndices = allLMconeIndices(idx);
    obj.mConeIndices = setdiff(allLMconeIndices, obj.lConeIndices);
    obj.sConeIndices = setdiff(1:conesNum, allLMconeIndices);
    
    % k-cones: randomized positions
    randomIndices = randperm(numel(allLMconeIndices));
    kConesNum = round(obj.coneDensities(obj.KCONE_ID) * conesNum);
    
    % These LM cones will be reassigned to k-cones
    kConeIDs = allLMconeIndices(randomIndices(1:kConesNum));
    
    % Remove them from the list of L-cones
    obj.lConeIndices = setdiff(obj.lConeIndices, kConeIDs);
    % Remove them from the list of M-cones
    obj.mConeIndices = setdiff(obj.mConeIndices, kConeIDs);
    % Add them to the list of K-cones
    obj.kConeIndices = kConeIDs;
    
   
    % Make sure all cones have been assigned an ID
    assert(conesNum==numel(obj.lConeIndices)+numel(obj.mConeIndices)+numel(obj.sConeIndices)+numel(obj.kConeIndices), ...
        'Indices do not sum up to total cones');
    
    % Assign cone types 
    obj.coneTypes = zeros(conesNum,1);
    obj.coneTypes(obj.lConeIndices) = obj.LCONE_ID;
    obj.coneTypes(obj.mConeIndices) = obj.MCONE_ID;
    obj.coneTypes(obj.sConeIndices) = obj.SCONE_ID;
    obj.coneTypes(obj.kConeIndices) = obj.KCONE_ID;
    
    % Reshape indices
    obj.lConeIndices = reshape(obj.lConeIndices, [numel(obj.lConeIndices) 1]);
    obj.mConeIndices = reshape(obj.mConeIndices, [numel(obj.mConeIndices) 1]);
    obj.sConeIndices = reshape(obj.sConeIndices, [numel(obj.sConeIndices) 1]);
    obj.kConeIndices = reshape(obj.kConeIndices, [numel(obj.kConeIndices) 1]);
    
    % Update cone densities with achieved ones
    obj.coneDensities = [...
        numel(obj.lConeIndices)/conesNum ...
        numel(obj.mConeIndices)/conesNum ...
        numel(obj.sConeIndices)/conesNum ...
        numel(obj.kConeIndices)/conesNum];
    
    fprintf('Achieved cone densities: L (%2.3f), M (%2.3f), S (%2.3f), K (%2.3f)\n', ...
        obj.coneDensities(1), ...
        obj.coneDensities(2), ...
        obj.coneDensities(3), ...
        obj.coneDensities(4));
    
end


function LMconeIndices = determineLMconeIndices(conePositions, coneSpacings, desiredSconeDensity)

    % Determine relative S-cone spacing from S-cone density
    testSconeSpacing = 1.5:0.1:10;
    sConeDensityFromNormalizedSconeSpacing = 0.01*(1 + 30*exp(-0.8*(testSconeSpacing-1).^0.9));
    [~,idx] = min(abs(desiredSconeDensity-sConeDensityFromNormalizedSconeSpacing));
    relativeSconeSpacing = testSconeSpacing(idx);
    
    % All cones num
    conesNum = size(conePositions,1);
    
    % Compute ecc of all cones
    ecc = sqrt(sum(conePositions.^2,2));
    
    % Compute distances between each cone and its closest 100 cones
    [d, i] = pdist2(conePositions, conePositions, 'euclidean', 'smallest', 200);
    
    % Remove the distance to the cone itself
    d = d(2:end,:);
    i = i(2:end,:);
    
    % Go through all cones assigning as S-cones those that are no closer
    % than coneSpacingsMicrons(coneIndex)*relativeSconeSpacing from each other
    coneIndex = 1;
    LMconeIndices = [];
    SconeIndices = [];
    remainingConeIndices = 1:conesNum;
    
    while (numel(remainingConeIndices)>0)
        % Leave the type of the current coneIndex as S.
        SconeIndices = cat(2, SconeIndices, coneIndex);
        % This means all cones around it within the exclusion radius must be non S
        currentExclusionRadius = coneSpacings(coneIndex)*relativeSconeSpacing;
        distancesToNearbyCones = d(:,coneIndex);
        idx = find(distancesToNearbyCones < currentExclusionRadius);
        LMconeIndices = cat(1, LMconeIndices, squeeze(i(idx, coneIndex)));
        % Keep a count of the remaining cone indices that need to be visited
        remainingConeIndices = setdiff(remainingConeIndices, SconeIndices);
        remainingConeIndices = setdiff(remainingConeIndices, LMconeIndices);
        % Next cone to visit
        [~,idx] = min(ecc(remainingConeIndices));
        coneIndex = remainingConeIndices(idx);
    end
    LMconeIndices = unique(LMconeIndices);
end