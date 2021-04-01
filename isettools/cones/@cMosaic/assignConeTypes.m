function assignConeTypes(obj, varargin)
    p = inputParser;
    p.addParameter('coneDensities', [], @(x)(isempty(x) || ((numel(x) == 3)||(numel(x)==4))));
    p.parse(varargin{:});
    
    % Ensure densities sum up to 1.0
    if (~isempty(p.Results.coneDensities))
        obj.coneDensities = abs(p.Results.coneDensities);
        obj.coneDensities = obj.coneDensities / sum(obj.coneDensities);
    end
    
    % Deal with missing 4-th cone density
    if (numel(obj.coneDensities) == 3)
        obj.coneDensities = [obj.coneDensities 0];
    end
    
    % Reserve all cones within the tritanopic area to be either L or M
    ecc = sqrt(sum(obj.coneRFpositionsDegs.^2,2));
    fovealLMconeIndices = find(ecc <= obj.tritanopicRadiusDegs);
    conesNum = numel(ecc);
    
    % See how to assign S-cones outside the tritanopic area
    peripheralConeIndices = setdiff(1:conesNum, fovealLMconeIndices);
    desiredSconesNum = round(obj.coneDensities(obj.SCONE_ID)*numel(peripheralConeIndices));
    
    if (obj.coneDensities(obj.SCONE_ID) > 0) && (obj.coneDensities(obj.SCONE_ID) <= 0.5)
        % Determine non-foveal LM cone indices, leaving room for regularly
        % spaced S-cones with density = obj.coneDensities(obj.SCONE_ID)
        idx = determineLMconeIndices(...
            obj.coneRFpositionsDegs(peripheralConeIndices,:), ...
            obj.coneRFspacingsDegs (peripheralConeIndices), ...
            obj.coneDensities(obj.SCONE_ID), desiredSconesNum);
        peripheralLMConeIndices = peripheralConeIndices(idx);
    else
        % S-cone density either 0 or too high, so no regularity
        idx = randperm(numel(peripheralConeIndices));
        if (desiredSconesNum >= numel(peripheralConeIndices))
            LMconeIndices = [];
        else
            LMconesNum = numel(peripheralConeIndices) - desiredSconesNum;
            LMconeIndices = idx(1:LMconesNum);
        end
        peripheralLMConeIndices = peripheralConeIndices(LMconeIndices);
    end
    
    % All LM cone indices
    allLMconeIndices = [fovealLMconeIndices(:); peripheralLMConeIndices(:)];
    
    % Assign L and M-cone indices based on relative ratio of L:M cone density
    LMratio = obj.coneDensities(obj.LCONE_ID) / obj.coneDensities(obj.MCONE_ID);
    
    p = rand(1,numel(allLMconeIndices));
    if (isinf(LMratio))
        idx = 1:numel(allLMconeIndices);
    elseif LMratio == 0
        idx = [];
    else
        idx = find(p<LMratio/(1+LMratio));
    end
    
    obj.lConeIndices = allLMconeIndices(idx);
    obj.mConeIndices = setdiff(allLMconeIndices, obj.lConeIndices);
    obj.sConeIndices = setdiff(1:conesNum, allLMconeIndices);
    
    if (obj.coneDensities(obj.KCONE_ID) > 0)
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
    end
   
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
    
    % Compute achieved cone densities
    achievedConeDensities = [...
        numel(obj.lConeIndices)/conesNum ...
        numel(obj.mConeIndices)/conesNum ...
        numel(obj.sConeIndices)/conesNum ...
        numel(obj.kConeIndices)/conesNum];
    
    % Update coneDensities
    obj.coneDensities = achievedConeDensities;
end


function LMconeIndices = determineLMconeIndices(conePositions, coneSpacings, desiredSconeDensity, desiredSconesNum)

    % All cones num
    conesNum = size(conePositions,1);
    
    % Determine relative S-cone spacing from S-cone density
    testSconeSpacing = 1.1:0.01:10;
    sConeDensityFromNormalizedSconeSpacing = 0.01*(1 + 30*exp(-0.8*(testSconeSpacing-1).^0.9));
    [~,idx] = min(abs(desiredSconeDensity-sConeDensityFromNormalizedSconeSpacing));
    relativeSconeSpacing = testSconeSpacing(idx);

    % Compute ecc of all cones
    ecc = sqrt(sum(conePositions.^2,2));
    
    % Compute distances between each cone and its closest 200 cones
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
    
    if (numel(SconeIndices) > desiredSconesNum)
        SconeIndices = SconeIndices(randperm(numel(SconeIndices)));
        SconeIndices = SconeIndices(1:desiredSconesNum);
        LMconeIndices = setdiff(1:conesNum, SconeIndices);
    elseif (numel(SconeIndices) < desiredSconesNum)
        LMconeIndices = LMconeIndices(randperm(numel(LMconeIndices)));
        desiredLMconesNum = conesNum - desiredSconesNum;
        LMconeIndices = LMconeIndices(1:desiredLMconesNum);
        SconeIndices = setdiff(1:conesNum, LMconeIndices);
    end
end