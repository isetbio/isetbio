function assignConeTypes(obj, src, ~)
% Needs comments from NC.  Not sure why 'src' is here.
%
% Called by the cMosaic constructor
% Looks like the cone positions are calculated by
% regenerateConePositions or initializeConePositions
%
% This routine then 'colors' each of the positions with a cone type.
% 
% See also
%   cMosaic.regenerateConePositions, cMosaic.initializeConePositions,
%   cMosaic.removeConesWithinOpticNerveHead
%

    coneDensities = obj.coneDensities;
    
    coneDensities = abs(coneDensities);
    coneDensities = coneDensities / sum(coneDensities);
    
    % Deal with missing 4-th cone density
    if (numel(coneDensities) == 3)
        coneDensities = [coneDensities 0];
    end
    
    % Reserve all cones within the tritanopic area to be either L or M
    ecc = sqrt(sum(obj.coneRFpositionsDegs.^2,2));

    fovealLMconeIndices = find(ecc <= obj.tritanopicRadiusDegs);
    conesNum = numel(ecc);
    
    % See how to assign S-cones outside the tritanopic area
    peripheralConeIndices = setdiff(1:conesNum, fovealLMconeIndices);
    desiredSconesNum = round(coneDensities(obj.SCONE_ID)*numel(peripheralConeIndices));
    
    if (coneDensities(obj.SCONE_ID) > 0) && (coneDensities(obj.SCONE_ID) <= 0.5)
        % Determine non-foveal LM cone indices, leaving room for regularly
        % spaced S-cones with density = coneDensities(obj.SCONE_ID)
        idx = determineLMconeIndices(...
            obj.coneRFpositionsDegs(peripheralConeIndices,:), ...
            obj.coneRFspacingsDegs (peripheralConeIndices), ...
            coneDensities(obj.SCONE_ID), desiredSconesNum);
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
    LMratio = coneDensities(obj.LCONE_ID) / coneDensities(obj.MCONE_ID);
    
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
    
    if (coneDensities(obj.KCONE_ID) > 0)
        % k-cones: randomized positions
        randomIndices = randperm(numel(allLMconeIndices));
        kConesNum = round(coneDensities(obj.KCONE_ID) * conesNum);

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
    obj.achievedConeDensities = achievedConeDensities;
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
    [distances, indices] = pdist2(conePositions, conePositions, 'euclidean', 'smallest', 200);
    
    % Remove the distance to the cone itself
    distances = distances(2:end,:);
    indices = indices(2:end,:);
    
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
        distancesToNearbyCones = distances(:,coneIndex);
        idx = find(distancesToNearbyCones < currentExclusionRadius);
        LMconeIndices = cat(1, LMconeIndices, squeeze(indices(idx, coneIndex)));
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