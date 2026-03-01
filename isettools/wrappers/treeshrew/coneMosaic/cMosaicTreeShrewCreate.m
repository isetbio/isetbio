function theConeMosaic = cMosaicTreeShrewCreate(varargin)
% Create a @cMosaic object for the TreeShrew retina
%
% Syntax:
%   theCMosaic = CMOSAICTREESHREWCREATE(varargin)
%
% Description:
%   Generate a @cMosaic object with treeshrew params
%
% History:
%    02/28/26  NPC  Rewrote it for isetbio/isetcam.

    p = inputParser;
    p.addParameter('fovDegs', [2 1], @isnumeric);
    p.addParameter('integrationTimeSeconds', 5/1000, @isnumeric);
    % Parse input
    p.parse(varargin{:});


    % Generate a treeshrew @cPigment object for the tree tree shrew
    thePhotopigment = treeShrewCPhotopigment();

    % Generate a treeshrew @Macula object for the three shrew
    theMacularPigment = macularPigmentTreeShrewCreate(thePhotopigment.wave);

    [coneLocationsMicrons, coneApertureDiametersMicrons, coneTypes, micronsPerDegree] = ...
        generateComponents(p.Results.fovDegs);

    % Generate struct with custom cone data
    theConeDataStruct = struct(...
        'positionUnits', 'microns', ...
        'positions', coneLocationsMicrons, ...
        'types', coneTypes,...
        'lightGatheringApertureDiameters', coneApertureDiametersMicrons ...                      
        );

    % Generate a @cMosaic object 
    theConeMosaic = cMosaic(...
        'coneData', theConeDataStruct, ...
        'micronsPerDegree', micronsPerDegree, ...
        'integrationTime', p.Results.integrationTimeSeconds, ...
        'pigment', thePhotopigment ,...
        'macular', theMacularPigment);

end


function [coneLocationsMicrons, coneApertureDiametersMicrons, coneTypes, micronsPerDegree] = ...
    generateComponents(mosaicFovDegs)

    mosaicSizeDegs = max(mosaicFovDegs(:));

    defaultParams = coneMosaicTreeShrewDefaultParams();

    theOI = oiTreeShrewCreate('pupilDiameterMM', 2.0);
    micronsPerDegree = theOI.optics.micronsPerDegree;

    mosaicRadiusMicrons = 0.5*mosaicSizeDegs*micronsPerDegree;

    coneLocationsMicrons = generateLatticeOfConePositions(...
        mosaicRadiusMicrons, ...
        defaultParams.centralRetinaConeSpacingMicrons);
    conesNum = size(coneLocationsMicrons,1);

    coneRFpositionsDegs = coneLocationsMicrons/micronsPerDegree;
    coneRFspacingsDegs = ones(1, conesNum) * defaultParams.centralRetinaConeSpacingMicrons/micronsPerDegree;
    coneApertureDiametersMicrons = ones(1, conesNum)*defaultParams.innerSegmentDiameterMicrons;

    coneTypes = assignConeTypes(defaultParams.coneDensities, coneRFpositionsDegs, coneRFspacingsDegs);


    idx = find(...
        abs(coneRFpositionsDegs(:,1)) <= 0.5*mosaicFovDegs(1)  & ...
        abs(coneRFpositionsDegs(:,2)) <= 0.5*mosaicFovDegs(2));

    coneTypes = coneTypes(idx);
    coneLocationsMicrons = coneLocationsMicrons(idx,:); 
    coneApertureDiametersMicrons = coneApertureDiametersMicrons(idx)
end


function conePositionsMicrons = generateLatticeOfConePositions(mosaicRadiusMicrons, lambdaMicrons)

    % Generate a mosaic that is 20% larger to minimize edge effects
    margin = 1.2;
    radius = margin*mosaicRadiusMicrons;
    rows = ceil(2 * radius);
    cols = rows;
    
    conePositionsMicrons = generateRegularHexGrid(rows, cols, lambdaMicrons);
end


function hexLocs = generateRegularHexGrid(rows, cols, lambda)
    
    scaleF = sqrt(3) / 2;
    extraCols = round(cols / scaleF) - cols;
    rectXaxis2 = (1:(cols + extraCols));
    [X2, Y2] = meshgrid(rectXaxis2, 1:rows);

    X2 = X2 * scaleF ;
    for iCol = 1:size(Y2, 2)
        Y2(:, iCol) = Y2(:, iCol) - mod(iCol - 1, 2) * 0.5;
    end

    % Scale for lambda
    X2 = X2 * lambda;
    Y2 = Y2 * lambda;
    marginInRFPositions = 0.1;
    indicesToKeep = (X2 >= -marginInRFPositions) & ...
                    (X2 <= cols+marginInRFPositions) &...
                    (Y2 >= -marginInRFPositions) & ...
                    (Y2 <= rows+marginInRFPositions);
    xHex = X2(indicesToKeep);
    yHex = Y2(indicesToKeep);
    
    % Center central cone at (0,0)
    hexLocs = [xHex(:) - mean(xHex(:)) yHex(:) - mean(yHex(:))];
    mins = min(abs(hexLocs),[],1);
    hexLocs = bsxfun(@minus, hexLocs, [mins(1) 0]);
end


function coneTypes = assignConeTypes(coneDensities, coneRFpositionsDegs, coneRFspacingsDegs)


    coneDensities = coneDensities / sum(coneDensities);
    
    conesNum = size(coneRFpositionsDegs,1);
    
    desiredSconesNum = round(coneDensities(cMosaic.SCONE_ID)*conesNum);

    allConeIndices = 1:conesNum;

    if (coneDensities(cMosaic.SCONE_ID) > 0) && (coneDensities(cMosaic.SCONE_ID) <= 0.5)
        % Determine non-foveal LM cone indices, leaving room for regularly
        % spaced S-cones with density = coneDensities(obj.SCONE_ID)
        idx = determineLMconeIndices(...
            coneRFpositionsDegs, ...
            coneRFspacingsDegs, ...
            coneDensities(cMosaic.SCONE_ID), desiredSconesNum);
        allLMconeIndices = allConeIndices(idx);
    else
        % S-cone density either 0 or too high, so no regularity
        idx = randperm(numel(allConeIndices));
        if (desiredSconesNum >= numel(allConeIndices))
            LMconeIndices = [];
        else
            LMconesNum = numel(allConeIndices) - desiredSconesNum;
            LMconeIndices = idx(1:LMconesNum);
        end
        allLMconeIndices = allConeIndices(LMconeIndices);
    end
    
   
    % Assign L and M-cone indices based on relative ratio of L:M cone density
    LMratio = coneDensities(cMosaic.LCONE_ID) / coneDensities(cMosaic.MCONE_ID);
    
    p = rand(1,numel(allLMconeIndices));
    if (isinf(LMratio))
        idx = 1:numel(allLMconeIndices);
    elseif LMratio == 0
        idx = [];
    else
        idx = find(p<LMratio/(1+LMratio));
    end
    
    lConeIndices = allLMconeIndices(idx);
    mConeIndices = setdiff(allLMconeIndices, lConeIndices);
    sConeIndices = setdiff(1:conesNum, allLMconeIndices);
    kConeIndices = [];
   
    % Make sure all cones have been assigned an ID
    assert(conesNum==numel(lConeIndices)+numel(mConeIndices)+numel(sConeIndices)+numel(kConeIndices), ...
        'Indices do not sum up to total cones');
    
    % Assign cone types 
    coneTypes = zeros(conesNum,1);
    coneTypes(lConeIndices) = cMosaic.LCONE_ID;
    coneTypes(mConeIndices) = cMosaic.MCONE_ID;
    coneTypes(sConeIndices) = cMosaic.SCONE_ID;
    coneTypes(kConeIndices) = cMosaic.KCONE_ID;
    
    % Reshape indices
    lConeIndices = reshape(lConeIndices, [numel(lConeIndices) 1]);
    mConeIndices = reshape(mConeIndices, [numel(mConeIndices) 1]);
    sConeIndices = reshape(sConeIndices, [numel(sConeIndices) 1]);
    kConeIndices = reshape(kConeIndices, [numel(kConeIndices) 1]);
    
    % Compute achieved cone densities
    achievedConeDensities = [...
        numel(lConeIndices)/conesNum ...
        numel(mConeIndices)/conesNum ...
        numel(sConeIndices)/conesNum ...
        numel(kConeIndices)/conesNum];
    
   
end


function LMconeIndices = determineLMconeIndices(conePositions, coneSpacings, desiredSconeDensity, desiredSconesNum)

    % All cones num
    conesNum = size(conePositions,1);
    
    % Determine relative S-cone spacing from S-cone density
    testSconeSpacing = 1.1: 0.01 :10;
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
