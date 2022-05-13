function  [RGCRFinputs, RGCRFweights, availableZeroInputRGCsNum] = ...
    connectEachConeToItsNearestRGC(coneRFpos, coneRFspacings, coneTypes, ...
                                   RGCRFpos, chromaticSpatialVarianceTradeoff, ...
                                   maxNearbyRGCsNum)
% Connect each of a set of input cones (L- or M-) to a set of RGC RF centers
%
% Syntax:
%   [RGCRFinputs, RGCRFweights, availableZeroInputRGCsNum] = ...
%        RGCRFconnector.connectEachConeToItsNearestRGC(...
%                    coneRFpos, coneRFspacings, coneTypes, ...
%                    RGCRFpos, chromaticSpatialVarianceTradeoff)
%
% Description:
%   Connect each of a set of input cones (L- or M-) to a set of RGC RF centers
%
% Inputs:
%    coneRFpos          - [N x 2] matrix of (x,y) positions of N input cones
%    coneRFspacings     - [N x 1] vector of spacings of N input cones
%    coneTypes          - [N x 1] vector of types of N input cones
%    RGCRFpos           - [M x 2] matrix of (x,y) positions of M target RGC RF centers
%    chromaticSpatialVarianceTradeoff  - Chromatic-SpatialVariance tradefoff
%    maxNearbyRGCsNum   - Scalar. Max number of nearby RGCs to look for assignment, if
%                                 a cone cannot be assigned to its closest
%                                 RGC because it already has an input cone
%
% Outputs:
%    RGCRFinputs                - Cell array with indices of the input cones, one cell per each target RGC RF center
%    RGCRFweights               - Cell array with weights of the input cones, one cell per each target RGC RF center   
%    availableZeroInputRGCsNum  - Number of the M target RGCs with zero inputs 
%
% Optional key/value pairs
%   none
%   
% History:
%   5/11/2022       NPC     Wrote it
%

    % Sort cones according to distance from the center of the cone mosaic
    mosaicCenter = mean(coneRFpos,1);
    d = sum((bsxfun(@minus, coneRFpos, mosaicCenter)).^2, 2);
    [~,sortedConeIndices] = sort(d);

    % Find the closest RGC to each cone
    [~,closestMRGCIndices] = pdist2(RGCRFpos, coneRFpos, '', 'smallest', 1);

    % Store the cone inputs in a temporary array
    tmpRGCRFinputs = cell(1,size(RGCRFpos,1));
    
    % Pass 1. Only connect a cone to its nearest RGC iff that RGC has 0 inputs

    % Keep a track of cone indices that were not connected because the
    % nearest RGC was already connected to another cone
    unconnectedConeIndices = [];
    for iCone = 1:size(coneRFpos,1)

        coneIndex = sortedConeIndices(iCone);
        if (coneTypes(coneIndex) == cMosaic.SCONE_ID)
            % Do not connect S-cones
            continue;
        end
        theRGCindex = closestMRGCIndices(coneIndex);
        inputConeIndices = tmpRGCRFinputs{theRGCindex};
        if (numel(inputConeIndices) > 0)
            unconnectedConeIndices(numel(unconnectedConeIndices)+1) = coneIndex;
            continue;
        end
        inputConeIndices(numel(inputConeIndices)+1) = coneIndex;
        tmpRGCRFinputs{theRGCindex} = inputConeIndices;
    end

    fprintf('%d out of %d L/M cones could not be connected to their closest RGC because the closest RGC already had one L/M cone attached\n', ...
        numel(unconnectedConeIndices), numel((coneTypes ~= cMosaic.SCONE_ID)));
    fprintf('Trying to connect them to nearby RGCs, minimizing the cost\n');
    
    % Pass 2. Deal with unconnected cone indices
    d = sum((bsxfun(@minus, coneRFpos(unconnectedConeIndices,:), mosaicCenter)).^2, 2);
    [~, idx] = sort(d, 'ascend');
    unconnectedConeIndices = unconnectedConeIndices(idx);

    for iCone = 1:numel(unconnectedConeIndices)
        theConeIndex = unconnectedConeIndices(iCone);
        theConeRFPosition = coneRFpos(theConeIndex,:);

        % Find neigboring RGCs 
        [~, nearestRGCindices] = pdist2(RGCRFpos, theConeRFPosition, '', 'smallest', maxNearbyRGCsNum);

        % Compute the cost of assigning this cone to each of the nearestRGCs
        projectedCost = zeros(1,numel(nearestRGCindices));

        for iRGC = 1:numel(nearestRGCindices)
            % Retrieve this cell's cone inputs
            theRGCindex = nearestRGCindices(iRGC);
            inputConeIndices = tmpRGCRFinputs{theRGCindex};

            % Add this cone as an input
            inputConeIndices(numel(inputConeIndices)+1) = theConeIndex;
            inputConePositions = coneRFpos(inputConeIndices,:);
            inputConeSpacings = coneRFspacings(inputConeIndices);
            inputConeTypes = coneTypes(inputConeIndices);
            inputConeWeights = ones(1, numel(inputConeIndices));
            % Compute cost, had this cone been assigned to this RGC
            projectedCost(iRGC) = RGCRFconnector.costToMaintainInputCones(chromaticSpatialVarianceTradeoff, ...
                inputConePositions, inputConeSpacings, ...
                inputConeTypes, inputConeWeights);
        end % iRGC

        % Find the RGC with the min projected cost
        [~,idx] = min(projectedCost);
        targetRGCindex = nearestRGCindices(idx);
        
        % Assign that previously non-assigned cone to this RGC
        inputConeIndices = tmpRGCRFinputs{targetRGCindex};
        inputConeIndices(numel(inputConeIndices)+1) = theConeIndex;
        tmpRGCRFinputs{targetRGCindex} = inputConeIndices;
    end

    % Find RGCs that end up having 0 cone inputs and remove them
    nonZeroInputRGCs = [];
    for iRGC = 1:numel(tmpRGCRFinputs)
        numberOfInputs = numel(tmpRGCRFinputs{iRGC});
        if (numberOfInputs > 0)
            nonZeroInputRGCs(numel(nonZeroInputRGCs)+1) = iRGC;
        end
    end
    availableZeroInputRGCsNum = numel(tmpRGCRFinputs)-numel(nonZeroInputRGCs);


    % Only keep the non-zero cone input RGCs
    RGCRFinputsTmp = cell(1,numel(nonZeroInputRGCs));
    RGCRFweightsTmp = cell(1, numel(nonZeroInputRGCs));
    
    centroids = zeros(numel(nonZeroInputRGCs),2);
    for kkk = 1:numel(nonZeroInputRGCs)
        RGCRFinputsTmp{kkk} = tmpRGCRFinputs{nonZeroInputRGCs(kkk)};
        RGCRFweightsTmp{kkk} = RGCRFinputsTmp{kkk}*0+1;
        centroids(kkk,:) = RGCRFconnector.centroidsFromConeInputs({RGCRFinputsTmp{kkk}}, {RGCRFweightsTmp{kkk}}, coneRFpos);
    end

    % Sort RFs according to distance from the center of the RGCmosaic
    mosaicCenter = mean(centroids,1);
    d = sum((bsxfun(@minus, centroids, mosaicCenter)).^2, 2);
    [~,idx] = sort(d);

    RGCRFinputs = cell(1,numel(nonZeroInputRGCs));
    RGCRFweights = cell(1, numel(nonZeroInputRGCs));
    for kkk = 1:numel(idx)
        RGCRFinputs{kkk} = RGCRFinputsTmp{idx(kkk)};
        RGCRFweights{kkk} = RGCRFweightsTmp{idx(kkk)};
    end

end