function divergeSourceRFsToNearbyDestinationRFs(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('destinationRFoverlapRatio', [], @(x)((isempty(x))||(isscalar(x)&&(x>=0))));
    p.parse(varargin{:});

    % Remove any negative weights (which indicate overlapping cone weights)
    obj.connectivityMatrix(obj.connectivityMatrix<0) = 0;

    if (~isempty(p.Results.destinationRFoverlapRatio))
        obj.wiringParams.destinationRFoverlapRatio = p.Results.destinationRFoverlapRatio;
    end

    destinationRFsNum = size(obj.connectivityMatrix,2);
    sourceRFsNum = size(obj.connectivityMatrix,1);

    if (obj.wiringParams.destinationRFoverlapRatio == 0.0)
        % Find indices of sourceRFs that are connected to destinationRFs
        connectedSourceRFIndices = find(sum(abs(obj.connectivityMatrix),2)>0);

        fprintf('There are %d source RFs (%d of which are connected to %d RGCs\n', sourceRFsNum , numel(connectedSourceRFIndices), destinationRFsNum);
        totalInputForDestinationRF = sum(abs(obj.connectivityMatrix(connectedSourceRFIndices,:)),1);
        totalOutputForEachConnectedSourceRF = sum(abs(obj.connectivityMatrix(connectedSourceRFIndices,:)),2);
    
        fprintf('max input across all %d destination RFs: %f\n', numel(totalInputForDestinationRF ), max(full(totalInputForDestinationRF )));
        fprintf('min input across all %d destination RFs: %f\n', numel(totalInputForDestinationRF), min(full(totalInputForDestinationRF)));
        fprintf('max output across all %d connected sourceRFs: %f\n', numel(totalOutputForEachConnectedSourceRF), max(full(totalOutputForEachConnectedSourceRF)));
        fprintf('min output across all %d connected sourceRFs: %f\n', numel(totalOutputForEachConnectedSourceRF), min(full(totalOutputForEachConnectedSourceRF)));
        return;
    end

    
    allSourceRFPositions = obj.sourceLattice.RFpositionsMicrons;
    theFullConnectivityMatrix = cell(destinationRFsNum,2);

    for iDestinationRF = 1:destinationRFsNum
        theDestinationRFcentroid = obj.destinationRFcentroidsFromInputs(iDestinationRF,:);
        theDestinationRFspacing = obj.destinationRFspacingsFromCentroids(iDestinationRF);
        overlapSigma = obj.wiringParams.destinationRFoverlapRatio * theDestinationRFspacing;
        overlapRadius = overlapSigma*3.0;

        % Indices of non-overlapping source RFs
        nonOverlappingSourceRFIndices = find(squeeze(obj.connectivityMatrix(:, iDestinationRF))>0);
        % Weights of non-overlapping source RFs
        nonOverlappingSourceRFWeights = full(obj.connectivityMatrix(nonOverlappingSourceRFIndices, iDestinationRF ));

        % Find sourceRF indices within overlapRadius distance from centroid
        distances = MosaicConnector.pdist2(allSourceRFPositions(obj.connectableSourceRFindices,:), theDestinationRFcentroid);

        % Find the indices of the overlapping cones
        idx = find(distances <= overlapRadius);
        sourceRFIndicesWithinOverlapRadius = obj.connectableSourceRFindices(idx);
        sourceRFDistancesWithinOverlapRadius = distances(idx);
        sourceRFWeightsWithinOverlapRadius = exp(-0.5*(sourceRFDistancesWithinOverlapRadius/overlapSigma).^2);

        % Find which of the sourceRFIndicesWithinOverlapRadius are the main (non-overlapping cones)
        [isMainCone, ia] = ismember(sourceRFIndicesWithinOverlapRadius, nonOverlappingSourceRFIndices);

        % Find the max Gaussian weight that would have been assigned to the non-overlapping  cones
        [overlappingConesFound,ib] = ismember(nonOverlappingSourceRFIndices, sourceRFIndicesWithinOverlapRadius);

        if (~any(overlappingConesFound))
            % No cones exist that are overlapping and non main ones
            theFullConnectivityMatrix{iDestinationRF}{1} = nonOverlappingSourceRFIndices;
            theFullConnectivityMatrix{iDestinationRF}{2} = nonOverlappingSourceRFWeights;
            continue;
        end

        scalingFactor = max(sourceRFWeightsWithinOverlapRadius(ib(ib>0)));

        for iCone = 1:numel(sourceRFIndicesWithinOverlapRadius)
            if (isMainCone(iCone))
                % Non-overlapping cone, so keep original weight
                sourceRFWeightsWithinOverlapRadius(iCone) = nonOverlappingSourceRFWeights(ia(iCone));
            else
                % Overlapping cone. Set the weights to negative polarity so as to
                % discriminate between main and overlapping cones
                w = min([1 sourceRFWeightsWithinOverlapRadius(iCone)/scalingFactor]);
                if (w > 1) || ( w < 0)
                    error('How can this be?')
                end
                sourceRFWeightsWithinOverlapRadius(iCone) = -w;
            end
        end

        % Modify connectivity matrix to encode the weights of the overlaping cone inputs. 
        theFullConnectivityMatrix{iDestinationRF}{1} = sourceRFIndicesWithinOverlapRadius;
        theFullConnectivityMatrix{iDestinationRF}{2} = sourceRFWeightsWithinOverlapRadius;
    end

    % Finalize connectivity matrix
    for iDestinationRF = 1:destinationRFsNum
        sourceRFIndices = theFullConnectivityMatrix{iDestinationRF}{1};
        sourceRFWeights = theFullConnectivityMatrix{iDestinationRF}{2};
        if (sum(sourceRFWeights (:)) == 0)
            error('How can this be?')
        end
        obj.connectivityMatrix(sourceRFIndices, iDestinationRF) = sourceRFWeights;
    end

    % Find indices of sourceRFs that are connected to destinationRFs
    connectedSourceRFIndices = find(sum(abs(obj.connectivityMatrix),2)>0);
    
    % Sum of the output of each sourceRF to all destinationRFs must equal 1.0;
    totalOutputForEachConnectedSourceRF = sum(abs(obj.connectivityMatrix(connectedSourceRFIndices,:)),2);
    obj.connectivityMatrix(connectedSourceRFIndices,:) = ...
        bsxfun(@times, ...
        obj.connectivityMatrix(connectedSourceRFIndices,:), ...
        1./totalOutputForEachConnectedSourceRF);

    fprintf('There are %d source RFs (%d of which are connected to %d RGCs\n', sourceRFsNum , numel(connectedSourceRFIndices), destinationRFsNum);
    totalInputForDestinationRF = sum(abs(obj.connectivityMatrix(connectedSourceRFIndices,:)),1);
    totalOutputForEachConnectedSourceRF = sum(abs(obj.connectivityMatrix(connectedSourceRFIndices,:)),2);

    fprintf('max input across all %d destination RFs: %f\n', numel(totalInputForDestinationRF ), max(full(totalInputForDestinationRF )));
    fprintf('min input across all %d destination RFs: %f\n', numel(totalInputForDestinationRF), min(full(totalInputForDestinationRF)));
    fprintf('max output across all %d connected sourceRFs: %f\n', numel(totalOutputForEachConnectedSourceRF), max(full(totalOutputForEachConnectedSourceRF)));
    fprintf('min output across all %d connected sourceRFs: %f\n', numel(totalOutputForEachConnectedSourceRF), min(full(totalOutputForEachConnectedSourceRF)));

end