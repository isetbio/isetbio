function divergeSourceRFsToNearbyDestinationRFs(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('destinationRFoverlapRatio', [], @(x)((isempty(x))||(isscalar(x)&&(x>=0)&&(x<1))));
    p.parse(varargin{:});

    if (obj.connectivityMatrixIsNonExclusiveAnyMore)
        fprintf(2, 'The connectivityMatrix is non-exclusive anymore. Will not diverge inputs again.')
        return;
    end
    
    if (~isempty(p.Results.destinationRFoverlapRatio)) || (p.Results.destinationRFoverlapRatio == 0)
        obj.wiringParams.destinationRFoverlapRatio = p.Results.destinationRFoverlapRatio;
    end

    destinationRFsNum = size(obj.connectivityMatrix,2);
    sourceRFsNum = size(obj.connectivityMatrix,1);

    if (obj.wiringParams.destinationRFoverlapRatio == 0.0)
        % Find indices of sourceRFs that are connected to destinationRFs
        connectedSourceRFIndices = find(sum(abs(obj.connectivityMatrix),2)>0);

        fprintf('There are %d source RFs (%d of which are connected to %d RGCs with ZERO RF OVERLAP)\n', sourceRFsNum , numel(connectedSourceRFIndices), destinationRFsNum);
        totalInputForDestinationRF = sum(abs(obj.connectivityMatrix(connectedSourceRFIndices,:)),1);
        totalOutputForEachConnectedSourceRF = sum(abs(obj.connectivityMatrix(connectedSourceRFIndices,:)),2);
    
        fprintf('max input across all %d destination RFs: %f\n', numel(totalInputForDestinationRF ), max(full(totalInputForDestinationRF )));
        fprintf('min input across all %d destination RFs: %f\n', numel(totalInputForDestinationRF), min(full(totalInputForDestinationRF)));
        fprintf('max output across all %d connected sourceRFs: %f\n', numel(totalOutputForEachConnectedSourceRF), max(full(totalOutputForEachConnectedSourceRF)));
        fprintf('min output across all %d connected sourceRFs: %f\n', numel(totalOutputForEachConnectedSourceRF), min(full(totalOutputForEachConnectedSourceRF)));
        return;
    end

    % Non-zero overlap, so set the connectivityMatrixIsNonExclusiveAnyMore
    obj.connectivityMatrixIsNonExclusiveAnyMore = true;
    
    % Find the nearby destination RFs for each destinationRF
    [~, nearbyDestinationRFIndices] = ...
        MosaicConnector.pdist2(...
            obj.destinationRFcentroidsFromInputs, ...
            obj.destinationRFcentroidsFromInputs, ...
            'smallest', obj.wiringParams.maxNeighborsNum+1);
    nearbyDestinationRFIndices = nearbyDestinationRFIndices(2:end,:);

    % Find the mean number (locally averaged across all neighboring destination RFs) of 
    % exclusive source RF inputs for each destination RF
    localMeanSourceRFinputsNum = zeros(1,destinationRFsNum);
    for iDestinationRF = 1:destinationRFsNum
        nearbyDestinationRFindicesForThisDestinationRF = nearbyDestinationRFIndices(:,iDestinationRF);
        exclusiveSourceRFnum = zeros(1,numel(nearbyDestinationRFindicesForThisDestinationRF));
        for ii = 1:numel(nearbyDestinationRFindicesForThisDestinationRF)
            exclusiveSourceRFnum(ii) = numel(find(squeeze(obj.connectivityMatrix(:, nearbyDestinationRFindicesForThisDestinationRF(ii)))>0));
        end
        localMeanSourceRFinputsNum(iDestinationRF) = mean(exclusiveSourceRFnum);
    end


    allSourceRFPositions = obj.sourceLattice.RFpositionsMicrons;
    for iDestinationRF = 1:destinationRFsNum 

        theDestinationRFspacing = obj.destinationRFspacingsFromCentroids(iDestinationRF);
        overlapRadius = MosaicConnector.radiusToAchieveOverlap(obj.wiringParams.destinationRFoverlapRatio, theDestinationRFspacing);
        overlapSigma = overlapRadius / 3.0;

        % Find sourceRF indices within overlapRadius distance from centroid
        theDestinationRFcentroid = obj.destinationRFcentroidsFromInputs(iDestinationRF,:);
        distances = MosaicConnector.pdist2(allSourceRFPositions(obj.connectableSourceRFindices,:), theDestinationRFcentroid);
        idx = find(distances <= overlapRadius);
        sourceRFIndicesWithinOverlapRadius = obj.connectableSourceRFindices(idx);

        % Indices of the exclusice source RFs
        exclusiveSourceRFIndices = find(squeeze(obj.connectivityMatrix(:, iDestinationRF))>0);
        
        % SourceRF indices within overlapRadius distance from centroid that are non-exclusive
        nonExclusiveSourceRFIndicesWithinOverlapRadius = setdiff(sourceRFIndicesWithinOverlapRadius,exclusiveSourceRFIndices);
        
        % Find the index of the non-exclusive sourceRF that is closest to the centroid
        [~, idx] = min(MosaicConnector.pdist2(allSourceRFPositions(nonExclusiveSourceRFIndicesWithinOverlapRadius,:), theDestinationRFcentroid));
        if (isempty(idx))
            continue;
        end
        closestNonExclusiveSourceRFindex = nonExclusiveSourceRFIndicesWithinOverlapRadius(idx);

        % Modify the centroid to move it toward the closestNonExclusiveSourceRFindex
        oldCentroidBias = 0.9;
        if (numel(exclusiveSourceRFIndices) == 1)
            oldCentroidBias = 0.9;
        end
        theUpdatedDestinationRFcentroid = oldCentroidBias     * theDestinationRFcentroid + ...
                                          (1-oldCentroidBias) * allSourceRFPositions(closestNonExclusiveSourceRFindex,:);

        
        % Distances of sourceRFs within the overlap radius from the updatedDestinationRFcentroid
        distances = MosaicConnector.pdist2(allSourceRFPositions(sourceRFIndicesWithinOverlapRadius,:), theUpdatedDestinationRFcentroid);

        % Updated source RF indices and weights
        weights = exp(-0.5*(distances/overlapSigma).^2);

        newlyConnectedSourceRFindices = sourceRFIndicesWithinOverlapRadius;
        newlyConnectedSourceRFweights = weights;

        % Only keep inputs with weight > maxWeight*threshold
        threshold = 1/100;
        idx = find(newlyConnectedSourceRFweights > threshold * max(newlyConnectedSourceRFweights(:)));
        newlyConnectedSourceRFindices = newlyConnectedSourceRFindices(idx);
        newlyConnectedSourceRFweights = newlyConnectedSourceRFweights(idx);


        % Update connectivity matrix
        obj.connectivityMatrix(:, iDestinationRF) = 0;
        obj.connectivityMatrix(newlyConnectedSourceRFindices, iDestinationRF) = newlyConnectedSourceRFweights;
        
        debugOverlap = false;
        if (debugOverlap)

            iidx = find(obj.connectivityMatrix(:, iDestinationRF)>0);
            [newRF, spatialSupportNewX, spatialSupportNewY] = ...
                generateRFoutline(allSourceRFPositions(iidx,:), full(obj.connectivityMatrix(iidx,iDestinationRF)), [],[]);
            oldRF = ...
                generateRFoutline(allSourceRFPositions(exclusiveSourceRFIndices,:), ones(1:numel(exclusiveSourceRFIndices)), spatialSupportNewX, spatialSupportNewY);

            figure(5555); clf;
            subplot(1,2,1);
            max(oldRF(:))
            imagesc(spatialSupportNewX, spatialSupportNewY, oldRF, [0 1]);
            hold on;
            plot(theDestinationRFcentroid(1), theDestinationRFcentroid(2), 'rx');
            plot([spatialSupportNewX(1) spatialSupportNewX(end)], mean(spatialSupportNewY)*[1 1], 'r-');
            plot(mean(spatialSupportNewX)*[1 1], [spatialSupportNewY(1) spatialSupportNewY(end)], 'r-');
    
            axis 'image'; axis 'xy'
            set(gca, 'XLim', [spatialSupportNewX(1) spatialSupportNewX(end)], ...
                     'YLim', [spatialSupportNewY(1) spatialSupportNewY(end)]);
    
            subplot(1,2,2);
            max(newRF(:))
            imagesc(spatialSupportNewX, spatialSupportNewY, newRF, [0 1]); hold('on');
            plot(theUpdatedDestinationRFcentroid(1), theUpdatedDestinationRFcentroid(2), 'rx');
            plot([spatialSupportNewX(1) spatialSupportNewX(end)], mean(spatialSupportNewY)*[1 1], 'r-');
            plot(mean(spatialSupportNewX)*[1 1], [spatialSupportNewY(1) spatialSupportNewY(end)], 'r-');
    
            axis 'image'; axis 'xy'
            set(gca, 'XLim', [spatialSupportNewX(1) spatialSupportNewX(end)], ...
                     'YLim', [spatialSupportNewY(1) spatialSupportNewY(end)]);
            colormap(gray(1024));
            pause
        end

    end


    % Find indices of sourceRFs that are connected to destinationRFs
    connectedSourceRFIndices = find(sum(obj.connectivityMatrix,2)>0);

    % Sum of the output of each sourceRF to all destinationRFs must equal 1.0;
    totalOutputForEachConnectedSourceRF = sum(obj.connectivityMatrix(connectedSourceRFIndices,:),2);
    obj.connectivityMatrix(connectedSourceRFIndices,:) = ...
        bsxfun(@times, ...
        obj.connectivityMatrix(connectedSourceRFIndices,:), ...
        1./totalOutputForEachConnectedSourceRF);

    % Update the destination RF centroids
    allDestinationRFindices = 1:destinationRFsNum;
    obj.updateDestinationCentroidsFromInputs(allDestinationRFindices);

    % Update the destination RF spacings
    obj.updateDestinationRFspacingsBasedOnCentroids();

    fprintf('There are %d source RFs (%d of which are connected to %d RGCs with %1.2f RF OVERLAP)\n', ...
        sourceRFsNum , numel(connectedSourceRFIndices), destinationRFsNum, obj.wiringParams.destinationRFoverlapRatio);
    totalInputForDestinationRF = sum(abs(obj.connectivityMatrix(connectedSourceRFIndices,:)),1);
    totalOutputForEachConnectedSourceRF = sum(abs(obj.connectivityMatrix(connectedSourceRFIndices,:)),2);

    fprintf('max input across all %d destination RFs: %f\n', numel(totalInputForDestinationRF ), max(full(totalInputForDestinationRF )));
    fprintf('min input across all %d destination RFs: %f\n', numel(totalInputForDestinationRF), min(full(totalInputForDestinationRF)));
    fprintf('max output across all %d connected sourceRFs: %f\n', numel(totalOutputForEachConnectedSourceRF), max(full(totalOutputForEachConnectedSourceRF)));
    fprintf('min output across all %d connected sourceRFs: %f\n', numel(totalOutputForEachConnectedSourceRF), min(full(totalOutputForEachConnectedSourceRF)));

end

function [RF, spatialSupportX, spatialSupportY] = generateRFoutline(positions, weights, spatialSupportX, spatialSupportY)

    if (isempty(spatialSupportX))
        minX = floor(min(positions(:,1)))-1;
        maxX = ceil(max(positions(:,1)))+1;
        minY = floor(min(positions(:,2)))-1;
        maxY = ceil(max(positions(:,2)))+1;
    
        rangeX = maxX-minX;
        rangeY = maxY-minY;
        if (rangeX > rangeY)
            range = rangeX;
        else
            range = rangeY;
        end

        xo = (minX+maxX)/2;
        yo = (minY+maxY)/2;
    
        spatialSupportX = linspace(xo-range/2,xo+range/2,16);
        spatialSupportY = linspace(yo-range/2,yo+range/2,16);
        
    end

    [X,Y] = meshgrid(spatialSupportX, spatialSupportY);
    RF = X * 0;

    for i = 1:size(positions,1)
        [~,ix] = min(abs(positions(i,1)-spatialSupportX));
        [~,iy] = min(abs(positions(i,2)-spatialSupportY));
        RF(iy,ix) = RF(iy,ix) + weights(i);
    end

end
