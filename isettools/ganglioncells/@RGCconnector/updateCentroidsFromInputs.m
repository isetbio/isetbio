function updateCentroidsFromInputs(obj, RGClist)
% Update the centroids of all RGCs in the RGClist
        
    cm = obj.coneConnectivityMatrix;
    allConePosMicrons = obj.inputConeMosaic.coneRFpositionsMicrons;
    
    centroids = zeros(numel(RGClist),2);

    parfor iRGC = 1:numel(RGClist)
        theRGCindex = RGClist(iRGC);
        connectedConeIndices = find(squeeze(cm(:, theRGCindex))>0);

        if (isempty(connectedConeIndices))
            centroids(iRGC,:) = inf(1,2);
        else

            % Weights of these input cones
            inputConeWeights = full(cm(connectedConeIndices, theRGCindex));
    
            % Positions and spacings of these input cones
            inputConePositions = allConePosMicrons(connectedConeIndices,:);
        
            % Compute centroids
            [~, centroids(iRGC,:)] = var(inputConePositions,inputConeWeights,1);
        end
    end

    obj.RGCRFcentroidsFromInputs(RGClist,:) = centroids;
end