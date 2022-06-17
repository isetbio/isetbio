function updateCentroidsFromInputs(obj, RGClist)
% Update the centroids of all RGCs in the RGClist. If an RGC has 0 inputs
% its centroid is set to Inf
        
    cm = obj.coneConnectivityMatrix;
    allConePosMicrons = obj.inputConeMosaic.coneRFpositionsMicrons;
    
    centroids = inf(numel(RGClist),2);

    parfor iRGC = 1:numel(RGClist)
        theRGCindex = RGClist(iRGC);
        connectedConeIndices = find(squeeze(cm(:, theRGCindex))>0);

        if (~isempty(connectedConeIndices))
            % Weights of these input cones
            inputConeWeights = full(cm(connectedConeIndices, theRGCindex));
    
            % Positions and spacings of these input cones
            inputConePositions = allConePosMicrons(connectedConeIndices,:);
        
            % Compute centroids
            % [~, centroids(iRGC,:)] = var(inputConePositions,inputConeWeights,1);
           
            p1 = any(isinf(inputConePositions(:)));
            p2 = any(isinf(inputConeWeights));
            if (p1 || p2)
                error('Oh oh');
            end
            
            if (sum(inputConeWeights(:)) == 0)
                error('inputConeWeights = 0')
            end
            
            centroids(iRGC,:) = RGCconnector.weightedMean(inputConePositions, inputConeWeights);
        end
    end

    obj.RGCRFcentroidsFromInputs(RGClist,:) = centroids;
end