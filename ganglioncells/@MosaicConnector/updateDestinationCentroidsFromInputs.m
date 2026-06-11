function updateDestinationCentroidsFromInputs(obj, destinationRFList)

    cm = obj.connectivityMatrix;
    sourceRFpositions = obj.sourceLattice.RFpositionsMicrons;
    
    % Initialize the centroids
    centroids = inf(numel(destinationRFList),2);  
    
    
    % Update the centroids of all destination RFs in the destinationRFList. 
    parfor iDestinationRF = 1:numel(destinationRFList)
        theDestinationRFindex = destinationRFList(iDestinationRF);
        connectedSourceRFIndices = find(squeeze(cm(:, theDestinationRFindex))>0);

        if (~isempty(connectedSourceRFIndices))
            % Connection weights of these input source RFs
            inputSourceRFweights = full(cm(connectedSourceRFIndices, theDestinationRFindex));
    
            % Positions and spacings of these input source RFs
            inputSourceRFpositions = sourceRFpositions(connectedSourceRFIndices,:);
        
            p1 = any(isinf(inputSourceRFpositions(:)));
            p2 = any(isinf(inputSourceRFweights));
            if (p1 || p2)
                error('Oh oh');
            end
            
            if (sum(inputSourceRFweights(:)) == 0)
                error('inputSourceRFweights = 0')
            end
            
            centroids(iDestinationRF,:) = MosaicConnector.weightedMean(inputSourceRFpositions, inputSourceRFweights);
        else
            error('No connectedSourceRFIndices');
        end
    end

    obj.destinationRFcentroidsFromInputs(destinationRFList,:) = centroids;
end

