function optimizeTransferOfConeInputsToZeroInputRGC(obj,...
             theSourceRGCindex, theSourceRGCinputConeIndices, ...
             theSourceRGCinputConeWeights, theDestinationRGCindex)

    % Optimize how many and which of theSourceRGCinputConeIndices will be transfered to a zero input RGC
    
    % All combinations of half of the current cones
    halfConesNum = round(numel(theSourceRGCinputConeIndices)/2);
    
    C = nchoosek(1:numel(theSourceRGCinputConeIndices),halfConesNum);
    combinationsNum = size(C,1);
    projectedCostsSourceRGC = zeros(combinationsNum,1);

    parfor iSourceConeComboIndex = 1:combinationsNum
        % The indices of the cones to be transfered
        theConeIndices = theSourceRGCinputConeIndices(C(iSourceConeComboIndex,:));

        % Remove the iSourceConeIndex from theSourceRGCinputIndices
        [newSourceRGCinputConeIndices, ia] = setdiff(theSourceRGCinputConeIndices, theConeIndices);
        newSourceRGCinputConeWeights = theSourceRGCinputConeWeights(ia);

        % Compute projected cost for the source RGC
        projectedCostsSourceRGC(iSourceConeComboIndex) = obj.costToMaintainInputs(...
            newSourceRGCinputConeIndices, newSourceRGCinputConeWeights, ...
            obj.localRGCRFspacingsMicrons(theSourceRGCindex));
    end

    % Find the (iSourceConeIndex, iNearbyRGCindex) pair that minimizes the total cost
    [~,iBestConeCombo] = min(projectedCostsSourceRGC(:));
    
   
    % The indices of the cones to be transfered
    theConeIndices = theSourceRGCinputConeIndices(C(iBestConeCombo,:));

    % Reassign the optical cone from the sourceRGC to the optimal nearbyRGC
    fprintf('Transfering % cones to a zero input RGC\n', numel(theConeIndices));
     
    % DISCONNECT theConeIndices from the current RGC
    obj.coneConnectivityMatrix(theConeIndices, theSourceRGCindex) = 0; % disconnect

    % And CONNECT the to the new RGC
    obj.coneConnectivityMatrix(theConeIndices, theDestinationRGCindex) = 1;
    
    % Update the centroids of the 2 RGCs
    affectedRGCindices = [theSourceRGCindex theDestinationRGCindex];
    obj.updateCentroidsFromInputs(affectedRGCindices);

end

