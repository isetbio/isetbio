function  RGCRFPositionsMicrons = alignRGCmosaicToConeMosaic(...
    conePositionsMicrons, coneSpacingsMicrons, ...
    RGCRFPositionsMicrons, RGCRFSpacingsMicrons,...
    desiredConesToRGCratios, visualizeProcess)
   
    if (visualizeProcess)
        X1 = []; X2 = [];
        Y1 = []; Y2 = [];
        
        visualizeEachAlignment = ~true;
        if (visualizeEachAlignment)
            xOutline = cosd(0:10:360);
            yOutline = sind(0:10:360);
        end
    end
    
    % Numbers of neurons
    conesNum = size(conePositionsMicrons,1);
    rgcsNum = size(RGCRFPositionsMicrons,1);
    
    % Keep a track of cones that have already been aligned to some RGC
    coneAlignedWithRGCalready = false(1, conesNum);
    
    % Align RGC with cones only for those RGCs whose desired cone-to-RGC
    % ratio is less than 2
    indicesOfRGCsrequiringAlignment = find(desiredConesToRGCratios < 2); 
    fprintf('Will align %d of %d RGCs in this patch, which had a cone-to-RGC ration < 2\n', ...
        numel(indicesOfRGCsrequiringAlignment), rgcsNum);
    
    % Go through all the RGCs one by one
    for iRGC = 1:numel(indicesOfRGCsrequiringAlignment)
        
        % Get RGC index
        rgcIndex = indicesOfRGCsrequiringAlignment(iRGC);
        
        % Compute distance of this RGC to all the cones
        rgcPMicrons = RGCRFPositionsMicrons(rgcIndex,:);
        distances = sqrt(sum((bsxfun(@minus, conePositionsMicrons, rgcPMicrons ).^2),2));
        
        % Find indices of cones in neighborhood of this RGC
        coneIndicesWithinReach = find(distances < 0.55*(RGCRFSpacingsMicrons(rgcIndex)+max(coneSpacingsMicrons)));
        
        % Sort neigboring cones according to their distance to the RGC
        [~, sortedIndices] = sort(distances(coneIndicesWithinReach));
        
        % Find which cone (within the above neigborhhod) to align to. This
        % will be the closest cone that has not already been aligned with
        % another RGC
        keepGoing = true; k = 0;
        alignmentConeIndex = nan;
        while (keepGoing) && (k < numel(sortedIndices))
            k = k + 1;
            theConeIndex = coneIndicesWithinReach(sortedIndices(k));
            if (~coneAlignedWithRGCalready(theConeIndex))
                keepGoing = false;
                alignmentConeIndex = theConeIndex;
                coneAlignedWithRGCalready(theConeIndex) = true;
            end
        end
        
        % If all neighboring cones are already aligned to other RGCs, keep
        % original position
        if (isnan(alignmentConeIndex))
            fprintf('Could not find a cone to align RGC %d that has not been aligned with another RGC already. Keeping original position\n', rgcIndex);
            continue;
        end
        
        % Position of the cone to be aligned to
        rgcPMicronsConeAligned = conePositionsMicrons(alignmentConeIndex,:); 
        
        % Update position of RGC to get closer to  rgcPMicronsConeAligned
        % depending on desiredConesToRGCratios(rgcIndex);
        g = min([1 desiredConesToRGCratios(rgcIndex)-1]);
        rgcPMicronsAligned = g * rgcPMicrons + (1-g)*rgcPMicronsConeAligned;
        
        % Update RGC mosaic
        RGCRFPositionMicronsBeforeAlignment = RGCRFPositionsMicrons(rgcIndex,:);
        RGCRFPositionsMicrons(rgcIndex,:) = rgcPMicronsAligned;
        
        if (visualizeProcess)
            
            if (visualizeEachAlignment)
                visualizeCurrentRGCalignment(rgcPMicrons, rgcPMicronsAligned, RGCRFSpacingsMicrons(rgcIndex), ...
                    xOutline, yOutline, coneIndicesWithinReach, conePositionsMicrons, coneSpacingsMicrons, desiredConesToRGCratios);
            end
            
            % Update alignment vectors
            X1 = cat(2, X1, RGCRFPositionMicronsBeforeAlignment(1));
            X2 = cat(2, X2, rgcPMicronsAligned(1));
            Y1 = cat(2, Y1, RGCRFPositionMicronsBeforeAlignment(2));
            Y2 = cat(2, Y2, rgcPMicronsAligned(2));
            
            if (iRGC == numel(indicesOfRGCsrequiringAlignment))
                visualizeFinalAlignment(conePositionsMicrons, RGCRFPositionsMicrons, X1, X2, Y1, Y2, desiredConesToRGCratios);
            end
        end % visualizeProcess
        
    end
end
