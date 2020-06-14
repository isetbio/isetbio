function  RGCRFPositionsMicrons = alignRGCmosaicToConeMosaic(...
    conePositionsMicrons, coneSpacingsMicrons, ...
    RGCRFPositionsMicrons, RGCRFSpacingsMicrons,...
    coneTypes, desiredConesToRGCratios, visualizeProcess)
   
    phaseString = 'Aligning: ';
    
    if (visualizeProcess)
        X1 = []; X2 = [];
        Y1 = []; Y2 = [];
        
        visualizeEachAlignment = ~true;
        if (visualizeEachAlignment)
            xOutline = cosd(0:10:360);
            yOutline = sind(0:10:360);
        end
    end
    
    
    % Align RGC with cones only for those RGCs whose desired cone-to-RGC
    % ratio is less than 2
    indicesOfRGCsrequiringAlignment = find(desiredConesToRGCratios < 2); 
    rgcsNum = size(RGCRFPositionsMicrons,1);
    fprintf('%s Will align %d of %d RGCs in this patch (those RGCs have a cone-to-RGC ratio < 2)\n', ...
        phaseString, numel(indicesOfRGCsrequiringAlignment), rgcsNum);

    % Sort according to ecc
    ecc = sqrt(sum(RGCRFPositionsMicrons(indicesOfRGCsrequiringAlignment,:).^2,2));
    [~,idx] = sort(ecc, 'ascend');
    indicesOfRGCsrequiringAlignment = indicesOfRGCsrequiringAlignment(idx);
    
    % We only align mRGCs to L/M cones, not S
    eligibleConeIndices = find(coneTypes < 4);
    conePositionsMicrons = conePositionsMicrons(eligibleConeIndices,:);
    
    % Numbers of eligibe cones
    conesNum = size(conePositionsMicrons,1);
    
    % Keep a track of cones that have already been aligned to some RGC
    coneAlignedWithRGCalready = false(1, conesNum);
      
    % Go through all the RGCs one by one
    notAlignedRGCs = 0;
    for iRGC = 1:numel(indicesOfRGCsrequiringAlignment)
        % Get RGC index
        rgcIndex = indicesOfRGCsrequiringAlignment(iRGC);

        % Compute distance of this RGC to all the cones
        rgcPMicrons = RGCRFPositionsMicrons(rgcIndex,:);
        distances = sqrt(sum((bsxfun(@minus, conePositionsMicrons, rgcPMicrons ).^2),2));
        
        % Find indices of cones in neighborhood of this RGC
        coneIndicesWithinReach = find(distances < 0.8*(RGCRFSpacingsMicrons(rgcIndex)));
        
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
            notAlignedRGCs = notAlignedRGCs+1;
            %fprintf('%s Could not find a cone to align RGC %d that has not been aligned with another RGC already. Keeping original position %2.1f,%2.1f)\n', phaseString, rgcIndex, rgcPMicrons(1), rgcPMicrons(2));
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
    end % iRGC 
    
    fprintf('%s Could not find a cone that was not already aligned to another RGC for %d out of %d RGCs\n', phaseString, notAlignedRGCs, numel(indicesOfRGCsrequiringAlignment));
           
end
