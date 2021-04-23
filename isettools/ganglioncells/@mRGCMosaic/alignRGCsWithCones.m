% Method to align RGC RFs to cones in the central retina.
% Called by obj.wireRFcenterToInputCones()
function alignRGCsWithCones(obj, coneRFpositionsMicrons, coneRFpositionsDegs, visualizeAlignment)

    if (visualizeAlignment)
        % Plot the mosaic before alignment
        hFig = figure(1);
        xyLims(1:2) = obj.eccentricityDegs(1) + obj.sizeDegs(1)/2 *[-1 1];
        xyLims(3:4) = obj.eccentricityDegs(2) + obj.sizeDegs(2)/2 *[-1 1];
    
        ax = subplot(2,2,1);
        obj.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', xyLims, ...
            'covisualizeInputConeMosaic', true, ...
            'plotTitle', 'before alignment');
    end
    
    phaseString = 'Aligning: ';
    
    % Set returned args
    rgcRFpositionsDegs = obj.rgcRFpositionsDegs;
    rgcRFpositionsMicrons = obj.rgcRFpositionsMicrons;
    rgcRFspacingsMicrons = obj.rgcRFspacingsMicrons;
    
    % Compute local ratio of cones to RGCs
    coneToRGCratios = obj.coneToRGCDensityRatios();
    
    % Align RGCs whose  cone-to-RGC ratio is less than 2
    indicesOfRGCsrequiringAlignment = find(coneToRGCratios < 2); 
    rgcsNum = size(obj.rgcRFpositionsMicrons,1);
    fprintf('%s Will align %d of %d RGCs in this patch (those RGCs have a cone-to-RGC ratio < 2)\n', ...
        phaseString, numel(indicesOfRGCsrequiringAlignment), rgcsNum);
    
    % Sort according to ecc
    ecc = sqrt(sum(obj.rgcRFpositionsMicrons(indicesOfRGCsrequiringAlignment,:).^2,2));
    [~,idx] = sort(ecc, 'ascend');
    indicesOfRGCsrequiringAlignment = indicesOfRGCsrequiringAlignment(idx);
    
    % We only align mRGCs to L/M cones
    eligibleConeIndices = find(...
        (obj.inputConeMosaic.coneTypes == cMosaic.LCONE_ID) | ...
        (obj.inputConeMosaic.coneTypes == cMosaic.MCONE_ID) );
    
    coneRFpositionsMicrons = coneRFpositionsMicrons(eligibleConeIndices,:);
    coneRFpositionsDegs = coneRFpositionsDegs(eligibleConeIndices,:);
    
    % Numbers of eligibe cones
    conesNum = size(coneRFpositionsMicrons,1);
    
    % Keep a track of cones that have already been aligned to some RGC
    coneAlignedWithRGCalready = false(1, conesNum);
    
    if (visualizeAlignment)
        X1 = []; X2 = [];
        Y1 = []; Y2 = [];
    end
    
    % Go through all the RGCs one by one
    notAlignedRGCs = 0;
    for iRGC = 1:numel(indicesOfRGCsrequiringAlignment)
        % Get RGC index
        rgcIndex = indicesOfRGCsrequiringAlignment(iRGC);
        
        % Compute distance of this RGC to all the cones
        rgcPDegs = rgcRFpositionsDegs(rgcIndex,:);
        rgcPMicrons = rgcRFpositionsMicrons(rgcIndex,:);
        distances = sqrt(sum((bsxfun(@minus, coneRFpositionsMicrons, rgcPMicrons).^2),2));
        
        % Find indices of cones in neighborhood of this RGC
        coneIndicesWithinReach = find(distances < 1.1*(rgcRFspacingsMicrons(rgcIndex)));
        
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
        
        % If still all neighboring cones are already aligned to other RGCs, keep
        % original position
        if (isnan(alignmentConeIndex))
            notAlignedRGCs = notAlignedRGCs+1;
            %fprintf('%s Could not find a cone to align RGC %d that has not been aligned with another RGC already. Keeping original position %2.1f,%2.1f)\n', phaseString, rgcIndex, rgcPMicrons(1), rgcPMicrons(2));
            continue;
        end
        
        % Position of the cone to be aligned to
        rgcPMicronsConeAligned = coneRFpositionsMicrons(alignmentConeIndex,:); 
        rgcPDegsConeAligned = coneRFpositionsDegs(alignmentConeIndex,:); 
        
        % Update position of RGC to get closer to  rgcPMicronsConeAligned
        % depending on desiredConesToRGCratios(rgcIndex);
        g = min([1 coneToRGCratios(rgcIndex)-1]);
        rgcPMicronsAligned = g * rgcPMicrons + (1-g)*rgcPMicronsConeAligned;
        rgcPDegsAligned = g * rgcPDegs + (1-g)*rgcPDegsConeAligned;
        
        % Update RGC mosaic
        rgcRFposBeforeAlignment = rgcRFpositionsMicrons(rgcIndex,:);
        rgcRFpositionsMicrons(rgcIndex,:) = rgcPMicronsAligned;
        rgcRFpositionsDegs(rgcIndex,:) = rgcPDegsAligned;  
        
        if (visualizeAlignment)
            % Update alignment vectors
            X1 = cat(2, X1, rgcRFposBeforeAlignment(1));
            X2 = cat(2, X2, rgcPMicronsAligned(1));
            Y1 = cat(2, Y1, rgcRFposBeforeAlignment(2));
            Y2 = cat(2, Y2, rgcPMicronsAligned(2));
        end
        
    end % iRGC
    

    % Update positions and spacings of RGC RFs
    obj.rgcRFpositionsDegs = rgcRFpositionsDegs;
    obj.rgcRFpositionsMicrons = rgcRFpositionsMicrons;
    obj.rgcRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(rgcRFpositionsMicrons);
    obj.rgcRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(rgcRFpositionsDegs);
    
    if (visualizeAlignment)
        % Plot the alignment vectors
        ax = subplot(2,2,[3 4]);
        scatter(ax, coneRFpositionsMicrons(:,1), coneRFpositionsMicrons(:,2), 'k'); hold on;
        scatter(ax, obj.rgcRFpositionsMicrons(:,1), obj.rgcRFpositionsMicrons(:,2), 'r');
        plot(ax, [X1; X2], [Y1; Y2], 'k-', 'LineWidth', 1.5);
    end
    
    %Feedback
    fprintf('%s Could not find a cone that was not already aligned to another RGC for %d out of %d RGCs\n', phaseString, notAlignedRGCs, numel(indicesOfRGCsrequiringAlignment));
    
    if (visualizeAlignment)
        % Plot the mosaic after the alignment
        ax = subplot(2,2,2);
        obj.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', xyLims, ...
            'covisualizeInputConeMosaic', true, ...
            'plotTitle', 'after alignment');
    end
    
end
