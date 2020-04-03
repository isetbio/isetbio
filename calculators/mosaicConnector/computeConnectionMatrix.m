function connectivityMatrix = computeConnectionMatrix(RGCRFPositionsMicrons, conePositionsMicrons, RGCRFSpacingsMicrons, desiredConesToRGCratios)

    % Define region of interest to work on
    roi.center = [1100 0];
    roi.size = [150 70];

    % Find cones within the roi
    idxCones = positionsWithinROI(roi, conePositionsMicrons);
    conesNum = numel(idxCones);
    conePositionsMicrons = conePositionsMicrons(idxCones,:);
    coneSpacingsMicrons = coneStats(conePositionsMicrons);
    
    
    % Find RGCs within the roi
    idxRGC = positionsWithinROI(roi, RGCRFPositionsMicrons);
    rgcsNum = numel(idxRGC);
    RGCRFPositionsMicrons = RGCRFPositionsMicrons(idxRGC,:);
    RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(idxRGC);
    desiredConesToRGCratios = desiredConesToRGCratios(idxRGC);
    
    % Step1. Align each RGC with its nearest cone. This ensure all RGC's
    % are connected to at least one cone. Since cones are more numerous
    % than RGCs some cones will not connect to an RGC at this step. This 
    % step occurs only for RGCs for which the cone-to-RGC ratio is [1..2]
    
    visualizeAlignmentProcess = true;
    RGCRFPositionsMicrons = alignRGCmosaicToConeMosaic(...
        conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        desiredConesToRGCratios, visualizeAlignmentProcess );
end

function  RGCRFPositionsMicrons = alignRGCmosaicToConeMosaic(...
    conePositionsMicrons, coneSpacingsMicrons, ...
    RGCRFPositionsMicrons, RGCRFSpacingsMicrons,...
    desiredConesToRGCratios, visualizeAlignmentProcess)
   
    if (visualizeAlignmentProcess)
        X1 = []; X2 = [];
        Y1 = []; Y2 = [];
        
        xOutline = cosd(0:10:360);
        yOutline = sind(0:10:360);
    
        % Instantiate a plotlab object
        plotlabOBJ = plotlab();
    
        % Apply the default plotlab recipe overriding 
        % the color order and the figure size
        plotlabOBJ.applyRecipe(...
            'colorOrder', [0 0 0; 1 0 0.5], ...
            'figureWidthInches', 25, ...
            'figureHeightInches', 15);
    end
    
    % Cones Num
    conesNum = size(conePositionsMicrons,1);
    rgcsNum = size(RGCRFPositionsMicrons,1);
    
    % Keep a track of cones that have already been aligned to some RGC
    coneAlignedWithRGCalready = false(1, conesNum);
    
    % Go through all the RGCs one by one
    indicesOfRGCsrequiringAlignment = find(desiredConesToRGCratios < 2);
    
    fprintf('Will align %d of %d RGCs in this patch, which had a cone-to-RGC ration < 2\n', ...
        numel(indicesOfRGCsrequiringAlignment), rgcsNum);
    
    for iRGC = 1:numel(indicesOfRGCsrequiringAlignment)
        
        % Get RGC index
        rgcIndex = indicesOfRGCsrequiringAlignment(iRGC);
        
        % Compute distance of this RGC to all the cones
        rgcPMicrons = RGCRFPositionsMicrons(rgcIndex,:);
        distances = sqrt(sum((bsxfun(@minus, conePositionsMicrons, rgcPMicrons ).^2),2));
        
        % Find indices of cones in neighborhood of this RGC
        coneIndicesWithinReach = find(distances < 0.55*(RGCRFSpacingsMicrons(rgcIndex)+max(coneSpacingsMicrons)));
        
        % Sort neigbors according to their distance to the RGC
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
        
        % If all neighboring cones are already aligned to other RGCs, just
        % align to the closest cone.
        if (isnan(alignmentConeIndex))
            fprintf(2, 'Could not find a cone to align this RGC that has not been aligned with another RGC already. Aligning to the closest cone.\n');
            alignmentConeIndex = coneIndicesWithinReach(sortedIndices(1));
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
        
        if (visualizeAlignmentProcess)
            visualizeEachAlignment = ~false;
            if (visualizeEachAlignment)
                hFig = figure(1); clf;
                theAxesGrid = plotlab.axesGrid(hFig, 'leftMargin', 0.03);

                xPts = rgcPMicrons(1)+0.5*RGCRFSpacingsMicrons(rgcIndex)*xOutline;
                yPts = rgcPMicrons(2)+0.5*RGCRFSpacingsMicrons(rgcIndex)*yOutline;
                plot(theAxesGrid{1,1}, rgcPMicrons(1)+0.5*RGCRFSpacingsMicrons(rgcIndex)*xOutline, ...
                     rgcPMicrons(2)+0.5*RGCRFSpacingsMicrons(rgcIndex)*yOutline, 'r-'); hold on; 

                plot(theAxesGrid{1,1}, rgcPMicronsAligned(1)+0.5*RGCRFSpacingsMicrons(rgcIndex)*xOutline, ...
                     rgcPMicronsAligned(2)+0.5*RGCRFSpacingsMicrons(rgcIndex)*yOutline, 'g-');
                hold(theAxesGrid{1,1}, 'on')
                for k = 1:numel(coneIndicesWithinReach)
                    coneIndex = coneIndicesWithinReach(k);
                    xPts = cat(2, xPts, conePositionsMicrons(coneIndex,1) + 0.5*coneSpacingsMicrons(coneIndex)*xOutline);
                    yPts = cat(2, yPts, conePositionsMicrons(coneIndex,2) + 0.5*coneSpacingsMicrons(coneIndex)*yOutline);
                    plot(theAxesGrid{1,1}, conePositionsMicrons(coneIndex,1) + 0.6*coneSpacingsMicrons(coneIndex)*xOutline, ...
                         conePositionsMicrons(coneIndex,2) + 0.6*coneSpacingsMicrons(coneIndex)*yOutline,'b-');
                end

                xMin = min(xPts); xMax = max(xPts); xRange = xMax-xMin;
                yMin = min(yPts); yMax = max(yPts); yRange = yMax-yMin;

                xLim = [xMin xMax];
                yLim = [yMin yMax];

                set(theAxesGrid{1,1}, 'XLim', xLim, 'YLim', yLim);
                axis(theAxesGrid{1,1}, 'equal')
                title(theAxesGrid{1,1},sprintf('cone-to-RGC ratio: %2.2f', desiredConesToRGCratios(rgcIndex)));
            end
            
            hFig = figure(2); clf;
            theAxesGrid = plotlab.axesGrid(hFig, ...
                'leftMargin', 0.04, ...
                'bottomMargin', 0.05);
            scatter(theAxesGrid{1,1}, conePositionsMicrons(:,1), conePositionsMicrons(:,2), 'b'); hold on;
            scatter(theAxesGrid{1,1},RGCRFPositionsMicrons(:,1), RGCRFPositionsMicrons(:,2), 300, 'g');
            X1 = cat(2, X1, RGCRFPositionMicronsBeforeAlignment(1));
            X2 = cat(2, X2, rgcPMicronsAligned(1));
            Y1 = cat(2, Y1, RGCRFPositionMicronsBeforeAlignment(2));
            Y2 = cat(2, Y2, rgcPMicronsAligned(2));
            plot(theAxesGrid{1,1},[X1; X2], ...
                 [Y1; Y2], 'k-', 'LineWidth', 1.5);
            title(theAxesGrid{1,1},sprintf('cone-to-RGC ratio: %2.2f', desiredConesToRGCratios(rgcIndex)))
            pause(0.1);
        end % visualizeAlignmentProcess
        
    end
    
    
end

    
function indices = positionsWithinROI(roi, positions)
    d = bsxfun(@minus,positions, roi.center);
    ecc = sqrt(sum(positions.^2,2));
    indices = find((abs(d(:,1)) <= 0.5*roi.size(1)) & (abs(d(:,2)) <= 0.5*roi.size(2)));
    
    % Re-order according to increasing eccentricity
    [~,sortedIdx] = sort(ecc(indices), 'ascend');
    indices = indices(sortedIdx);
end

function coneSpacings = coneStats(conePositions)
    p = pdist2(conePositions, conePositions, 'euclidean', 'Smallest', 3);
    p = p(2:end,:);
    coneSpacings = mean(p,1);
end