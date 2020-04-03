function connectionMatrix = computeConnectionMatrix(RGCRFPositionsMicrons, conePositionsMicrons, RGCRFSpacingsMicrons, desiredConesToRGCratios)

    % Define region of interest to work on
    roi.center = [100 0];
    roi.size = [100 70];

    % Find cones within the roi
    idxCones = positionsWithinROI(roi, conePositionsMicrons);
    conePositionsMicrons = conePositionsMicrons(idxCones,:);
    coneSpacingsMicrons = coneStats(conePositionsMicrons);
    
    % Instantiate a plotlab object
    plotlabOBJ = plotlab();

    % Apply the default plotlab recipe overriding 
    % the color order and the figure size
    plotlabOBJ.applyRecipe(...
        'colorOrder', [0 0 0; 1 0 0.5], ...
        'figureWidthInches', 25, ...
        'figureHeightInches', 15);
        
    
    % Find RGCs within the roi
    idxRGC = positionsWithinROI(roi, RGCRFPositionsMicrons);
    RGCRFPositionsMicrons = RGCRFPositionsMicrons(idxRGC,:);
    RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(idxRGC);
    desiredConesToRGCratios = desiredConesToRGCratios(idxRGC);
    
    % Step1. Align each RGC with its nearest cone. This ensure all RGC's
    % are connected to at least one cone. Since cones are more numerous
    % than RGCs some cones will not connect to an RGC at this step. This 
    % step occurs only for RGCs for which the cone-to-RGC ratio is [1..2]
    
    visualizeProcess = true;
    RGCRFPositionsMicrons = alignRGCmosaicToConeMosaic(...
        conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        desiredConesToRGCratios, visualizeProcess);
    
    % Step 2. Connect each cone to its neigboring RGC. Connection weights
    % depend on 2 factors: cone-to-RGC ration and proximity
    connectionMatrix = connectConesToRGC(conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        desiredConesToRGCratios, visualizeProcess);
end

function connectionMatrix = connectConesToRGC(conePositionsMicrons, coneSpacingsMicrons, ...
        RGCRFPositionsMicrons, RGCRFSpacingsMicrons, ...
        desiredConesToRGCratios, visualizeProcess)
    
    % Numbers of neurons
    conesNum = size(conePositionsMicrons,1);
    rgcsNum = size(RGCRFPositionsMicrons,1);
    
    
    % First pass. Connect each cone to its closest RGC. Since there are more cones than RGCs, some
    % RGCs will receive inputs from more than 1 cone in this pass.
    connectionMatrix = zeros(conesNum, rgcsNum,2);
    numberOfConeInputs = zeros(1,rgcsNum);
    for iCone = 1:conesNum
        % Compute distance of this cone to all RGCs
        conePosMicrons = conePositionsMicrons(iCone,:);
        distances = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, conePosMicrons).^2),2));
        [~, connectedRGCIndex] = min(distances);
        
        connectionMatrix(iCone, connectedRGCIndex,1) = 1;                            % strength
        connectionMatrix(iCone, connectedRGCIndex,2) = distances(connectedRGCIndex); % distance
        numberOfConeInputs(connectedRGCIndex) = numberOfConeInputs(connectedRGCIndex)+1;
    end % for iCone
    
     if (visualizeProcess)
        visualizeConnectivity(3, conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, mean(desiredConesToRGCratios))
     end
    
    % Second pass
    % Some RGCs will have N > 2 cone inputs. See if we can assign inputs
    % 3 ... N to other RGCs in the neigborhood that have the minimal number of
    % inputs <= N-1
    % 
    doSecondPass = true;
    if (doSecondPass)
        rgcIDsWithMoreThanTwoInputs = find(squeeze(sum(squeeze(connectionMatrix(:, :,1)),1)) > 2);
        fprintf('There are %d out of a total of %d RGCs that receive more than 2 cone inputs 1\n', ...
            numel(rgcIDsWithMoreThanTwoInputs), size(connectionMatrix,2));
        pause
    
        for iRGC = 1:numel(rgcIDsWithMoreThanTwoInputs)
            
            % Get RGC index
            rgcIndex = rgcIDsWithMoreThanTwoInputs(iRGC);
        
            % Get the indices of cones connected to this RGC
            indicesOfMultipleConesConnectedToThisRGC = find(squeeze(connectionMatrix(:, rgcIndex,1)) == 1);
            
            % Exclude the cone with the smallest distance. We dont want to
            % assign this one to another RGC
            [~,closestConeIdx] = min(squeeze(connectionMatrix(indicesOfMultipleConesConnectedToThisRGC, rgcIndex,2)));
            indicesOfMultipleConesConnectedToThisRGC = setdiff(indicesOfMultipleConesConnectedToThisRGC, ...
                indicesOfMultipleConesConnectedToThisRGC(closestConeIdx));
            
            % Go throuch each of these cones and see if we can assign some
            % to other RGCs
            iCone = 0;
            while (numberOfConeInputs(rgcIndex) > 2) && (iCone < numel(indicesOfMultipleConesConnectedToThisRGC))
                iCone = iCone+ 1;
                coneIndex = indicesOfMultipleConesConnectedToThisRGC(iCone);
                % Compute distance of this cone to all RGCs
                conePosMicrons = conePositionsMicrons(coneIndex,:);
                distances = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, conePosMicrons).^2),2));
                % Find indices of RGCs in neighborhood of thiscone
                rgcIndicesWithinReach = find(distances < 0.6*(max(RGCRFSpacingsMicrons)+coneSpacingsMicrons(coneIndex)));
                % Sort neigboring RGCs according to their distance to the cone
                [~, sortedIndices] = sort(distances(rgcIndicesWithinReach), 'ascend');
                
                % Find which RGC to connect to. This will be the closest RGC that
                % has the minimal # of cones already connected to it.
                % The minimum number of cone connections for the RGCs within reach 
                neighboringRGCConnectionsNum = sum(squeeze(connectionMatrix(:, rgcIndicesWithinReach,1)),1);
                minConnections = min(neighboringRGCConnectionsNum);
                
                keepGoing = true; k = 0;
                betterAlternativeRGCindex = nan;
                while (keepGoing) && (k < numel(sortedIndices))
                    k = k + 1;
                    theNeighboringRGCIndex = rgcIndicesWithinReach(sortedIndices(k));
                    if (neighboringRGCConnectionsNum(sortedIndices(k)) == minConnections)
                        keepGoing = false;
                        betterAlternativeRGCindex = theNeighboringRGCIndex;
                    end
                end
                
                if (~isnan(betterAlternativeRGCindex))
                    % Of DISCONNECT coneIndex from rgcIndex
                    connectionMatrix(coneIndex, rgcIndex, 1) = 0; % strength
                    connectionMatrix(coneIndex, rgcIndex, 2) = 0; % distance
                    numberOfConeInputs(rgcIndex) = numberOfConeInputs(rgcIndex)-1;
                    % And CONNECT coneIndex to betterAlternativeRGCindex
                    connectionMatrix(coneIndex, betterAlternativeRGCindex, 1) = 1; % strength
                    connectionMatrix(coneIndex, betterAlternativeRGCindex, 2) = distances(betterAlternativeRGCindex); % distance
                    numberOfConeInputs(betterAlternativeRGCindex) = numberOfConeInputs(betterAlternativeRGCindex)+1;
                end
            end % while
        end % for iRGC
    
    
    
        % Third pass. Some RGCs will have N = 2 cone inputs. If those inputs are from cones
        % of the same type we do nothing. If the inputs are from different cone
        % types see if we can assign the second cone to another RGC that has a
        % SINGLE input of the same cone type
        % NEED TO IMPLEMENT THIS
    
    end
    
    
    
    
    normalizeNetWeights = true;
    if (normalizeNetWeights)
        % Some RGCs will have more than one cones connected to them. Adjust the
        % connection strengths of those cones depending on proximity and
        % cone-to-RGC ratio for each RGC
        RGCswithZeroInputs = find(numberOfConeInputs == 0);
        for k = 1:numel(RGCswithZeroInputs)
            pos = RGCRFPositionsMicrons(RGCswithZeroInputs(k),:);
            fprintf('RGC at position %2.1f, %2.1f has zero cone inputs\n', pos(1), pos(2));
        end
        
        assert(isempty(RGCswithZeroInputs), 'There are RGCs with zero cone inputs');

        % Connection matrix for RGCs with more than 1 cone inputs
        RGCswithMoreThanOneInputs = find(numberOfConeInputs > 1);
        for k = 1:numel(RGCswithMoreThanOneInputs)
            theRGCindex = RGCswithMoreThanOneInputs(k);
            % Find cones connected to this RGC
            coneIndicesConnectedToThisRGC = find(squeeze(connectionMatrix(:, theRGCindex, 1)) == 1);
            
            % Adjust connection strengths
            distancesOfConesConnectedToThisRGC = squeeze(connectionMatrix(coneIndicesConnectedToThisRGC, theRGCindex,2));
            originalConnectionStrengths = squeeze(connectionMatrix(coneIndicesConnectedToThisRGC, theRGCindex, 1));
            weights = exp(-(distancesOfConesConnectedToThisRGC/RGCRFSpacingsMicrons(theRGCindex)));
            adjustedConeInputWeights = originalConnectionStrengths .* weights;
            % Net weight of all inputs should be 1.
            adjustedConeInputWeights = adjustedConeInputWeights / sum(adjustedConeInputWeights);
            % Update connection matrix
            connectionMatrix(coneIndicesConnectedToThisRGC, theRGCindex, 1) = adjustedConeInputWeights;
        end
    end
    
    
    if (visualizeProcess)
        visualizeConnectivity(4, conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, mean(desiredConesToRGCratios))
    end
    
    % Make sure that the marginal sums are all 1
    weights = squeeze(connectionMatrix(:,:,1));
    netWeights = sum(weights,1);
    if (any(netWeights>1+eps))
        disp('Some RGCs have net cone weights > 1')
    end
    if (any(netWeights<1-eps))
        disp('Some RGCs have net cone weights < 1')
    end
    
end

function visualizeConnectivity(figNo, conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, meanConesToRGCratio)
    hFig = figure(figNo); clf;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.04, ...
            'bottomMargin', 0.05);
    scatter(theAxesGrid{1,1}, conePositionsMicrons(:,1), conePositionsMicrons(:,2), 'b'); hold on;
    scatter(theAxesGrid{1,1}, RGCRFPositionsMicrons(:,1), RGCRFPositionsMicrons(:,2), 300, 'g');

    rgcsNum = size(RGCRFPositionsMicrons,1);
    
    for theRGCindex = 1:rgcsNum  
        rgcPos = RGCRFPositionsMicrons(theRGCindex,:);

        % Find cones connected to this RGC
        coneWeights = squeeze(connectionMatrix(:, theRGCindex, 1));
        coneIndicesConnectedToThisRGC = find(coneWeights > 0);
        coneWeights = coneWeights(coneIndicesConnectedToThisRGC);
        lineColors = coneWeights / max(coneWeights);

        for iCone = 1:numel(coneIndicesConnectedToThisRGC)
            theConeIndex = coneIndicesConnectedToThisRGC(iCone);
            conePos = conePositionsMicrons(theConeIndex,:);
            line([rgcPos(1) conePos(1)], [rgcPos(2) conePos(2)], ...
                'LineWidth', 1.5, 'Color', (1-lineColors(iCone))*[1 1 1]);
        end
    end
    title(sprintf('mean cone-to-RGC ratio: %2.2f', meanConesToRGCratio));
    
end






function  RGCRFPositionsMicrons = alignRGCmosaicToConeMosaic(...
    conePositionsMicrons, coneSpacingsMicrons, ...
    RGCRFPositionsMicrons, RGCRFSpacingsMicrons,...
    desiredConesToRGCratios, visualizeProcess)
   
    if (visualizeProcess)
        X1 = []; X2 = [];
        Y1 = []; Y2 = [];
        
        visualizeEachAlignment = false;
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
        
        if (visualizeProcess)
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
            
            X1 = cat(2, X1, RGCRFPositionMicronsBeforeAlignment(1));
            X2 = cat(2, X2, rgcPMicronsAligned(1));
            Y1 = cat(2, Y1, RGCRFPositionMicronsBeforeAlignment(2));
            Y2 = cat(2, Y2, rgcPMicronsAligned(2));
            
            if (iRGC == numel(indicesOfRGCsrequiringAlignment))
                hFig = figure(2); clf;
                theAxesGrid = plotlab.axesGrid(hFig, ...
                    'leftMargin', 0.04, ...
                    'bottomMargin', 0.05);

                scatter(theAxesGrid{1,1}, conePositionsMicrons(:,1), conePositionsMicrons(:,2), 'b'); hold on;
                scatter(theAxesGrid{1,1},RGCRFPositionsMicrons(:,1), RGCRFPositionsMicrons(:,2), 300, 'g');
            
                plot(theAxesGrid{1,1},[X1; X2], ...
                     [Y1; Y2], 'k-', 'LineWidth', 1.5);
                title(theAxesGrid{1,1},sprintf('cone-to-RGC ratio: %2.2f', mean(desiredConesToRGCratios)));
                drawnow;
                pause(0.1);
            end
        end % visualizeProcess
        
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