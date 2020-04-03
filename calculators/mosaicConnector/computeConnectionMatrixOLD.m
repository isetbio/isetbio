function connectivityMatrix = computeConnectionMatrixOLD(RGCRFPositionsMicrons, conePositionsMicrons, RGCRFSpacingsMicrons, desiredConesToRGCratios)

    % Define region of interest to work on
    roi.center = [300 0];
    roi.radius = 30;

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
    
    
    reUseSynapses = ~true;
    
    % Generate cone synapses
    synapsesPerCone = 30;
    synapseSpreadFactor = 0.15;  % synapse spread sigma = synapseSpreadFactor * coneSpacing
    allConeSynapseLocs = generateConeSynapses(conePositionsMicrons, coneSpacingsMicrons, ...
        desiredConesToRGCratios, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, synapsesPerCone, synapseSpreadFactor);
    
    % Generate RGC pooling outlines
    poolingSpreadFactor = 0.5;  % pooling radius = poolingSpreadFactor * rgcSpacing
    rgcPoolingOutlines = generateRGCPooolingAreas(RGCRFPositionsMicrons, RGCRFSpacingsMicrons, poolingSpreadFactor);
    
    % Connect RGCs to cone synapses
    RGCconnectionsToSynapses = cell(1, rgcsNum);
    availableConeSynapseLocs = allConeSynapseLocs;
    
    for k = 1:rgcsNum 
        % Find the indices of available cone synapses that this RGC will
        % connect to
        [indicesOfConnectedConeSynapses, availableConeSynapseLocs] = ...
            indicesOfConnectedSynapses(availableConeSynapseLocs, squeeze(rgcPoolingOutlines(k,:,:)), ...
            desiredConesToRGCratios(k), synapsesPerCone, reUseSynapses);
        
        % Save cone synapses connected to this RGC
        RGCconnectionsToSynapses{k} = indicesOfConnectedConeSynapses;
        
    end
    
    % Constraint #2. Minimize overlap of connections from one cone to
    % multiple RGCs
    for k = 1:conesNum
        
    end
    
    % Compute connectivity matrix
    connectivityMatrix = zeros(rgcsNum, conesNum);
    for k = 1:rgcsNum
        indicesOfConnectedConeSynapses = RGCconnectionsToSynapses{k};
        % Find the cones to which each of the closest synapse belongs to
        coneIDs = floor((indicesOfConnectedConeSynapses-1)/synapsesPerCone)+1;

        % Form connectivity matrix
        for ll = 1:numel(coneIDs)
            connectivityMatrix(k, coneIDs(ll)) = connectivityMatrix(k, coneIDs(ll)) + 1;
        end
    end
    
    % Total connection weighgts to each RGC is 1
    connectivityMatrix = bsxfun(@times, connectivityMatrix,  1./ sum(connectivityMatrix,2));
    
     
    % Plot results
    % Instantiate a plotlab object
    plotlabOBJ = plotlab();
    
    % Apply the default plotlab recipe overriding 
    % the color order and the figure size
    plotlabOBJ.applyRecipe(...
        'colorOrder', [0 0 0; 1 0 0.5], ...
        'figureWidthInches', 15, ...
        'figureHeightInches', 15);
    
    % Visualize how connection of cones to RGC happens
    plotConnectionProcedure(2, allConeSynapseLocs, conePositionsMicrons, desiredConesToRGCratios, RGCRFPositionsMicrons, RGCconnectionsToSynapses, rgcPoolingOutlines, reUseSynapses, roi);
    
    % New figure
    hFig = figure(1); clf;
    
    plot(allConeSynapseLocs(:,1), allConeSynapseLocs(:,2), '.'); hold on;
    plotConeDendrite(conePositionsMicrons,  allConeSynapseLocs);
    
    for k = 1:numel(idxRGC)
  
        indicesOfConnectedConeSynapses = RGCconnectionsToSynapses{k};
        scatter(allConeSynapseLocs(indicesOfConnectedConeSynapses,1), allConeSynapseLocs(indicesOfConnectedConeSynapses,2), 50, 'r');
           
        x = RGCRFPositionsMicrons(k,1);
        y = RGCRFPositionsMicrons(k,2);
        scatter(x,y, 'g');
        text(x,y,sprintf('%2.2f', desiredConesToRGCratios(k)));
  
    end
    title(sprintf('cones/RGC ratio in patch: %2.2f', numel(idxCones)/numel(idxRGC)));
    axis 'square'
    set(gca, 'XLim', roi.center(1)+1.1*roi.radius*[-1 1], 'YLim', roi.center(2)+1.1*roi.radius*[-1 1]);
    
    
    % Plot RF's
    % Generate RGC RFs based on connectivityMatrix
    [RGCrfs, RGCRFsupportX, RGCRFsupportY] = generateRGCRFsFromConnectivityMatrix(...
        conePositionsMicrons, coneSpacingsMicrons, connectivityMatrix, roi);
  
    figure(3); clf;
    scatter(conePositionsMicrons(:,1), conePositionsMicrons(:,2)); hold on
    scatter(RGCRFPositionsMicrons(:,1), RGCRFPositionsMicrons(:,2), 200);
    
    zLevels = exp([-2 -1]);
    whichLevelsToContour = [2];

    zLevels = 0.25 : 0.25 : 1.0;
    whichLevelsToContour = [1];
    
    for RGCindex = 1:numel(idxRGC)
        
        C = contourc(RGCRFsupportX, RGCRFsupportY,squeeze(RGCrfs(RGCindex,:,:)), zLevels);
        k = 1;
        while k < size(C,2)
            level = C(1,k);
            points = C(2,k);
            for kLevel = 1:numel(zLevels)
                if (level == zLevels(kLevel))
                    xRGCEnsembleOutline.level{kLevel} = C(1,k+(1:points));
                    yRGCEnsembleOutline.level{kLevel} = C(2,k+(1:points));
                end
            end
            k = k+points+1;
        end
    
        for iLevel = 1:numel(whichLevelsToContour)
            theLevel = whichLevelsToContour(iLevel);
            maxLevel = max(whichLevelsToContour);
            f = 1:numel(xRGCEnsembleOutline.level{theLevel});
            v = [xRGCEnsembleOutline.level{theLevel}(:) yRGCEnsembleOutline.level{theLevel}(:)];
            patch('Faces', f, 'Vertices', v, 'FaceColor', (1-0.8*theLevel/maxLevel)*[0.8 0.8 0.8], ...
                'FaceAlpha', 0.2, 'EdgeColor', [0 0 0.6], 'EdgeAlpha', 0.6, 'LineWidth', 1.0);
        end
    	
        maxConnectionStrength = max(squeeze(connectivityMatrix(RGCindex,:)));
        
        for coneIndex = 1:numel(idxCones)
            relativeConnectionStrength = connectivityMatrix(RGCindex, coneIndex)/maxConnectionStrength;
            if (relativeConnectionStrength>=0.05)
                plot([RGCRFPositionsMicrons(RGCindex,1) conePositionsMicrons(coneIndex,1)], ...
                 [RGCRFPositionsMicrons(RGCindex,2) conePositionsMicrons(coneIndex,2)], ...
                 'r-', 'LineWidth', 10*relativeConnectionStrength);
            end
        end
        
    end
     set(gca, 'XLim', roi.center(1)+1.1*roi.radius*[-1 1], 'YLim', roi.center(2)+1.1*roi.radius*[-1 1]);
    
     
end


function [RGCrfs, xAxis, yAxis] = generateRGCRFsFromConnectivityMatrix(conePositions, coneSpacings, connectivityMatrix, roi)
    % Sampling for contours
    deltaX = 0.1;
    xAxis = (roi.center(1)-roi.radius): deltaX: (roi.center(1)+roi.radius);
    yAxis = (roi.center(2)-roi.radius): deltaX: (roi.center(2)+roi.radius);
    [X,Y] = meshgrid(xAxis,yAxis);
    
    conesNum = size(connectivityMatrix,2);
    rgcsNum = size(connectivityMatrix,1);
    
    coneProfiles = zeros(conesNum, size(X,1), size(X,2));
    
    for coneIndex = 1:conesNum 
        cP = squeeze(conePositions(coneIndex,:));
        coneSigma = coneSpacings(coneIndex)/3;
        
        coneProfile = exp(-0.5*((X-cP(1))/coneSigma).^2) .* exp(-0.5*((Y-cP(2))/coneSigma).^2);
        coneProfile = coneProfile.^0.9;
        coneProfiles(coneIndex,:,:) = coneProfile/max(coneProfile(:));
    end
    
    RGCrfs = zeros(rgcsNum, size(X,1), size(X,2));
    for RGCindex = 1:rgcsNum
        RGCrf = squeeze(RGCrfs(RGCindex,:,:));
        for coneIndex = 1:conesNum
            if (connectivityMatrix(RGCindex, coneIndex)>0)
                RGCrf = RGCrf + squeeze(coneProfiles(coneIndex,:,:)) * connectivityMatrix(RGCindex, coneIndex);
            end
        end
        RGCrfs(RGCindex,:,:) = RGCrf / max(RGCrf(:));
    end
end





function [indicesOfConnectedConeSynapses, availableConeSynapseLocs] = ...
            indicesOfConnectedSynapses(availableConeSynapseLocs, rgcPoolingOutline, ...
            desiredConeToRGCratio, synapsesPerCone, reUseSynapses)
        
    % Find all cone synapses that are within the rgcPoolingOutline
    [insideIndices, onIndices] = inpolygon(...
        availableConeSynapseLocs(:,1), availableConeSynapseLocs(:,2), ...
        squeeze(rgcPoolingOutline(1,:)), squeeze(rgcPoolingOutline(2,:)));
        
    % Combine synapses that are inside and on perimeter
    allIndices = insideIndices | onIndices; 
    indicesOfConnectedConeSynapses = find(allIndices(:));
    
    % Find the cones to which each of the synapses belong to
    coneIDs = floor((indicesOfConnectedConeSynapses-1)/synapsesPerCone)+1;
    uniqueConeIDs = unique(coneIDs);
    
    % Contraint #1. Adhere to desiredConeToRGC ratio.
    % Drop inputs from cones depending based on this RGC's desiredConeToRGCratio
%    connectedConesNum = numel(uniqueConeIDs);
%     if (connectedConesNum > desiredConeToRGCratio)
%         
%         % Sort cones according to the synapses from them
%         connectionStrengths = zeros(1,numel(uniqueConeIDs));
%         %fprintf('\n--------------\n');
%         for k = 1:numel(uniqueConeIDs)
%             connectionStrengths(k) = numel(find(coneIDs == uniqueConeIDs(k)));
%             %fprintf('Connection with coneID %d is via %d synapses\n', uniqueConeIDs(k), connectionStrengths(k));
%         end
%         % Reorder connection strengths with decreasing order
%         [~, idx] = sort(connectionStrengths, 'descend');
%         connectionStrengths = connectionStrengths(idx);
%         uniqueConeIDs = uniqueConeIDs(idx);
%         
%         % We will probabilistically remove synapses from cones > desiredConeToRGCratio
%         conesIDsToBeProbabilisticallySampled = uniqueConeIDs(floor(desiredConeToRGCratio)+1:numel(uniqueConeIDs));
%         for c = 1:numel(conesIDsToBeProbabilisticallySampled)
%             for cID = 1:numel(coneIDs)
%                 if (coneIDs(cID)== conesIDsToBeProbabilisticallySampled(c))
%                     % See if we will keep this synapse
%                     pKeep = desiredConeToRGCratio-floor(desiredConeToRGCratio);
%                     if (rand(1,1) > pKeep)
%                         % Label this synapse for removal
%                         % fprintf('removing synapse from coneID: %d\n', coneIDs(cID));
%                         indicesOfConnectedConeSynapses(cID) = nan;
%                     end
%                 end
%             end
%         end
%    
%         % Update the indices of connected cone synapses
%         indicesOfConnectedConeSynapses = indicesOfConnectedConeSynapses(find(~isnan(indicesOfConnectedConeSynapses)));
% 
%         % Find the cones to which each of the synapses belong to
%         coneIDs = floor((indicesOfConnectedConeSynapses-1)/synapsesPerCone)+1;
%         uniqueConeIDs = unique(coneIDs);
%     end % (connectedConesNum > desiredConeToRGCratio)

    % Determine connection strength to all connected cones based on number of 
    % synapses with those cones
    connectionStrengths = zeros(1,numel(uniqueConeIDs));
    for k = 1:numel(uniqueConeIDs)
        connectionStrengths(k) = numel(find(coneIDs == uniqueConeIDs(k)));
        %fprintf('After adjusting for coneToRGC ratio, connection with coneID %d is via %d synapses\n', uniqueConeIDs(k), connectionStrengths(k));
    end
        
    
    % Make the positions of these synapses Inf, so no other RGC can connect to them
    if (~reUseSynapses)
        availableConeSynapseLocs(indicesOfConnectedConeSynapses,:) = Inf;
    end
end

function poolingOutlines = generateRGCPooolingAreas(positions, spacings, spreadFactor)
    unitOutline(:,1)= cosd(0:30:360);
    unitOutline(:,2) = sind(0:30:360);
    rgcNum = size(positions,1);
    for k = 1:rgcNum
        poolingOutlines(k,1,:) = positions(k,1) + spreadFactor*spacings(k)*unitOutline(:,1);
        poolingOutlines(k,2,:) = positions(k,2) + spreadFactor*spacings(k)*unitOutline(:,2);
    end
    
end


function plotConnectionProcedure(figNo, allConeSynapseLocs, conePositionsMicrons, desiredConesToRGCratios, RGCRFPositionsMicrons, RGCconnectionsToSynapses, rgcPoolingOutlines, reUseSynapses, roi)
    figure(figNo); clf; 

    availableSynapseLocs = allConeSynapseLocs;
    plotConeDendrite(conePositionsMicrons,  allConeSynapseLocs);
    scatter(RGCRFPositionsMicrons(:,1),RGCRFPositionsMicrons(:,2), 'g');
    
    % Plot the available synapses
    p1 = plot(availableSynapseLocs(:,1), availableSynapseLocs(:,2), 'k.');
    
    rgcsNum = numel(RGCconnectionsToSynapses);
    for k = 1:rgcsNum
        
        x = RGCRFPositionsMicrons(k,1);
        y = RGCRFPositionsMicrons(k,2);
        text(x,y,sprintf('%2.2f', desiredConesToRGCratios(k)));
        
        % Plot the RGC pooling outline in red
        fill(squeeze(rgcPoolingOutlines(k,1,:)), squeeze(rgcPoolingOutlines(k,2,:)), 'r', 'FaceAlpha', 0.4);
        
        % Remove used synapses from plotting
        indicesOfConnectedConeSynapses = RGCconnectionsToSynapses{k};
        if (~reUseSynapses)
            availableSynapseLocs(indicesOfConnectedConeSynapses,:) = Inf;
        end
        
        % Plot the connected synapses in blue
        p2 = scatter(allConeSynapseLocs(indicesOfConnectedConeSynapses,1), ...
                allConeSynapseLocs(indicesOfConnectedConeSynapses,2), 64, 'bo', 'LineWidth', 1.0);
        set(gca, 'XLim', roi.center(1) +1.1*roi.radius*[-1 1], ...
                 'YLim', roi.center(2) +1.1*roi.radius*[-1 1]);
             
        drawnow
        pause(0.2);
        set(p1, 'XData', availableSynapseLocs(:,1), 'YData', availableSynapseLocs(:,2));
        set(p2, 'XData', [], 'YData', []);
    end
end
    
function plotConeDendrite(conePositions,  allConeSynapsePositions)
    conesNum = size(conePositions,1);
    coneSynapsesNum = size(allConeSynapsePositions,1);
    synapsesPerCone = coneSynapsesNum/conesNum;
    hold on;
    for l = 1:conesNum
        kkk = (l-1)*synapsesPerCone;
        X1 = repmat(conePositions(l,1), [1 synapsesPerCone]);
        X2 = allConeSynapsePositions(kkk+(1:synapsesPerCone),1);
        Y1 = repmat(conePositions(l,2), [1 synapsesPerCone]);
        Y2 = allConeSynapsePositions(kkk+(1:synapsesPerCone),2);
        plot([X1; X2'], [Y1; Y2'], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    end
end


function allConeSynapseLocs = generateConeSynapses(conePositionsMicrons, coneSpacingsMicrons, desiredConesToRGCratios, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, synapsesPerCone, synapseSpreadFactor)
    conesNum = size(conePositionsMicrons,1);
    allConeSynapseLocs = zeros(synapsesPerCone*conesNum,2);
    rgcsNum = size(RGCRFPositionsMicrons,1);
    
    rgcHasAlreadyAttractedConeSynapses = zeros(1,rgcsNum);
    
    figure(100); clf;  
    % Generate Gaussian cloud of synapses
    for k = 1:conesNum
        
        theConePosition = conePositionsMicrons(k,:);
        theNearestRGCBiasedPosition = nan;
        

        scatter(theConePosition(1), theConePosition(2), 'k^'); hold on
            
        % Decide if we will adjust position of synapses
        biasConeSynapsePositionsTowardsClosestRGC = true;
        
        if (biasConeSynapsePositionsTowardsClosestRGC)
            % Sort RGCs according to their distance from the cone position
            diffs = sqrt(sum((bsxfun(@minus, RGCRFPositionsMicrons, theConePosition).^2),2));
            [distancesToNearestRGCs, nearestRGCindices] = sort(diffs, 'ascend');

            % Determine search radius for selecting the nearby RGCs
            searchRadius = 1.2*RGCRFSpacingsMicrons(nearestRGCindices(1));
            idx = find(distancesToNearestRGCs < searchRadius);
            nearestRGCindices = nearestRGCindices(idx);
           
            scatter(RGCRFPositionsMicrons(nearestRGCindices,1), RGCRFPositionsMicrons(nearestRGCindices,2), 'r');
            for kkkk = 1:numel(nearestRGCindices)
                theRGCindex = nearestRGCindices(kkkk);
                text(RGCRFPositionsMicrons(theRGCindex ,1), RGCRFPositionsMicrons(theRGCindex ,2), sprintf('%2.2f',desiredConesToRGCratios(theRGCindex)));
            end
            
            % Find nearest RGCs that have not already attracted synapses from some cone
            kk = 0;
            keepGoing = true;
            while (kk < numel(nearestRGCindices))&&(keepGoing)
                kk = kk + 1;
                theRGCindex = nearestRGCindices(kk);
                if (rgcHasAlreadyAttractedConeSynapses(theRGCindex)==0)
                    g = 0.998;
                    theNearestRGCBiasedPosition = (1-g)*theConePosition + g*RGCRFPositionsMicrons(theRGCindex,:);
                    rgcHasAlreadyAttractedConeSynapses(theRGCindex) = 1;
                    keepGoing = false;
                else
                   % depending on the current cone-to-RGC ratio, decide if
                   % we can add some of these cone inputs to an RGC that already had cone
                   % inputs
                   threshold = desiredConesToRGCratios(theRGCindex)-rgcHasAlreadyAttractedConeSynapses(theRGCindex);
                   if (rand(1,1) < threshold)
                        g = 0.998;
                        theNearestRGCBiasedPosition = (1-g)*theConePosition + g*RGCRFPositionsMicrons(theRGCindex,:);
                        rgcHasAlreadyAttractedConeSynapses(theRGCindex) = rgcHasAlreadyAttractedConeSynapses(theRGCindex) + 1;
                        keepGoing = false;
                   end
                end
            end
            
            for kk = 1:numel(nearestRGCindices)
                theRGCindex = nearestRGCindices(kk);
                if (rgcHasAlreadyAttractedConeSynapses(theRGCindex)> 0)
                   scatter(RGCRFPositionsMicrons(theRGCindex,1), RGCRFPositionsMicrons(theRGCindex,2), rgcHasAlreadyAttractedConeSynapses(theRGCindex)*150, 'y');
                end
            end
            
        end
        
        if (isnan(theNearestRGCBiasedPosition))
            theNearestRGCBiasedPosition = theConePosition;
        end
        
        scatter(theNearestRGCBiasedPosition(1),theNearestRGCBiasedPosition(2), 'g');
        
        hold off;
        set(gca, 'YLim', [-30 30], 'XLim', [270 330]);
        pause
        % Generate cloud of synapses around (possibly adjusted) cone position 
        synapseSpread = synapseSpreadFactor*coneSpacingsMicrons(k);
        
        allConeSynapseLocs((k-1)*synapsesPerCone+(1:synapsesPerCone),:) = ...
                bsxfun(@plus, theNearestRGCBiasedPosition, randn(synapsesPerCone,2) * synapseSpread);
        
    end
end

function indices = positionsWithinROI(roi, positions)
    d = bsxfun(@minus,positions, roi.center);
    d = sum(d.^2,2)/(roi.radius)^2;
    % Re-order according to increasing eccentricity
    [d,sortedIdx] = sort(d, 'ascend');
    idx = find(d < 1);
    indices = sortedIdx(idx);
end