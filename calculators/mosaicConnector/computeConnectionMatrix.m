function connectionMatrix = computeConnectionMatrix(RGCRFPositionsMicrons, conePositionsMicrons, RGCRFSpacingsMicrons, desiredConesToRGCratios)

    xo = 300;
    yo = 300;
    radius = 20;

    d = bsxfun(@minus,conePositionsMicrons, [xo yo]);
    d = sqrt(sum(d.^2,2))/radius;
    idxCones = find(d < 1);
    conePositionsMicrons = conePositionsMicrons(idxCones,:);
    coneSpacingsMicrons = coneStats(conePositionsMicrons);
    
    d = bsxfun(@minus,RGCRFPositionsMicrons, [xo yo]);
    d = sqrt(sum(d.^2,2))/radius;
    idxRGC = find(d < 1);
    RGCRFPositionsMicrons = RGCRFPositionsMicrons(idxRGC,:);
    RGCRFSpacingsMicrons = RGCRFSpacingsMicrons(idxRGC);
    desiredConesToRGCratios = desiredConesToRGCratios(idxRGC);
    
    
    
    synapsesPerCone = 100;
    allConeSynapseLocs = zeros(synapsesPerCone*numel(idxCones),2);
    % Generate Gaussian cloud of 100 cone synapses  with a sigma = 0.2*cone spacing
    for k = 1:numel(idxCones)
        allConeSynapseLocs((k-1)*synapsesPerCone+(1:synapsesPerCone),:) = ...
            bsxfun(@plus, conePositionsMicrons(k,:), randn(synapsesPerCone,2) * 0.2*coneSpacingsMicrons(k));
    end
    
    allConeSynapseLocsOriginal = allConeSynapseLocs;
    
    connectivityMatrix = zeros(numel(idxRGC), numel(idxCones));
    RGCconnections = cell(1, numel(idxRGC));

    figure(2); clf;
    for k = 1:numel(idxRGC)
        plot(allConeSynapseLocs(:,1), allConeSynapseLocs(:,2), 'k.'); hold on;
        
        drawnow;
        
        connectionApproach = 'new';
        if (connectionApproach == 'new')
            
            % Pooling area for cone synapses
            rgcPoolingOutline.x = RGCRFPositionsMicrons(k,1) + 0.6*RGCRFSpacingsMicrons(k)*cosd(0:10:360);
            rgcPoolingOutline.y = RGCRFPositionsMicrons(k,2) + 0.6*RGCRFSpacingsMicrons(k)*sind(0:10:360);
            fill(rgcPoolingOutline.x, rgcPoolingOutline.y, 'r', 'FaceAlpha', 0.4);
            drawnow;
            
            % find all cone synapses that are within the rgcPoolingOutline
            
        end
        
        if (connectionApproach == 'old')
            % Generate a Gaussian cloud of RGC synapses (N = cone synapses * cone-to-RGC ratio) with a sigma = 0.2*RGCspacing
            synapsesPerRGC = round(synapsesPerCone * desiredConesToRGCratios(k));
            thisRGCSynapseLocs = bsxfun(@plus, RGCRFPositionsMicrons(k,:), randn(synapsesPerRGC,2) * 0.2*RGCRFSpacingsMicrons(k));
            
            % Find closest cone synapse to each of the RGC synapses
            [~,indicesOfConnectedConeSynapses] = pdist2(allConeSynapseLocs, thisRGCSynapseLocs,   'euclidean', 'Smallest', 1);
        
            % Find the cone to which each of the closest synapse belongs to
            coneIDs = floor((indicesOfConnectedConeSynapses-1)/synapsesPerCone)+1;
            fprintf('RGC % d connects to %d cones\n', k, numel(unique(coneIDs)));
       
            % Form connectivity matrix
            for ll = 1:numel(coneIDs)
                connectivityMatrix(k, coneIDs(ll)) = connectivityMatrix(k, coneIDs(ll)) + 1;
            end
        
            % Save cone synapses connected to this RGC
            RGCconnections{k} = indicesOfConnectedConeSynapses;
            %allConeSynapseLocs(indicesOfConnectedConeSynapses,:) = Inf;
        end  % OLD connection approach
        
    end
    
    % Total connection weighgts to each RGC is 1
    connectivityMatrix = bsxfun(@times, connectivityMatrix,  1./ sum(connectivityMatrix,2));
    
    for k = 1:numel(idxRGC)
        x = RGCRFPositionsMicrons(k,1);
        y = RGCRFPositionsMicrons(k,2);
        text(x,y,sprintf('%2.2f', desiredConesToRGCratios(k)));
    end
    
    % Instantiate a plotlab object
    plotlabOBJ = plotlab();
    
    % Apply the default plotlab recipe overriding 
    % the color order and the figure size
    plotlabOBJ.applyRecipe(...
        'colorOrder', [0 0 0; 1 0 0.5], ...
        'figureWidthInches', 12, ...
        'figureHeightInches', 12);
    
    % New figure
    hFig = figure(1); clf; hold on;
    
    % Cones
    scatter(conePositionsMicrons(:,1), conePositionsMicrons(:,2));
    scatter(RGCRFPositionsMicrons(:,1), RGCRFPositionsMicrons(:,2));
    for k = 1:numel(idxCones)
        scatter(conePositionsMicrons(k,1), conePositionsMicrons(k,2), 'b');
       % coneSynapseIndices = (k-1)*synapsesPerCone+(1:synapsesPerCone);
       % plot(allConeSynapseLocs(coneSynapseIndices,1), allConeSynapseLocs(coneSynapseIndices,2), 'k.');
    end
    
    plot(allConeSynapseLocsOriginal(:,1), allConeSynapseLocsOriginal(:,2), 'k.');
    for k = 1:numel(idxRGC)
        
        x = RGCRFPositionsMicrons(k,1);
        y = RGCRFPositionsMicrons(k,2);
        scatter(x,y, 'g');
        pause
        indicesOfConnectedConeSynapses = RGCconnections{k};
        scatter(allConeSynapseLocsOriginal(indicesOfConnectedConeSynapses,1), allConeSynapseLocsOriginal(indicesOfConnectedConeSynapses,2), 50);
        
        text(x,y,sprintf('%2.2f', desiredConesToRGCratios(k)));
        pause
    end
    title(sprintf('cones/RGC ratio in patch: %2.2f', numel(idxCones)/numel(idxRGC)));
    axis 'square'
    set(gca, 'XLim', xo+1.1*radius*[-1 1], 'YLim', yo+1.1*radius*[-1 1]);
    
    connectivityMatrix
    figure(3);
    scatter(conePositionsMicrons(:,1), conePositionsMicrons(:,2)); hold on
    scatter(RGCRFPositionsMicrons(:,1), RGCRFPositionsMicrons(:,2), 200);
    for RGCindex = 1:numel(idxRGC)
        for coneIndex = 1:numel(idxCones)
            if (connectivityMatrix(RGCindex, coneIndex)>0)
            plot([RGCRFPositionsMicrons(RGCindex,1) conePositionsMicrons(coneIndex,1)], ...
                 [RGCRFPositionsMicrons(RGCindex,2) conePositionsMicrons(coneIndex,2)], ...
                 'r-', 'LineWidth', 4*connectivityMatrix(RGCindex, coneIndex));
            end
        end
        
    end
    
end

function coneSpacings = coneStats(conePositions)
    p = pdist2(conePositions, conePositions, 'euclidean', 'Smallest', 6);
    p = p(2:end,:);
    coneSpacings = mean(p,1);
end
