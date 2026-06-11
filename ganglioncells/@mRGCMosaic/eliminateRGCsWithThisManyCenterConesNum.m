function eliminateRGCsWithThisManyCenterConesNum(obj, targetCenterConesNum)

   if (isempty(obj.rgcRFcenterConePoolingMatrix))
       indicesOfRGCsToEliminate = obj.indicesOfRGCsWithThisManyCenterCones(targetCenterConesNum);
       fprintf('Will eliminate %d RGCs which all have %d cones in their RF center', numel(indicesOfRGCsToEliminate), targetCenterConesNum)
   else
       error('RGC elimination can only be done before cones are connected to surrounds\n');
   end
   

   % Visualize mosaic before RGC elimination
   identifyPooledCones = true;
   identifyInputCones = true;
   plotRFoutlines = true;
   indicesOfRGCsIdentified = indicesOfRGCsToEliminate;

   hFig = figure(10); clf;
   set(hFig, 'Position', [10 10 1800 1000], 'Color', [1 1 1], ...
        'Name', sprintf('Before elimination of RGCs with %d center cones', targetCenterConesNum));
   obj.visualize(...
            'figureHandle', hFig, ...
            'identifiedConeApertureThetaSamples', 32, ...
            'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
            'identifyPooledCones', identifyPooledCones, ...
            'identifyInputCones',identifyInputCones, ...
            'labelRGCsWithIndices', indicesOfRGCsIdentified , ...
            'plotRFoutlines', plotRFoutlines, ...
            'domainVisualizationLimits', [], ...
            'domainVisualizationTicks', [], ...
            'backgroundColor', [0 0 0], ...
            'plotTitle', sprintf('Before elimination of RGCs with %d center cones', targetCenterConesNum));


    elimininateRGCs = true;
    if (targetCenterConesNum == 2)
        % Special case if the # of center cones is 2, ask the user whether
        % to add single cone RGCs so as to split those 2-cone RGCs into 2
        % single cone RGCs
        theChoice = input('Add RGCs to split 2-input RGCs into 2 single-input RGCs [y=YES]: ', 's');
        if (strcmpi(theChoice, 'y'))
            addSingleConeRGCs = true;
        else
            addSingleConeRGCs = false;
        end

        if (addSingleConeRGCs)
            elimininateRGCs = false;
            for iRGC = 1:numel(indicesOfRGCsToEliminate)
                connectedConeIndices = find(obj.rgcRFcenterConeConnectivityMatrix(:,indicesOfRGCsToEliminate(iRGC))>0.001);
                % Add RGC
                currentRGCsNum = size(obj.rgcRFcenterConeConnectivityMatrix,2);
                addedRGCindex = currentRGCsNum + 1;
                fprintf('Adding RGC #%d\n', addedRGCindex);
                
                obj.rgcRFcenterConeConnectivityMatrix(connectedConeIndices(1), addedRGCindex) = 1;
                % Disconnect cone from 2-cone input RGC
                obj.rgcRFcenterConeConnectivityMatrix(connectedConeIndices(1), indicesOfRGCsToEliminate(iRGC)) = 0;
                
                % Update positions
                obj.rgcRFpositionsMicrons(indicesOfRGCsToEliminate(iRGC),:) = obj.inputConeMosaic.coneRFpositionsMicrons(connectedConeIndices(2),:);
                obj.rgcRFpositionsDegs(indicesOfRGCsToEliminate(iRGC),:) = obj.inputConeMosaic.coneRFpositionsDegs(connectedConeIndices(2),:);
                obj.rgcRFpositionsMicrons(addedRGCindex,:) = obj.inputConeMosaic.coneRFpositionsMicrons(connectedConeIndices(1),:);
                obj.rgcRFpositionsDegs(addedRGCindex,:) = obj.inputConeMosaic.coneRFpositionsDegs(connectedConeIndices(1),:);
                obj.rgcRFspacingsMicrons(addedRGCindex) = obj.rgcRFspacingsMicrons(indicesOfRGCsToEliminate(iRGC));
                obj.rgcRFspacingsDegs(addedRGCindex) = obj.rgcRFspacingsDegs(indicesOfRGCsToEliminate(iRGC));
                % Reset the visualization cache
                obj.visualizationCache = [];
            end
        end
    end

    if (elimininateRGCs)
        % --------- Do the RGC elimination -------
        indicesOfRGCsToKeep = setdiff(1:obj.rgcsNum, indicesOfRGCsToEliminate);
        obj.rgcRFcenterConeConnectivityMatrix = obj.rgcRFcenterConeConnectivityMatrix(:,indicesOfRGCsToKeep);
        obj.rgcRFpositionsMicrons = obj.rgcRFpositionsMicrons(indicesOfRGCsToKeep,:);
        obj.rgcRFpositionsDegs = obj.rgcRFpositionsDegs(indicesOfRGCsToKeep,:);
        obj.rgcRFspacingsMicrons = obj.rgcRFspacingsMicrons(indicesOfRGCsToKeep);
        obj.rgcRFspacingsDegs = obj.rgcRFspacingsDegs(indicesOfRGCsToKeep);
        % Reset the visualization cache
        obj.visualizationCache = [];
        % ------------------------------------------
    end


    % Visualize mosaic after RGC elimination
    indicesOfRGCsIdentified = obj.indicesOfRGCsWithThisManyCenterCones(targetCenterConesNum);
    hFig = figure(20); clf;
    set(hFig, 'Position', [10 10 1800 1000], 'Color', [1 1 1], ...
       'Name', sprintf('After elimination of RGCs with %d center cones', targetCenterConesNum));
    obj.visualize(...
            'figureHandle', hFig, ...
            'identifiedConeApertureThetaSamples', 32, ...
            'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
            'identifyPooledCones', identifyPooledCones, ...
            'identifyInputCones',identifyInputCones, ...
            'labelRGCsWithIndices', indicesOfRGCsIdentified, ...
            'plotRFoutlines', plotRFoutlines, ...
            'domainVisualizationLimits', [], ...
            'domainVisualizationTicks', [], ...
            'backgroundColor', [0 0 0], ...
            'plotTitle', sprintf('After elimination of RGCs with %d center cones', targetCenterConesNum));

    % Report stats of the updated mosaic
    obj.centerConePoolingStats();
end
