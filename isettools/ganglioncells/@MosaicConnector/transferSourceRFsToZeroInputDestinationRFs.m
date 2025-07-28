function transferSourceRFsToZeroInputDestinationRFs(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.addParameter('debugPhase', false, @islogical);
    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;
    debugPhase = p.Results.debugPhase;
    
     % Video setup
    if (generateProgressVideo)
        videoOBJ = VideoWriter('Step4', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    else
        videoOBJ = [];
    end

    for pass = 1:10    
        % Zone destinationRFs based on their input numerosities
        eccZones = []; inputNumerosityGroups = [];
        [eccZones, inputNumerosityGroups, indicesOfZeroInputDestinationRFs, ...
         inputsNumToAllDestinationRFs, indicesOfDestinationRFs, ...
         eccForAllDestinationRFs, eccZoneForAllDestinationRFs] = zoneDestinationRFs(obj, eccZones, inputNumerosityGroups);

        remainingZeroInputDestinationRFsNum = numel(indicesOfZeroInputDestinationRFs);
    
        if (remainingZeroInputDestinationRFsNum == 0)
            fprintf('There are no more zero input RFs for pass %d\n', pass);
            continue;
        end

        fprintf('---> Transfer to zero RGCs - PASS %d with %d remaining zero-input destination RFs\n', ...
            pass, numel(indicesOfZeroInputDestinationRFs));

        % For each zone, treat destination RFs based on their inputNumerosityGroup
        % beginning with highest numerosity to lower
        for currentEccZone = 1:numel(eccZones)

            indicesOfDestinationRFsWithSuccesfulTransferFromLastGroup = [];

            for iGroup = 1:numel(inputNumerosityGroups)
                if (inputNumerosityGroups(iGroup) < 2)
                    % We dont do anything with a 1-input destination RF
                    continue
                end

                remainingZeroInputDestinationRFsNum = numel(indicesOfZeroInputDestinationRFs);
                if (remainingZeroInputDestinationRFsNum == 0)
                    fprintf('There are no more zero input RFs to deal with any %d-input RFs in zone %d\n', ...
                        inputNumerosityGroups(iGroup), ...
                        currentEccZone);
                    continue;
                end
    
                % All destination RFs that have the same # of inputs and are in current zone
                idx = find(...
                        (inputsNumToAllDestinationRFs(:) == inputNumerosityGroups(iGroup)) & ...
                        (eccZoneForAllDestinationRFs(:) == currentEccZone));

                if ( ...
                        (numel(idx)>0) || ...
                        ((~isempty(indicesOfDestinationRFsWithSuccesfulTransferFromLastGroup)) && (numel(indicesOfDestinationRFsWithSuccesfulTransferFromLastGroup)>0))...
                    )

                    %fprintf('Breaking down %d %d-input RFs using %d remaining zero input RFs\n', ...
                    % numel(indicesOfDestinationRFsInThisGroup), inputNumerosityGroups(iGroup), remainingZeroInputDestinationRFsNum);
    
                    % Transfer inputs to any remaining zero input destination RFs
                    deltaEcc =  eccZones(2)-eccZones(1);
                    phaseDescriptor = sprintf('transfer to zero input destination RFs (zone %d): %d %d-input RFs', ...
                        currentEccZone, ...
                        numel(idx),...
                        inputNumerosityGroups(iGroup));
    
                    if (debugPhase)
                        figNo = 998;
                        obj.visualizeCurrentConnectivity(figNo, ...
                            'titleString', sprintf('BEFORE %s', phaseDescriptor), ...
                            'visualizedDestinationRFindices', indicesOfDestinationRFs(idx), ...
                            'onlyVisualizeConnectivity', true);
                    end

                    indicesOfDestinationRFsToWorkOn = indicesOfDestinationRFs(idx);
                    if (~isempty(indicesOfDestinationRFsWithSuccesfulTransferFromLastGroup)) && (numel(indicesOfDestinationRFsWithSuccesfulTransferFromLastGroup)>0)
                        fprintf('Adding %d destinationRFs that were split during the last transfer to zero input destination RFs\n', numel(indicesOfDestinationRFsWithSuccesfulTransferFromLastGroup));
                        indicesOfDestinationRFsToWorkOn = cat(1, indicesOfDestinationRFsToWorkOn(:), indicesOfDestinationRFsWithSuccesfulTransferFromLastGroup(:));
                    end

                    if (debugPhase)
                        figNo = 999;
                        obj.visualizeCurrentConnectivity(figNo, ...
                            'titleString', sprintf('BEFORE %s', phaseDescriptor), ...
                            'visualizedDestinationRFindices', indicesOfDestinationRFsToWorkOn, ...
                            'onlyVisualizeConnectivity', true);
                    end

                    [indicesOfZeroInputDestinationRFs, indicesOfDestinationRFsWithSuccesfulTransferFromLastGroup] = ...
                        attemptToTrasfterInputsToRemainingZeroInputDestinationRFs(obj, ...
                            indicesOfDestinationRFsToWorkOn, indicesOfZeroInputDestinationRFs, ...
                            inputNumerosityGroups(iGroup), ...
                            phaseDescriptor, currentEccZone, videoOBJ);
    
                    if (debugPhase)
                        figNo = 1000;
                        obj.visualizeCurrentConnectivity(figNo, ...
                            'titleString', sprintf('AFTER %s', phaseDescriptor), ...
                            'visualizedDestinationRFindices', indicesOfDestinationRFs(idx), ...
                            'onlyVisualizeConnectivity', true);
                    end

                    fprintf('Used %d of %d zero input RFs in zone %d (%2.2f-%2.2f)to break down %d RFs with %d inputs\n', ...
                        remainingZeroInputDestinationRFsNum-numel(indicesOfZeroInputDestinationRFs), ...
                        remainingZeroInputDestinationRFsNum, currentEccZone, min(eccForAllDestinationRFs(idx(:))), max(eccForAllDestinationRFs(idx(:))),...
                        numel(idx), inputNumerosityGroups(iGroup));
               end %  if (numel(indicesOfDestinationRFsInThisGroup)>0)
    
            end  % for currentEccZone = 1:numel(eccZonesNum)
        end % iGroup
    end % Pass
    
    if (generateProgressVideo)
        videoOBJ.close();
    end

    % If there are zero-input destination RFs still available, remove them
    % Find out how many destination RFs have zero-inputs
    
    inputsNumToAllDestinationRFs = squeeze(sum(obj.connectivityMatrix,1));
    indicesOfZeroInputDestinationRFs = find(inputsNumToAllDestinationRFs == 0);
    if ~isempty(indicesOfZeroInputDestinationRFs)
        fprintf('There are STILL %d destination RFs with zero inputs. Will remove them.\n', numel(indicesOfZeroInputDestinationRFs));
        obj.removeZeroInputDestinationRFs(indicesOfZeroInputDestinationRFs);
    end
    
    % Visualize connectivity
    if (obj.visualizeConnectivityAtIntermediateStages)
        obj.intermediateFigureHandles{numel(obj.intermediateFigureHandles)+1} = ...
            obj.visualizeCurrentConnectivity(1004);
    end
end


function [indicesOfZeroInputDestinationRFs, indicesOfDestinationRFsWithSuccesfulTransfer] = attemptToTrasfterInputsToRemainingZeroInputDestinationRFs(obj, ...
    indicesOfDestinationRFs, indicesOfZeroInputDestinationRFs, ...
    inputNumerosity, phaseDescriptor, currentEccZone, videoOBJ)

    % Feedback
    fprintf('Zone %d: Transfering inputs from %d input destination RFs to zero input destination RFs. Please wait ...\n', currentEccZone, inputNumerosity);

    indicesOfDestinationRFsWithSuccesfulTransfer = [];
    for iRGC = 1:numel(indicesOfDestinationRFs)
        
        fprintf('.');
        if (mod((iRGC-1),100) == 0)
           fprintf('  [%d/%d]\n', iRGC , numel(indicesOfDestinationRFs));
        end

        if (~isempty(indicesOfZeroInputDestinationRFs))
            % The multi-input destination RF 
            theMultiInputDestinationRFindex = indicesOfDestinationRFs(iRGC);
            
            % The inputs of the multi-input destinationRF 
            theMultiInputDestinationRFinputIndices = find(squeeze(obj.connectivityMatrix(:,theMultiInputDestinationRFindex))>0);
            theInputNumerosity = numel(theMultiInputDestinationRFinputIndices);
            
            % Must be at least 2 inputs
            if (theInputNumerosity < 2)
                theInputNumerosity
                error('How can this be? Less than 2 inputs. No transfer to zero input')
                continue;
            end
    
            % Find zero-input destination RF to use
            theZeroInputDestinationRF = indicesOfZeroInputDestinationRFs(1);

            % Remove it from the pool of zero-input destination RFs
            indicesOfZeroInputDestinationRFs = indicesOfZeroInputDestinationRFs(2:end);
    
            % Transfer inputs from theMultiInputDestinationRFindex to theZeroInputDestinationRF
            optimizeTransferOfInputRFsToZeroInputDestinationRF(obj,...
                 theMultiInputDestinationRFindex, theMultiInputDestinationRFinputIndices, ...
                 theZeroInputDestinationRF);
                 
            % Add theMultiInputDestinationRFindex to the pool of destinationRFs with successful transfer
            indicesOfDestinationRFsWithSuccesfulTransfer(numel(indicesOfDestinationRFsWithSuccesfulTransfer)+1) = ...
                       theMultiInputDestinationRFindex;

            if (~isempty(videoOBJ))
                % Visualize current connectivity
                hFig = obj.visualizeCurrentConnectivity(1004);
                videoOBJ.writeVideo(getframe(hFig));
            end
        end % (~isempty(indicesOfZeroInputDestinationRFs))

    end % iRGC

    % Feedback
    fprintf('  [%d/%d]\n', iRGC, numel(indicesOfDestinationRFs));
       
    % Update the destinationRF spacings based on the updated connectivity
    obj.updateDestinationRFspacingsBasedOnCentroids();

    % Save the metaDataStuct for this stage
    if (obj.saveIntermediateConnectivityStagesMetaData)
        obj.updateIntermediateMetaDataStructs(phaseDescriptor, [], []);
    end
    
    % Visualize connectivity at this stage
    if (obj.visualizeConnectivityAtIntermediateStages)
        obj.visualizeCurrentConnectivity(1004);
    end
end


function [eccZones, inputNumerosityGroups, indicesOfZeroInputDestinationRFs, ...
         inputsNumToAllDestinationRFs, indicesOfDestinationRFs, ...
         eccForAllDestinationRFs, eccZoneForAllDestinationRFs] = zoneDestinationRFs(obj, eccZones, inputNumerosityGroups)

    % Compute the # of cone inputs for all destinationRFs
    inputsNumToAllDestinationRFs = squeeze(full(sum(obj.connectivityMatrix,1)));
    rgcsNum = size(obj.connectivityMatrix,2);

    % Find out how many destination RFs have zero-inputs
    indicesOfZeroInputDestinationRFs = find(inputsNumToAllDestinationRFs == 0);
    fprintf('There are %d destination RFs with zero inputs.\n', numel(indicesOfZeroInputDestinationRFs));
    if (isempty(indicesOfZeroInputDestinationRFs))
         inputsNumToAllDestinationRFs = [];
         indicesOfDestinationRFs = [];
         eccForAllDestinationRFs = [] ;
         eccZoneForAllDestinationRFs = [];
        return;
    end

    if (isempty(eccZones))
        % Sort the zero input destinationRFs according to their position in the destinationRF lattice
        idx = sortZeroInputDestinationRFs(obj,indicesOfZeroInputDestinationRFs);
        indicesOfZeroInputDestinationRFs = indicesOfZeroInputDestinationRFs(idx);
    end

    indicesOfDestinationRFs = setdiff(1:rgcsNum, indicesOfZeroInputDestinationRFs);
    inputsNumToAllDestinationRFs = inputsNumToAllDestinationRFs(indicesOfDestinationRFs);

    % Sort the non-zero input destinationRFs according to their position in the destinationRF lattice
    [idx, eccForAllDestinationRFs] = obj.sortDestinationRFsBasedOnOptimizationCenter(...
        indicesOfDestinationRFs, ...
        'ignoreInfCentroids', true, ...
        'averageEcc', ~true);

    inputsNumToAllDestinationRFs = inputsNumToAllDestinationRFs(idx);
    indicesOfDestinationRFs = indicesOfDestinationRFs(idx);

    % Zone destination RFs into eccZonesNum
    if (isempty(eccZones))
        eccZonesNum = 10;
        minEcc = min(eccForAllDestinationRFs(:));
        maxEcc = max(eccForAllDestinationRFs(:));
        eccZones = linspace(minEcc, maxEcc,eccZonesNum);
    end

    deltaEcc = eccZones(2)-eccZones(1);
    eccZoneForAllDestinationRFs = zeros(1, numel(eccForAllDestinationRFs));
    for iZone = 1:numel(eccZones)
        idx = find(...
            (eccForAllDestinationRFs(:) >= eccZones(iZone)-deltaEcc/2) & ...
            (eccForAllDestinationRFs(:) < eccZones(iZone)+deltaEcc));
        eccZoneForAllDestinationRFs(idx) = iZone;
        fprintf('Assigned %d RFs to zone %d.\n', numel(idx), iZone);
    end

    if (isempty(inputNumerosityGroups))
        % Group destination RFs based on their # of inputs
        inputNumerosityGroups = unique(inputsNumToAllDestinationRFs);
        inputNumerosityGroups = sort(inputNumerosityGroups, 'descend');
    end

end


function sortedIndices = sortZeroInputDestinationRFs(obj,unsortedIndices)

    centroidsOfDestinationRFsInThisGroup = ...
        obj.destinationLattice.RFpositionsMicrons(unsortedIndices,:);

    if (isempty(obj.sourceLatticeCenter))
        obj.sourceLatticeCenter = mean(obj.sourceLattice.RFpositionsMicrons,1);
    end

    diff = bsxfun(@minus, centroidsOfDestinationRFsInThisGroup, obj.sourceLatticeCenter);
    ecc = sqrt(sum(diff.^2,2));
    
    % Compute sorted indices of destination RFs in increasing eccentricity
    [~, sortedIndices] = sort(ecc, 'ascend');
end

function optimizeTransferOfInputRFsToZeroInputDestinationRF(obj,...
             theMultiInputDestinationRFindex, theMultiInputDestinationRFinputIndices, ...
             theZeroInputDestinationRF)

    projectedCosts = zeros(numel(theMultiInputDestinationRFinputIndices),1);

    parfor idx = 1:numel(theMultiInputDestinationRFinputIndices)
        theSourceRFindex = theMultiInputDestinationRFinputIndices(idx);

        % Remove theInputRFindices from theMultiInputDestinationRFinputIndices
        theRemainingSourceRFindices = setdiff(theMultiInputDestinationRFinputIndices, theSourceRFindex);
       
        projectedCosts(idx) = computeCost(...
            obj.sourceLattice.metaData.coneTypes(theRemainingSourceRFindices), ...
            obj.sourceLattice.RFpositionsMicrons(theRemainingSourceRFindices,:), ...
            obj.wiringParams.spatialChromaticUniformityTradeoff);
    end

    % Find the sourceRFindex that minimizes the total cost
    [~,iBestSourceRFindex] = min(projectedCosts(:));
    
    % The indices of the inputRFs to be transfered
    theSourceRFindexToBeTransferredToZeroInputDestinationRF = theMultiInputDestinationRFinputIndices(iBestSourceRFindex);

    % DISCONNECT theInputRFindices from theMultiInputDestinationRFindex
    obj.connectivityMatrix(theSourceRFindexToBeTransferredToZeroInputDestinationRF, theMultiInputDestinationRFindex) = 0;

    % And CONNECT the to theZeroInputDestinationRF
    obj.connectivityMatrix(theSourceRFindexToBeTransferredToZeroInputDestinationRF, theZeroInputDestinationRF) = 1;
    
    % Update the centroids of the 2 destination RFs
    destinationRFList = [theMultiInputDestinationRFindex theZeroInputDestinationRF];
    obj.updateDestinationCentroidsFromInputs(destinationRFList);
end


function theCost = computeCost(sourceRFtypes, sourceRFpositions, spatialChromaticUniformityTradeoff)

    theSourceTypeUniformityCost = coneToMidgetRGCConnector.spectralUniformityCost(sourceRFtypes, []);
    theSpatialUniformityCost = sqrt(sum(var(sourceRFpositions,0,1)));

    theCost = spatialChromaticUniformityTradeoff * theSpatialUniformityCost + ...
            (1-spatialChromaticUniformityTradeoff) * theSourceTypeUniformityCost;

end
