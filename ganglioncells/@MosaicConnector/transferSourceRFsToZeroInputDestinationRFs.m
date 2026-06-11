function transferSourceRFsToZeroInputDestinationRFs(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;

    % Compute the # of cone inputs for all destinationRFs
    inputsNumToAllDestinationRFs = squeeze(full(sum(obj.connectivityMatrix,1)));

    % Find out how many destination RFs have zero-inputs
    indicesOfZeroInputDestinationRFs = find(inputsNumToAllDestinationRFs == 0);
    fprintf('There are %d destination RFs with zero inputs.\n', numel(indicesOfZeroInputDestinationRFs));

    if (isempty(indicesOfZeroInputDestinationRFs))
        return;
    end

    % Sort the zero input destinationRFs according to their position in the destinationRF lattice
    idx = sortZeroInputDestinationRFsBasedOnOptimizationCenter(obj,indicesOfZeroInputDestinationRFs);
    indicesOfZeroInputDestinationRFs = indicesOfZeroInputDestinationRFs(idx);

    % Video setup
    if (generateProgressVideo)
        videoOBJ = VideoWriter('Step4', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    else
        videoOBJ = [];
    end

    % Group destination RFs based on their # of inputs
    inputNumerosityGroups = unique(inputsNumToAllDestinationRFs);
    inputNumerosityGroups = sort(inputNumerosityGroups, 'descend');
    
    % Start from most populous group going down to least populous group
    for iGroup = 1:numel(inputNumerosityGroups)

        if (inputNumerosityGroups(iGroup) < 2)
            % We dont do anything with a 1-input destination RF
            continue
        end

        if (~isempty(indicesOfZeroInputDestinationRFs))
            % All destination RFs that have the same # of inputs
            indicesOfDestinationRFsInThisGroup = ...
                find(inputsNumToAllDestinationRFs == inputNumerosityGroups(iGroup));

            % Sort destinationRF indices (based on their eccentricity &
            % optimization center)
            idx = obj.sortDestinationRFsBasedOnOptimizationCenter(indicesOfDestinationRFsInThisGroup);
            indicesOfDestinationRFsInThisGroup = indicesOfDestinationRFsInThisGroup(idx);
    
            % Transfer inputs to any remaining zero input destination RFs
            indicesOfZeroInputDestinationRFs = attemptToTrasfterInputsToRemainingZeroInputDestinationRFs(obj, ...
                indicesOfDestinationRFsInThisGroup, indicesOfZeroInputDestinationRFs, ...
                inputsNumToAllDestinationRFs, inputNumerosityGroups(iGroup), videoOBJ);

        end % if (~isempty(indicesOfZeroInputDestinationRFs))
    end % iGroup

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


function indicesOfZeroInputDestinationRFs = attemptToTrasfterInputsToRemainingZeroInputDestinationRFs(obj, ...
    indicesOfDestinationRFs, indicesOfZeroInputDestinationRFs, ...
    inputsNumToAllDestinationRFs, inputNumerosity, videoOBJ)

    % Feedback
    fprintf('Transfering inputs from %d input destination RFs to zero input destination RFs. Please wait ...', inputNumerosity);

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

            % Sanity check
            if (inputsNumToAllDestinationRFs(theMultiInputDestinationRFindex) ~= theInputNumerosity)
                error('mismatch here')
            end
    
            
            % Find zero-input destination RF to use
            theZeroInputDestinationRF = indicesOfZeroInputDestinationRFs(1);

            % Remove it from the pool of zero-input destination RFs
            indicesOfZeroInputDestinationRFs = indicesOfZeroInputDestinationRFs(2:end);
    
            % Transfer inputs from theMultiInputDestinationRFindex to theZeroInputDestinationRF
            theMultiInputDestinationRFinputWeights = full(obj.connectivityMatrix(theMultiInputDestinationRFinputIndices,theMultiInputDestinationRFindex));
            obj.optimizeTransferOfInputRFsToZeroInputDestinationRF(...
                 theMultiInputDestinationRFindex, theMultiInputDestinationRFinputIndices, ...
                 theMultiInputDestinationRFinputWeights, theZeroInputDestinationRF);
                        
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
        obj.updateIntermediateMetaDataStructs();
    end
    
    % Visualize connectivity at this stage
    if (obj.visualizeConnectivityAtIntermediateStages)
        obj.visualizeCurrentConnectivity(1004);
    end
end


function sortedIndices = sortZeroInputDestinationRFsBasedOnOptimizationCenter(obj,unsortedIndices)

    centroidsOfDestinationRFsInThisGroup = ...
        obj.destinationLattice.RFpositionsMicrons(unsortedIndices,:);

    switch (obj.wiringParams.optimizationCenter)
        case 'origin'
            ecc = sum(centroidsOfDestinationRFsInThisGroup.^2,2);
        case 'latticeCenter'
            if (isempty(obj.sourceLatticeCenter))
                obj.sourceLatticeCenter = mean(obj.sourceLattice.RFpositionsMicrons,1);
            end

            diff = bsxfun(@minus, centroidsOfDestinationRFsInThisGroup, obj.sourceLatticeCenter);
            ecc = sum(diff.^2,2);
    end % switch

   % Compute sorted indices of destination RFs in increasing eccentricity
   [~, sortedIndices] = sort(ecc, 'ascend');

end
