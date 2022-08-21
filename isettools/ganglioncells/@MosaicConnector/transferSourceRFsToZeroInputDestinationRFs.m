function transferSourceRFsToZeroInputDestinationRFs(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;

    % Compute the # of cone inputs for all destinationRFs
    inputsNumToAllDestinationRFs = squeeze(sum(obj.connectivityMatrix,1));

    % Find out how many destination RFs have zero-inputs
    indicesOfZeroInputDestinationRFs = find(inputsNumToAllDestinationRFs == 0);
    fprintf('There are %d destination RFs with zero inputs.\n', numel(indicesOfZeroInputDestinationRFs));

    if (isempty(indicesOfZeroInputDestinationRFs))
        return;
    end

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
                inputsNumToAllDestinationRFs, videoOBJ);

        end % if (~isempty(indicesOfZeroInputDestinationRFs))
    end % iGroup

    if (generateProgressVideo)
        videoOBJ.close();
    end
end


function indicesOfZeroInputDestinationRFs = attemptToTrasfterInputsToRemainingZeroInputDestinationRFs(obj, ...
    indicesOfDestinationRFs, indicesOfZeroInputDestinationRFs, ...
    inputsNumToAllDestinationRFs, videoOBJ)

    for iRGC = 1:numel(indicesOfDestinationRFs)
        if (~isempty(indicesOfZeroInputDestinationRFs))
            % The multi-input destination RF 
            theMultiInputDestinationRFindex = indicesOfDestinationRFs(iRGC);
            
            % The inputs of the multi-input destinationRF 
            theMultiInputDestinationRFinputIndices = find(squeeze(obj.connectivityMatrix(:,theMultiInputDestinationRFindex))>0);
            theInputNumerosity = numel(theMultiInputDestinationRFinputIndices);
            
            % Must be at least 2 inputs
            if (theInputNumerosity < 2)
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


    % Update the destinationRF spacings based on the updated connectivity
    obj.updateDestinationRFspacingsBasedOnCentroids();

    % Visualize connectivity
    if (obj.visualizeConnectivityAtIntermediateStages)
        obj.visualizeCurrentConnectivity(1004);
    end
end
