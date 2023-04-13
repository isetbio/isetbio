function computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, opticsParams, ...
            maxSFcyclesPerDegree, stimSizeDegs, posIncrementDegs, ...
            coneMosaicSubspaceResponsesFileName, ...
            mRGCMosaicSubspaceResponsesFileName, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            onlyVisualizeOptimallyMappedRFmaps, ...
            varargin)

   
    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('visualizedResponses', false, @islogical);
    
    p.parse(varargin{:});
    parPoolSize = p.Results.parPoolSize;
    visualizedResponses = p.Results.visualizedResponses;

    % Generate the subspace mapping grid
    [X,Y, gridNodeIndicesToCompute] = generateSubspaceRFmappingSamplingGrid(...
        theComputeReadyMRGCmosaic.inputConeMosaic.sizeDegs, ...
        theComputeReadyMRGCmosaic.inputConeMosaic.eccentricityDegs, ...
        posIncrementDegs, stimSizeDegs);

    stimXYpositionGridDegs = [X(:) Y(:)];
    optimallyMappedVisualRFmaps = [];
    totalOptimallyMappedRFs = 0;
    indicesOfOptimallyMappedRGCs = [];

    for theGridNodeIndex = 1:size(stimXYpositionGridDegs,1)

        if (~isempty(gridNodeIndicesToCompute))&&(~ismember(theGridNodeIndex, gridNodeIndicesToCompute))
            fprintf(2, 'Skipping grid node #%d\n', theGridNodeIndex);
            continue;
        end

        % stimulus position within the mRGC mosaic
        gridNodeXYpositionDegs = stimXYpositionGridDegs(theGridNodeIndex,:);
    
        if (reComputeInputConeMosaicSubspaceRFmappingResponses || ...
            reComputeMRGCMosaicSubspaceRFmappingResponses || ...
            reComputeRFs)

            % Generate and set the optics at the examined grid node position
            opticsParamsAtThisPosition = opticsParams;
            opticsParamsAtThisPosition.positionDegs = gridNodeXYpositionDegs;
            theComputeReadyMRGCmosaic.setTheOptics(opticsParamsAtThisPosition);

            % Retrieve the optics at the examined grid node position
            theOptics = theComputeReadyMRGCmosaic.theCustomOptics;

            % Compute the RF maps for ALL cells using stimuli at this grid
            % position. Only some cells will be optimally mapped at each
            % position. These optimally derived RF maps are extracted by
            % the optimalyMappedRFsAtThisGridPosition() function
            computeRFmapsForAllCellsUsingStimuliAThisGridPosition(...
                theComputeReadyMRGCmosaic, theOptics, ...
                maxSFcyclesPerDegree, stimSizeDegs, gridNodeXYpositionDegs, ...
                parPoolSize, ...
                coneMosaicSubspaceResponsesFileName, ...
                mRGCMosaicSubspaceResponsesFileName, ...
                reComputeInputConeMosaicSubspaceRFmappingResponses, ...
                reComputeMRGCMosaicSubspaceRFmappingResponses, ...
                reComputeRFs, ...
                visualizedResponses);

        end

        if (onlyVisualizeOptimallyMappedRFmaps)
            % Extract visual RF maps for cells that are optimally mapped at this grid position
            [optimallyMappedVisualRFmaps, indicesOfOptimallyRGCsAtThisPosition] = optimalyMappedRFsAtThisGridPosition(...
                optimallyMappedVisualRFmaps, theComputeReadyMRGCmosaic, ...
                stimXYpositionGridDegs, theGridNodeIndex, mRGCMosaicSubspaceResponsesFileName);
            
            if (~isempty(indicesOfOptimallyRGCsAtThisPosition))
                totalOptimallyMappedRFs = totalOptimallyMappedRFs + numel(indicesOfOptimallyRGCsAtThisPosition);
                indicesOfOptimallyMappedRGCs = cat(1,indicesOfOptimallyMappedRGCs(:), indicesOfOptimallyRGCsAtThisPosition(:));
                if (totalOptimallyMappedRFs ~= numel(unique(indicesOfOptimallyMappedRGCs)))
                    error('Multiply-optimally mapped RGCs');
                end

            end
        end
    end % for theGridNodeIndex
    
    % Save all the optimally mapped visual RF maps
    if (onlyVisualizeOptimallyMappedRFmaps)&&(~isempty(optimallyMappedVisualRFmaps))
        fprintf('Saving optimally mapped subspace RF maps to %s\n', optimallyMappedSubspaceRFmapsFileName);
        save(optimallyMappedSubspaceRFmapsFileName, 'optimallyMappedVisualRFmaps', '-v7.3');
    end

end


function [optimallyMappedVisualRFmaps, indicesOfOptimallyRGCsAtThisPosition] = optimalyMappedRFsAtThisGridPosition(...
    optimallyMappedVisualRFmaps, theComputeReadyMRGCmosaic, ...
    stimXYpositionGridDegs, theGridNodeIndex, mRGCMosaicSubspaceResponsesFileName)

    % Encode examined spatial position
    gridNodeXYpositionDegs = stimXYpositionGridDegs(theGridNodeIndex,:);
    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f.mat', gridNodeXYpositionDegs(1), gridNodeXYpositionDegs(2));
    mRGCMosaicSubspaceResponsesAtThisPositionGridFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', positionPostFix);

    % Check to see if the subspace responses have been generated at this
    % grid position
    fileExists = isfile(mRGCMosaicSubspaceResponsesAtThisPositionGridFileName);
    if (~fileExists)
        fprintf('File %s does not exist.\nNot extracting optimally mapped visual RF maps at (x,y)=(%2.1f,%2.1f)\n', ...
            mRGCMosaicSubspaceResponsesAtThisPositionGridFileName, ...
            stimXYpositionGridDegs(theGridNodeIndex,1), stimXYpositionGridDegs(theGridNodeIndex,2));
        % Return empty RGC indices
        indicesOfOptimallyRGCsAtThisPosition = [];
        return;
    end

    % Check to see if subspace RF maps have been computed at this grid
    % positions
    s = whos('-file', mRGCMosaicSubspaceResponsesAtThisPositionGridFileName);
    computedVariableNames = cell(1, numel(s));
    for i = 1:numel(s)
        computedVariableNames{i} = s(i).name;
    end
    
    if (~ismember('theMRGCMosaicVisualRFmaps', computedVariableNames))
        fprintf('Subspace responses are computed, however RF maps have not been computed yet.\nNot extracting optimally mapped visual RF maps at (x,y) = (%2.1f,%2.1f)\n', ...
            stimXYpositionGridDegs(theGridNodeIndex,1), stimXYpositionGridDegs(theGridNodeIndex,2));
        % Return empty RGC indices
        indicesOfOptimallyRGCsAtThisPosition = [];
        return;
    end
    
    fprintf('Loading data. Please wait ...\n');
    % All good. Extract the computed subspace RF maps
    load(mRGCMosaicSubspaceResponsesAtThisPositionGridFileName, ...
        'spatialSupportDegs', ...
        'theMRGCMosaicVisualRFmaps');

    % All the cells in the MRGC mosaic
    cellsNum = numel(theMRGCMosaicVisualRFmaps);

    if (isempty(optimallyMappedVisualRFmaps))
        optimallyMappedVisualRFmaps = cell(1, cellsNum);
    end

    % Extract the optimally subspace RF maps
    optimallyMappedRFsForThisGridPosition = 0;
    indicesOfOptimallyRGCsAtThisPosition = [];

    fprintf('\nExtracting optimally mapped RGCs for desired grid node ...\n\t')
    counter = 0;
    for iCell = 1:cellsNum
        if (mod(iCell,50) == 1)
            fprintf('.');
            counter = counter + 1;
        end
        if (counter == 20)
            counter = 0;
            fprintf('\n\t');
        end
        % Find the grid position that this cell is closest to
        d = sqrt(sum((bsxfun(@minus, stimXYpositionGridDegs, theComputeReadyMRGCmosaic.rgcRFpositionsDegs(iCell,:))).^2,2));
        [~, theClosestGridNodeIndex] = min(d);

        if (theClosestGridNodeIndex ~= theGridNodeIndex)
            % This grid position is not optimal for mapping this cell's RF
            continue;
        end

        indicesOfOptimallyRGCsAtThisPosition(numel(indicesOfOptimallyRGCsAtThisPosition)+1) = iCell;
        optimallyMappedRFsForThisGridPosition = optimallyMappedRFsForThisGridPosition + 1;
        optimallyMappedVisualRFmaps{iCell} = struct(...
            'theRFmap', theMRGCMosaicVisualRFmaps{iCell}, ...
            'spatialSupportDegsX', spatialSupportDegs+stimXYpositionGridDegs(theGridNodeIndex, 1), ...
            'spatialSupportDegsY', spatialSupportDegs+stimXYpositionGridDegs(theGridNodeIndex, 2) ...
            );
    end

    fprintf('\nThere were %d optimally mapped RF maps at (x,y) = %2.1f,%2.1f\n', optimallyMappedRFsForThisGridPosition, ...
        stimXYpositionGridDegs(theGridNodeIndex, 1), stimXYpositionGridDegs(theGridNodeIndex, 2));
end


function [X,Y, gridNodeIndicesToCompute] = generateSubspaceRFmappingSamplingGrid(...
    inputConeMosaicSizeDegs, inputConeMosaicEccDegs, posIncrementDegs, stimSizeDegs)

    k = round(0.5*inputConeMosaicSizeDegs(1)/posIncrementDegs)-1;
    xCoords = inputConeMosaicEccDegs(1) + ...
            (-k:1:k)*posIncrementDegs;

    k = round(0.5*inputConeMosaicSizeDegs(2)/posIncrementDegs)-1;
    yCoords = inputConeMosaicEccDegs(2) + ...
            (-k:1:k)*posIncrementDegs;
    
    [X,Y] = meshgrid(xCoords, yCoords);
    X = X(:); Y = Y(:);

    hFig = figure(2); clf;
    set(hFig, 'Name', 'Subspace RF mapping grid', 'Position', [10 10 1000 1000]);
    ax = subplot('Position', [0.1 0.1 0.85 0.85]);
    for iNode = 1:numel(X)
        if (iNode == 1) || (iNode == numel(X))
            lineWidth = 3.0;
            lineColor = [0 0 1];
        else
            lineWidth = 1.0;
            lineColor = [0.5 0.5 0.5];
        end
        rect.x = X(iNode) + 0.5*stimSizeDegs(1)*[-1 -1 1  1 -1];
        rect.y = Y(iNode) + 0.5*stimSizeDegs(1)*[-1  1 1 -1 -1];
        plot(ax,rect.x, rect.y, 'r-', 'LineWidth', lineWidth, 'Color', lineColor);
        hold(ax, 'on')
    end
    XX = sort(unique(X), 'ascend');
    dx = 0.1*(XX(2)-XX(1));
    for iNode = 1:numel(X)
        if (iNode == 1) || (iNode == numel(X))
            fontWeight = 'Bold';
            lineColor = [0 0 1];
        else
            fontWeight = 'normal';
            lineColor = [0.5 0.5 0.5];
        end
        text(ax, X(iNode)-dx, Y(iNode), sprintf('%d', iNode), 'Color', lineColor, 'FontSize', 20, 'FontWeight', fontWeight, 'BackgroundColor', [1 1 1]);
    end

    set(ax, 'TickDir', 'both', 'FontSize', 16);
    drawnow;
    gridNodeIndicesToCompute = input('Enter grid node indices to compute :' );

    for iNode = 1:numel((gridNodeIndicesToCompute))
        rect.x = X(gridNodeIndicesToCompute(iNode)) + 0.5*stimSizeDegs(1)*[-1 -1 1  1 -1];
        rect.y = Y(gridNodeIndicesToCompute(iNode)) + 0.5*stimSizeDegs(1)*[-1  1 1 -1 -1];
        plot(ax,rect.x, rect.y, 'r-', 'LineWidth', 3, 'Color', [1 0 0]);
        text(ax, X(gridNodeIndicesToCompute(iNode))-dx, Y(gridNodeIndicesToCompute(iNode)), ...
            sprintf('%d', gridNodeIndicesToCompute(iNode)), 'Color', [1 0 0], ...
            'FontSize', 20, 'FontWeight', fontWeight, 'BackgroundColor', [1 1 1]);
    end

    drawnow;


end


function computeRFmapsForAllCellsUsingStimuliAThisGridPosition( ...
            theComputeReadyMRGCmosaic, theOptics, ...
            maxSFcyclesPerDegree, stimSizeDegs, gridNodeXYpositionDegs, ...
            parPoolSize,...
            coneMosaicSubspaceResponsesFileName, ...
            mRGCMosaicSubspaceResponsesFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            visualizedResponses)

    % Encode examined spatial position
    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f.mat', gridNodeXYpositionDegs(1), gridNodeXYpositionDegs(2));
    coneMosaicSubspaceResponsesAtThisPositionGridFileName = strrep(coneMosaicSubspaceResponsesFileName, '.mat', positionPostFix);
    mRGCMosaicSubspaceResponsesAtThisPositionGridFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', positionPostFix);

    % Compute responses of the input cone mosaic to the subspace RF mapping stimuli
    if (reComputeInputConeMosaicSubspaceRFmappingResponses)
        MosaicPoolingOptimizer.generateInputConeMosaicSubspaceRFmappingLinearResponses(...
            theComputeReadyMRGCmosaic, theOptics, ...
            maxSFcyclesPerDegree, stimSizeDegs, gridNodeXYpositionDegs, ...
            coneMosaicSubspaceResponsesAtThisPositionGridFileName, ...
            'visualizedResponses', visualizedResponses, ...
            'parPoolSize', parPoolSize);
    end
    fprintf('\nLoading cone mosaic subspace responses and Hartley spatial modulation patterns ...');
    % Load the previously computed responses
    load(coneMosaicSubspaceResponsesAtThisPositionGridFileName, ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceLinearResponses');
    HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);
    fprintf('Done loading !\n');
    

    % Compute mRGC mosaic responses to theConeMosaicSubspaceResponses
    % Compute RF maps of cells in the mRGC mosaic
    [HartleyStimNum, nCones] = size(theConeMosaicSubspaceLinearResponses);
    
    if (reComputeMRGCMosaicSubspaceRFmappingResponses)
        fprintf('MRGC mosaic subspace RF maps and responses will be saved to %s \n', mRGCMosaicSubspaceResponsesAtThisPositionGridFileName);
        fprintf('Computing visual subspace RF mapping responses for all RGCs in the mosaic ... \n');
        nTimePoints = 1;
        nTrials = 1;
        theConeMosaicResponseTemporalSupportSeconds = [0];
        
        theMRGCMosaicSubspaceRFmappingLinearResponses = zeros(...
                     HartleyStimNum, theComputeReadyMRGCmosaic.rgcsNum, 'single');

        % Use all processors
        [shutdownParPoolOnceCompleted, numWorkers] = MosaicPoolingOptimizer.resetParPool([]);
    
        parfor iStim = 1:HartleyStimNum
             fprintf('Computing mRGC mosaic response for Hartley pattern %d of %d (using %d parallel processes).\n', ...
                 iStim, HartleyStimNum, numWorkers);
             theConeMosaicResponse = squeeze(theConeMosaicSubspaceLinearResponses(iStim,:));
             theConeMosaicResponse = reshape(theConeMosaicResponse, [nTrials nTimePoints nCones]);

             % Compute !
             [theMRGCMosaicResponse, theMRGCresponseTemporalSupportSeconds] = ...
                    theComputeReadyMRGCmosaic.compute(theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicSubspaceRFmappingLinearResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(1, :,:)));
        end
    
        fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses to %s ...', mRGCMosaicSubspaceResponsesAtThisPositionGridFileName);
        save(mRGCMosaicSubspaceResponsesAtThisPositionGridFileName, ...
            'theMRGCMosaicSubspaceRFmappingLinearResponses', ...
            'spatialSupportDegs', 'lIndices', 'mIndices', ...
            '-v7.3');
    end

    if (reComputeRFs)
        % Load theMRGCMosaicSubspaceRFmappingLinearResponses
        load(mRGCMosaicSubspaceResponsesAtThisPositionGridFileName, ...
                'theMRGCMosaicSubspaceRFmappingLinearResponses');
    
        % Compute RF maps of all cells in the MRGC mosaic for stimuli delivered
        % in this position grid. Only some cells are optimally mapped in each
        % position grid, as the subspace stim size is usually small not covering the
        % entire RGC mosaic. These are selected by the
        % optimalyMappedRFsAtThisGridPosition() function later on
        theMRGCMosaicVisualRFmaps = computeRFs(...
            theMRGCMosaicSubspaceRFmappingLinearResponses, ...
            HartleySpatialModulationPatterns);
    
        fprintf('\nSaving computed visual RFs to %s ...', mRGCMosaicSubspaceResponsesAtThisPositionGridFileName);
        save(mRGCMosaicSubspaceResponsesAtThisPositionGridFileName, ...
            'theMRGCMosaicVisualRFmaps', '-append');
        fprintf('Done saving! \n');
    end

end


function theRFmaps = computeRFs( ...
    theSubspaceRFmappingLinearResponses, ...
    HartleySpatialModulationPatterns)

    nStim = size(theSubspaceRFmappingLinearResponses,1);
    cellsNum = size(theSubspaceRFmappingLinearResponses,2);
    pixelsNum = size(HartleySpatialModulationPatterns,2);

   
    m = max(abs(theSubspaceRFmappingLinearResponses),[],1);
    cellsWithNonZeroResponse = find(m > 0);

    theRFmaps = cell(cellsNum, 1);
    parfor iCell = 1:cellsNum
        fprintf('Computing visual RF by accumulating Hartley patterns for RGC #%d of %d ... \n', iCell, cellsNum)
        if (ismember(iCell, cellsWithNonZeroResponse))
            theRFmap = zeros(pixelsNum, pixelsNum, 'single');
            for iStim = 1:nStim
                r = theSubspaceRFmappingLinearResponses(iStim,iCell);
                theRFmap = theRFmap + ...
                        single(squeeze(HartleySpatialModulationPatterns(iStim,:,:)) * r);
            end
            theRFmaps{iCell} = theRFmap / max(abs(theRFmap(:)));
        end
    end

end