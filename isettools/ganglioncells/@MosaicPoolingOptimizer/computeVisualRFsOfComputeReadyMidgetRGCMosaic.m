function computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            maxSFcyclesPerDegree, stimSizeDegs, stimPositionDegs, ...
            coneMosaicSubspaceResponsesFileName, ...
            mRGCMosaicSubspaceResponsesFileName, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            visualizeOptimallyMappedRFmapLocations, ...
            varargin)

   
    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('visualizedResponses', false, @islogical);
    
    p.parse(varargin{:});
    parPoolSize = p.Results.parPoolSize;
    visualizedResponses = p.Results.visualizedResponses;

    optimallyMappedVisualRFmaps = [];
    totalOptimallyMappedRFs = 0;
    indicesOfOptimallyMappedRGCs = [];


    if (isempty(stimPositionDegs))
        stimPositionDegs = theComputeReadyMRGCmosaic.eccentricityDegs;
    end


    if (reComputeInputConeMosaicSubspaceRFmappingResponses || reComputeMRGCMosaicSubspaceRFmappingResponses || reComputeRFs)

        % Compute the RF maps for ALL cells using stimuli at this 
        % position. Only some cells will be optimally mapped at each
        % position. These optimally derived RF maps are extracted by
        % the optimalyMappedRFsAtThisGridPosition() function
        computeRFmapsForAllCellsUsingStimuliAtTargetPosition(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            maxSFcyclesPerDegree, stimSizeDegs, stimPositionDegs, ...
            parPoolSize, ...
            coneMosaicSubspaceResponsesFileName, ...
            mRGCMosaicSubspaceResponsesFileName, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            visualizedResponses);

    end

    if (visualizeOptimallyMappedRFmapLocations)
        visualizeAllOptimallyMappedRFmapLocations(optimallyMappedSubspaceRFmapsFileName, stimPositionDegs);
    end

end

function visualizeAllOptimallyMappedRFmapLocations(optimallyMappedSubspaceRFmapsFileName, stimPositionDegs)
    % Save all the optimally mapped visual RF maps

    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f.mat', stimPositionDegs(1), stimPositionDegs(2));
    optimallyMappedSubspaceRFmapsFileName = strrep(optimallyMappedSubspaceRFmapsFileName, '.mat', positionPostFix);


    fprintf('Loading optimally mapped subspace RF maps from %s\n', optimallyMappedSubspaceRFmapsFileName);
    load(optimallyMappedSubspaceRFmapsFileName, 'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCsAtThisPosition');

    figure(22); clf
    for iCell = 1:numel(indicesOfOptimallyMappedRGCsAtThisPosition)
        theRGCindex = indicesOfOptimallyMappedRGCsAtThisPosition(iCell);
        d = optimallyMappedVisualRFmaps{iCell};
        imagesc(d.spatialSupportDegsX,  d.spatialSupportDegsY, d.theRFmap);
        axis 'image'
        set(gca, 'CLim', 50*[-1 1]);
        colormap(gray(1024));
        title(sprintf('RGC %d: maxRF = %f', theRGCindex, max(d.theRFmap(:))));
        drawnow;
    end

end


function computeRFmapsForAllCellsUsingStimuliAtTargetPosition( ...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            maxSFcyclesPerDegree, stimSizeDegs, stimPositionDegs, ...
            parPoolSize,...
            coneMosaicSubspaceResponsesFileName, ...
            mRGCMosaicSubspaceResponsesFileName, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            visualizedResponses)

    % Encode examined spatial position
    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f.mat', stimPositionDegs(1), stimPositionDegs(2));
    coneMosaicSubspaceResponsesFileName = strrep(coneMosaicSubspaceResponsesFileName, '.mat', positionPostFix);
    mRGCMosaicSubspaceResponsesFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', positionPostFix);
    optimallyMappedSubspaceRFmapsFileName = strrep(optimallyMappedSubspaceRFmapsFileName, '.mat', positionPostFix);

    % Compute responses of the input cone mosaic to the subspace RF mapping stimuli
    if (reComputeInputConeMosaicSubspaceRFmappingResponses)
        MosaicPoolingOptimizer.generateInputConeMosaicSubspaceRFmappingLinearResponses(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            maxSFcyclesPerDegree, stimSizeDegs, stimPositionDegs, ...
            coneMosaicSubspaceResponsesFileName, ...
            'visualizedResponses', visualizedResponses, ...
            'parPoolSize', parPoolSize);
    end

    fprintf('\nLoading cone mosaic subspace responses and Hartley spatial modulation patterns ...');
    % Load the previously computed responses
    load(coneMosaicSubspaceResponsesFileName, ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceLinearResponses', 'theConeMosaicNullResponses');
    HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);
    fprintf('Done loading !\n');
    

    [HartleyStimNum, nCones] = size(theConeMosaicSubspaceLinearResponses);

    % Compute mRGC mosaic responses to theConeMosaicSubspaceResponses
    
    if (reComputeMRGCMosaicSubspaceRFmappingResponses)
        fprintf('MRGC mosaic subspace RF maps and responses will be saved to %s \n', mRGCMosaicSubspaceResponsesFileName);
        fprintf('Computing visual subspace RF mapping responses for all RGCs in the mosaic ... \n');
        nTimePoints = 1;
        nTrials = 1;
        theConeMosaicResponseTemporalSupportSeconds = [0];
        
        theMRGCMosaicSubspaceRFmappingLinearResponses = zeros(...
                     HartleyStimNum, theComputeReadyMRGCmosaic.rgcsNum, 'single');

        % Use all processors
        [shutdownParPoolOnceCompleted, numWorkers] = MosaicPoolingOptimizer.resetParPool([]);
    
        % Noise-free responses
        theComputeReadyMRGCmosaic.noiseFlag = 'none';

        parfor iStim = 1:HartleyStimNum
             fprintf('Computing mRGC mosaic response for Hartley pattern %d of %d (using %d parallel processes).\n', ...
                 iStim, HartleyStimNum, numWorkers);
             theConeMosaicResponse = squeeze(theConeMosaicSubspaceLinearResponses(iStim,:));
             theConeMosaicResponse = reshape(theConeMosaicResponse, [nTrials nTimePoints nCones]);

             % Compute !
             [theMRGCMosaicResponse, ~, theMRGCresponseTemporalSupportSeconds] = ...
                    theComputeReadyMRGCmosaic.compute(theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicSubspaceRFmappingLinearResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(nTrials, nTimePoints ,:)));
        end
    
        fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses to %s ...', mRGCMosaicSubspaceResponsesFileName);
        save(mRGCMosaicSubspaceResponsesFileName, ...
            'theMRGCMosaicSubspaceRFmappingLinearResponses', ...
            'spatialSupportDegs', 'lIndices', 'mIndices', ...
            '-v7.3');
    end

    if (reComputeRFs)
        % Load theMRGCMosaicSubspaceRFmappingLinearResponses
        load(mRGCMosaicSubspaceResponsesFileName, ...
                'theMRGCMosaicSubspaceRFmappingLinearResponses');
    
        % Compute RF maps of all cells in the MRGC mosaic for stimuli delivered
        % in this position. Only some cells are optimally mapped in each
        % position grid, as the subspace stim size is usually small not covering the
        % entire RGC mosaic. These are selected by the
        % optimalyMappedRFsAtThisGridPosition() function later on
        theMRGCMosaicVisualRFmaps = computeRFs(...
            theMRGCMosaicSubspaceRFmappingLinearResponses, ...
            HartleySpatialModulationPatterns);
    
        fprintf('\nSaving computed visual RFs to %s ...', mRGCMosaicSubspaceResponsesFileName);
        save(mRGCMosaicSubspaceResponsesFileName, ...
            'theMRGCMosaicVisualRFmaps', '-append');
        fprintf('Done saving! \n');

        % Extract visual RF maps for cells that are optimally mapped at this grid position
        optimalyMappedRFsAtThisGridPosition(...
            theComputeReadyMRGCmosaic, ...
            stimPositionDegs,  mRGCMosaicSubspaceResponsesFileName, optimallyMappedSubspaceRFmapsFileName);
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
            allResponses = squeeze(theSubspaceRFmappingLinearResponses(:,iCell));
            for iStim = 1:nStim
                theRFmap = theRFmap + ...
                        single(squeeze(HartleySpatialModulationPatterns(iStim,:,:)) * allResponses(iStim));
            end
            theRFmaps{iCell} = theRFmap;
        end
    end
end


function optimalyMappedRFsAtThisGridPosition(theComputeReadyMRGCmosaic, ...
    stimPositionDegs, mRGCMosaicSubspaceResponsesFileName, optimallyMappedSubspaceRFmapsFileName)

    % Encode examined spatial position
    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f.mat', stimPositionDegs(1), stimPositionDegs(2));
    optimallyMappedSubspaceRFmapsFileName = strrep(optimallyMappedSubspaceRFmapsFileName, '.mat', positionPostFix);

    % Check to see if the subspace responses have been generated at this
    % grid position
    fileExists = isfile(mRGCMosaicSubspaceResponsesFileName);
    if (~fileExists)
        error('File %s does not exist.\n Not extracting optimally mapped visual RF maps at (x,y)=(%2.1f,%2.1f)\n', ...
            mRGCMosaicSubspaceResponsesFileName, ...
            stimPositionDegs(1), stimPositionDegs(2));
    end

    % Check to see if subspace RF maps have been computed at this grid position
    s = whos('-file', mRGCMosaicSubspaceResponsesFileName);
    computedVariableNames = cell(1, numel(s));
    for i = 1:numel(s)
        computedVariableNames{i} = s(i).name;
    end
    
    if (~ismember('theMRGCMosaicVisualRFmaps', computedVariableNames))
        error('Subspace responses are computed, however RF maps have not been computed yet.\nNot extracting optimally mapped visual RF maps at (x,y) = (%2.1f,%2.1f)\n', ...
            stimPositionDegs(1), stimPositionDegs(2));
    end
    
    fprintf('Loading RF map data. Please wait ...\n');
    % All good. Extract the computed subspace RF maps
    load(mRGCMosaicSubspaceResponsesFileName, ...
        'spatialSupportDegs', ...
        'theMRGCMosaicVisualRFmaps');

    % All the cells in the MRGC mosaic
    cellsNum = theComputeReadyMRGCmosaic.rgcsNum;
    
    % Extract the optimally subspace RF maps
    fprintf('\nExtracting optimally mapped RGCs for stimulus position ...\n\t')
    counter = 0;
    rfAmplitudes = nan(1, cellsNum);

    for iCell = 1:cellsNum
        if (mod(iCell,50) == 1)
            fprintf('.');
            counter = counter + 1;
        end
        if (counter == 20)
            counter = 0;
            fprintf('\n\t');
        end

        theRFmap = theMRGCMosaicVisualRFmaps{iCell};
        rfAmplitudes(iCell) = max(abs(theRFmap(:)));
    end

    maxRFamplitude = max(rfAmplitudes(:));
    indicesOfOptimallyMappedRGCsAtThisPosition = find(rfAmplitudes > 0.1*maxRFamplitude);
    optimallyMappedRFsForThisGridPosition = numel(indicesOfOptimallyMappedRGCsAtThisPosition);
    optimallyMappedVisualRFmaps = cell(1, optimallyMappedRFsForThisGridPosition);

    for iCell = 1:optimallyMappedRFsForThisGridPosition
        theRGCindex = indicesOfOptimallyMappedRGCsAtThisPosition(iCell);
        optimallyMappedVisualRFmaps{iCell} = struct(...
            'theRFmap', theMRGCMosaicVisualRFmaps{theRGCindex}, ...
            'spatialSupportDegsX', spatialSupportDegs+stimPositionDegs(1), ...
            'spatialSupportDegsY', spatialSupportDegs+stimPositionDegs(2) ...
            );
    end

    fprintf('\nThere were %d/%d optimally mapped RF maps at (x,y) = %2.1f,%2.1f\n', optimallyMappedRFsForThisGridPosition, cellsNum, ...
        stimPositionDegs(1), stimPositionDegs(2));

    % Save all the optimally mapped visual RF maps
    fprintf('Saving optimally mapped subspace RF maps to %s\n', optimallyMappedSubspaceRFmapsFileName);
    save(optimallyMappedSubspaceRFmapsFileName, 'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCsAtThisPosition', '-v7.3');

end