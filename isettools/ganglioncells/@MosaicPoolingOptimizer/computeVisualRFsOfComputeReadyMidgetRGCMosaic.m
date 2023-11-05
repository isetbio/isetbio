function computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, stimulusChromaticity,...
            maxSFLimit, rfMappingPixelMagnificationFactor, ...
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

    if (isempty(stimPositionDegs))
        stimPositionDegs = theComputeReadyMRGCmosaic.eccentricityDegs;
    end


    if (reComputeInputConeMosaicSubspaceRFmappingResponses || reComputeMRGCMosaicSubspaceRFmappingResponses || reComputeRFs)

        % Compute the RF maps for ALL cells using stimuli at this 
        % position. Only some cells will be optimally mapped at each
        % position. These optimally derived RF maps are extracted by
        % the optimalyMappedRFsAtMosaicPosition() function
        computeRFmapsForAllCellsUsingStimuliAtTargetPosition(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, stimulusChromaticity, ...
            maxSFLimit, rfMappingPixelMagnificationFactor, ...
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
        visualizeAllOptimallyMappedRFmapLocations(optimallyMappedSubspaceRFmapsFileName, ...
            stimPositionDegs, rfMappingPixelMagnificationFactor);
    end

end

function visualizeAllOptimallyMappedRFmapLocations(optimallyMappedSubspaceRFmapsFileName,...
    stimPositionDegs, rfMappingPixelMagnificationFactor)
    % Save all the optimally mapped visual RF maps

    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f_PixelMagnification_%2.3f.mat', stimPositionDegs(1), stimPositionDegs(2), rfMappingPixelMagnificationFactor);
    optimallyMappedSubspaceRFmapsFileName = strrep(optimallyMappedSubspaceRFmapsFileName, '.mat', positionPostFix);

    fprintf('Loading optimally mapped subspace RF maps from %s\n', optimallyMappedSubspaceRFmapsFileName);
    load(optimallyMappedSubspaceRFmapsFileName, 'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs');


    figure(22); clf
    for iCell = 1:numel(indicesOfOptimallyMappedRGCs)
        theRGCindex = indicesOfOptimallyMappedRGCs(iCell);
        d = optimallyMappedVisualRFmaps{iCell};
        imagesc(d.spatialSupportDegsX,  d.spatialSupportDegsY, d.theRFmap);
        axis 'image'; axis 'xy';
        set(gca, 'CLim', max(abs(d.theRFmap(:)))*[-1 1]);
        colormap(gray(1024));
        title(sprintf('RGC %d: maxRF = %f', theRGCindex, max(d.theRFmap(:))));
        drawnow;
    end

end


function computeRFmapsForAllCellsUsingStimuliAtTargetPosition( ...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, stimulusChromaticity,...
            maxSFLimit, rfMappingPixelMagnificationFactor, ...
            parPoolSize,...
            coneMosaicSubspaceResponsesFileName, ...
            mRGCMosaicSubspaceResponsesFileName, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            visualizedResponses)

    % Encode examined spatial position
    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f_PixelMagnification_%2.3f.mat', stimPositionDegs(1), stimPositionDegs(2), rfMappingPixelMagnificationFactor);
    coneMosaicSubspaceResponsesFileName = strrep(coneMosaicSubspaceResponsesFileName, '.mat', positionPostFix);
    mRGCMosaicSubspaceResponsesFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', positionPostFix);
    optimallyMappedSubspaceRFmapsFileName = strrep(optimallyMappedSubspaceRFmapsFileName, '.mat', positionPostFix);

    % Compute responses of the input cone mosaic to the subspace RF mapping stimuli
    if (reComputeInputConeMosaicSubspaceRFmappingResponses)
        MosaicPoolingOptimizer.generateInputConeMosaicSubspaceRFmappingLinearResponses(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, stimulusChromaticity, ...
            coneMosaicSubspaceResponsesFileName, ...
            'maxSFLimit', maxSFLimit, ...
            'rfMappingPixelMagnificationFactor', rfMappingPixelMagnificationFactor, ...
            'visualizedResponses', visualizedResponses, ...
            'parPoolSize', parPoolSize);
    end

    fprintf('\nLoading cone mosaic subspace modulation responses and Hartley spatial modulation patterns ...');
    % Load the previously computed responses
    load(coneMosaicSubspaceResponsesFileName, ...
        'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'stimParams', 'lIndices', 'mIndices', ...
        'theConeMosaicSubspaceLinearModulationResponses');
    HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);
    fprintf('Done loading !\n');

    [HartleyStimNum, nCones] = size(theConeMosaicSubspaceLinearModulationResponses);

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
             theConeMosaicModulationResponse = squeeze(theConeMosaicSubspaceLinearModulationResponses(iStim,:));
             theConeMosaicModulationResponse = reshape(theConeMosaicModulationResponse, [nTrials nTimePoints nCones]);

             % Compute mRGC mosaic responses based on the cone mosaic modulation responses
             [theMRGCMosaicResponse, ~, theMRGCresponseTemporalSupportSeconds] = ...
                    theComputeReadyMRGCmosaic.compute(theConeMosaicModulationResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicSubspaceRFmappingLinearResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(nTrials, nTimePoints ,:)));
        end
    
        fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses to %s ...', mRGCMosaicSubspaceResponsesFileName);
        save(mRGCMosaicSubspaceResponsesFileName, ...
            'stimParams', ...
            'theMRGCMosaicSubspaceRFmappingLinearResponses', ...
            'spatialSupportDegs', 'lIndices', 'mIndices', ...
            '-v7.3');
    end

    if (reComputeRFs)
        % Load theMRGCMosaicSubspaceRFmappingLinearResponses
        load(mRGCMosaicSubspaceResponsesFileName, ...
                'theMRGCMosaicSubspaceRFmappingLinearResponses');
    
        indicesOfOptimallyMappedRGCs = indicesOfOptimallyMappedRGCsAtThisPosition(theComputeReadyMRGCmosaic, ...
            stimPositionDegs, stimSizeDegs);

        % Compute RF maps of all cells in the MRGC mosaic for stimuli delivered
        % in this position. Only some cells are optimally mapped in each
        % position grid, as the subspace stim size is usually small not covering the
        % entire RGC mosaic. These are selected by the
        % optimalyMappedRFsAtMosaicPosition() function later on
        theMRGCMosaicOptimallyMappedVisualRFmaps = computeRFs(...
            indicesOfOptimallyMappedRGCs, ...
            theMRGCMosaicSubspaceRFmappingLinearResponses, ...
            HartleySpatialModulationPatterns);
    
        fprintf('\nSaving computed visual RFs to %s ...', mRGCMosaicSubspaceResponsesFileName);
        save(mRGCMosaicSubspaceResponsesFileName, ...
            'theMRGCMosaicOptimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs', '-append');
        fprintf('Done saving! \n');

        % Export visual RF maps for cells that are optimally mapped at this  position
        exportOptimalyMappedRFmaps(stimPositionDegs, mRGCMosaicSubspaceResponsesFileName, optimallyMappedSubspaceRFmapsFileName);
    end

end


function theRFmaps = computeRFs( ...
    indicesOfOptimallyMappedRGCs, ...
    theSubspaceRFmappingLinearResponses, ...
    HartleySpatialModulationPatterns)

    nStim = size(theSubspaceRFmappingLinearResponses,1);
    cellsNum = size(theSubspaceRFmappingLinearResponses,2);
    pixelsNum = size(HartleySpatialModulationPatterns,2);

    m = max(abs(theSubspaceRFmappingLinearResponses),[],1);
    cellsWithNonZeroResponse = find(m > 0);
    theRFmaps = cell(cellsNum, 1);

    parfor iCell = 1:numel(indicesOfOptimallyMappedRGCs)
        fprintf('Computing visual RF by accumulating Hartley patterns for the %d of %d optimally mapped RGC... \n', iCell, numel(indicesOfOptimallyMappedRGCs));
        theRGCindex = indicesOfOptimallyMappedRGCs(iCell);
        if (ismember(theRGCindex, cellsWithNonZeroResponse))
            theRFmap = zeros(pixelsNum, pixelsNum, 'single');
            allResponses = squeeze(theSubspaceRFmappingLinearResponses(:,theRGCindex));
            for iStim = 1:nStim
                theRFmap = theRFmap + ...
                        single(squeeze(HartleySpatialModulationPatterns(iStim,:,:)) * allResponses(iStim));
            end
            theRFmaps{iCell} = theRFmap;
        end
    end
end


function indicesOfOptimallyMappedRGCs = indicesOfOptimallyMappedRGCsAtThisPosition(theMRGCmosaic, ...
    stimPositionDegs, stimSizeDegs)

    % Find RGCs within the stimulus region
    if (numel(stimSizeDegs) == 2)
        widthDegs = stimSizeDegs(1);
        heightDegs = stimSizeDegs(2);
    else
        widthDegs = stimSizeDegs(1);
        heightDegs = stimSizeDegs(1);
    end

    theStimulusRegion = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', stimPositionDegs, ...
            'width', widthDegs*0.75, ...
            'height', heightDegs*0.75 , ...
            'rotation', 0.0...
        ));

    indicesOfOptimallyMappedRGCs = theStimulusRegion.indicesOfPointsInside(theMRGCmosaic.rgcRFpositionsDegs);
    
    d = sqrt(sum(theMRGCmosaic.rgcRFpositionsDegs(indicesOfOptimallyMappedRGCs,:).^2,2));
    [~,idx] = sort(d, 'ascend');
    indicesOfOptimallyMappedRGCs = indicesOfOptimallyMappedRGCs(idx);

    % All the cells in the MRGC mosaic
    cellsNum = theMRGCmosaic.rgcsNum;
    
    fprintf('\nThere were %d/%d optimally mapped RF maps at (x,y) = %2.1f,%2.1f within the %2.2f x %2.2f mapping window\n', ...
        numel(indicesOfOptimallyMappedRGCs), cellsNum, ...
        stimPositionDegs(1), stimPositionDegs(2), widthDegs, heightDegs);

end

function exportOptimalyMappedRFmaps(stimPositionDegs, mRGCMosaicSubspaceResponsesFileName, optimallyMappedSubspaceRFmapsFileName)

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
    
    if (~ismember('theMRGCMosaicOptimallyMappedVisualRFmaps', computedVariableNames))
        error('Subspace responses are computed, however RF maps have not been computed yet.\nNot extracting optimally mapped visual RF maps at (x,y) = (%2.1f,%2.1f)\n', ...
            stimPositionDegs(1), stimPositionDegs(2));
    end
    
    fprintf('Loading RF map data from %s. Please wait ...\n', mRGCMosaicSubspaceResponsesFileName);
    % All good. Extract the computed subspace RF maps
    load(mRGCMosaicSubspaceResponsesFileName, ...
        'spatialSupportDegs', ...
        'indicesOfOptimallyMappedRGCs', ...
        'theMRGCMosaicOptimallyMappedVisualRFmaps');


    optimallyMappedRFsForThisGridPosition = numel(indicesOfOptimallyMappedRGCs);
    optimallyMappedVisualRFmaps = cell(1, optimallyMappedRFsForThisGridPosition);

    for iCell = 1:numel(theMRGCMosaicOptimallyMappedVisualRFmaps)
        optimallyMappedVisualRFmaps{iCell} = struct(...
            'theRFmap', theMRGCMosaicOptimallyMappedVisualRFmaps{iCell}, ...
            'spatialSupportDegsX', spatialSupportDegs+stimPositionDegs(1), ...
            'spatialSupportDegsY', spatialSupportDegs+stimPositionDegs(2) ...
            );
    end


    % Save all the optimally mapped visual RF maps
    fprintf('Saving optimally mapped subspace RF maps to %s\n', optimallyMappedSubspaceRFmapsFileName);
    save(optimallyMappedSubspaceRFmapsFileName, ...
        'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs', '-v7.3');

end