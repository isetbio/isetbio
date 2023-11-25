function computeVisualRFsUsingMSequenceMapping(...
        theComputeReadyMRGCmosaic, opticsToEmploy, ...
        stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
        stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
        ternaryInsteadOfBinaryMsequence, ...
        coneMosaicResponsesFileName, ...
        mRGCMosaicResponsesFileName, ...
        optimallyMappedRFmapsFileName, ...
        reComputeInputConeMosaicResponses, ...
        reComputeMRGCMosaicResponses, ...
        reComputeRFs, ...
        visualizeOptimallyMappedRFmapLocations, varargin)

    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.parse(varargin{:});
    parPoolSize = p.Results.parPoolSize;

    visualizedResponses = false;

    if (isempty(stimPositionDegs))
        stimPositionDegs = theComputeReadyMRGCmosaic.eccentricityDegs;
    end

    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f_PixelsAcross_%d.mat', stimPositionDegs(1), stimPositionDegs(2), rfPixelsAcross);

    coneMosaicResponsesFileName = strrep(coneMosaicResponsesFileName, '.mat', positionPostFix);
    mRGCMosaicResponsesFileName = strrep(mRGCMosaicResponsesFileName, '.mat', positionPostFix);
    optimallyMappedRFmapsFileName = strrep(optimallyMappedRFmapsFileName, '.mat', positionPostFix);

    if (reComputeInputConeMosaicResponses || reComputeMRGCMosaicResponses || reComputeRFs)

        % Compute the RF maps for ALL cells using stimuli at this 
        % position. Only some cells will be optimally mapped at each
        % position. These optimally derived RF maps are extracted by
        % the optimalyMappedRFsAtMosaicPosition() function
        computeRFmapsForAllCellsUsingStimuliAtTargetPosition(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            ternaryInsteadOfBinaryMsequence, ...
            parPoolSize, ...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            optimallyMappedRFmapsFileName, ...
            reComputeInputConeMosaicResponses, ...
            reComputeMRGCMosaicResponses, ...
            reComputeRFs, ...
            visualizedResponses);
    end

    if (visualizeOptimallyMappedRFmapLocations)
        visualizeAllOptimallyMappedRFmapLocations(optimallyMappedRFmapsFileName, ...
            theComputeReadyMRGCmosaic, visualizedRGCindex, stimulusChromaticity);
    end

end

function computeRFmapsForAllCellsUsingStimuliAtTargetPosition( ...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            ternaryInsteadOfBinaryMsequence, ...
            parPoolSize,...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            optimallyMappedRFmapsFileName, ...
            reComputeInputConeMosaicResponses, ...
            reComputeMRGCMosaicResponses, ...
            reComputeRFs, ...
            visualizedResponses)

    % Compute responses of the input cone mosaic to the subspace RF mapping stimuli
    if (reComputeInputConeMosaicResponses)
        fprintf('Cone mosaic m-sequence responses will be saved to %s \n', coneMosaicResponsesFileName);
        MosaicPoolingOptimizer.generateInputConeMosaicMSequenceRFmappingLinearResponses(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            ternaryInsteadOfBinaryMsequence, ...
            coneMosaicResponsesFileName, ...
            'visualizedResponses', visualizedResponses, ...
            'parPoolSize', parPoolSize);
    end

    if (reComputeMRGCMosaicResponses)
        fprintf('\nLoading cone mosaic m-sequence modulation responses and spatial modulation patterns ...');
        % Load the previously computed responses
        load(coneMosaicResponsesFileName, ...
            'mSequenceSpatialModulationPatterns', 'spatialSupportDegs', 'stimParams', ...
            'theConeMosaicMSequenceLinearModulationResponses');
        mSequenceSpatialModulationPatterns = single(mSequenceSpatialModulationPatterns);
        fprintf('Done loading !\n');

        [mSequenceStimsNum, nCones] = size(theConeMosaicMSequenceLinearModulationResponses);

        fprintf('MRGC mosaic m-sequence RF maps and responses will be saved to %s \n', mRGCMosaicResponsesFileName);
        fprintf('Computing visual subspace RF mapping responses for all RGCs in the mosaic ... \n');
        nTimePoints = 1;
        nTrials = 1;
        theConeMosaicResponseTemporalSupportSeconds = [0];
        
        theMRGCMosaicMSequenceRFmappingLinearResponses = zeros(...
                     mSequenceStimsNum, theComputeReadyMRGCmosaic.rgcsNum, 'single');

        % Use all processors
        [shutdownParPoolOnceCompleted, numWorkers] = MosaicPoolingOptimizer.resetParPool([]);
    
        % Noise-free responses
        theComputeReadyMRGCmosaic.noiseFlag = 'none';

        parfor iStim = 1:mSequenceStimsNum
             fprintf('Computing mRGC mosaic response for m-sequence frame %d of %d (using %d parallel processes).\n', ...
                 iStim, mSequenceStimsNum, numWorkers);
             theConeMosaicModulationResponse = squeeze(theConeMosaicMSequenceLinearModulationResponses(iStim,:));
             theConeMosaicModulationResponse = reshape(theConeMosaicModulationResponse, [nTrials nTimePoints nCones]);

             % Compute mRGC mosaic responses based on the cone mosaic modulation responses
             [theMRGCMosaicResponse, ~, theMRGCresponseTemporalSupportSeconds] = ...
                    theComputeReadyMRGCmosaic.compute(theConeMosaicModulationResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicMSequenceRFmappingLinearResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(nTrials, nTimePoints ,:)));
        end
    
        fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses to %s ...', mRGCMosaicResponsesFileName);
        save(mRGCMosaicResponsesFileName, ...
            'stimParams', ...
            'theMRGCMosaicMSequenceRFmappingLinearResponses', ...
            'spatialSupportDegs', ...
            '-v7.3');
    end


    if (reComputeRFs)
        load(coneMosaicResponsesFileName, 'mSequenceSpatialModulationPatterns', 'spatialSupportDegs', 'stimParams');
        mSequenceSpatialModulationPatterns = single(mSequenceSpatialModulationPatterns);

        % Load theMRGCMosaicMSequenceRFmappingLinearResponses
        load(mRGCMosaicResponsesFileName, 'theMRGCMosaicMSequenceRFmappingLinearResponses');
    
  
        % Determine indices of RGCs whose RF lie within the stimulus region
        indicesOfOptimallyMappedRGCs = MosaicPoolingOptimizer.indicesOfOptimallyMappedRGCsAtThisPosition(theComputeReadyMRGCmosaic, ...
            stimPositionDegs, stimSizeDegs);
        
        % Compute RF maps of all cells within the stimulus region
        theMRGCMosaicOptimallyMappedVisualRFmaps = computeRFs(...
            indicesOfOptimallyMappedRGCs, ...
            theMRGCMosaicMSequenceRFmappingLinearResponses, ...
            mSequenceSpatialModulationPatterns, ...
            spatialSupportDegs);
    
        fprintf('\nSaving computed visual RFs to %s ...', mRGCMosaicResponsesFileName);
        save(mRGCMosaicResponsesFileName, ...
            'theMRGCMosaicOptimallyMappedVisualRFmaps', ...
            'indicesOfOptimallyMappedRGCs', '-append');
        fprintf('Done saving! \n');

        % Export visual RF maps for cells that are optimally mapped at this  position
        exportOptimalyMappedRFmaps(stimPositionDegs, mRGCMosaicResponsesFileName, optimallyMappedRFmapsFileName);
    end

end

function exportOptimalyMappedRFmaps(stimPositionDegs, mRGCMosaicResponsesFileName, optimallyMappedRFmapsFileName)
    
    fprintf('Loading RF map data from %s. Please wait ...\n', mRGCMosaicResponsesFileName);
    % All good. Extract the computed subspace RF maps
    load(mRGCMosaicResponsesFileName, ...
        'spatialSupportDegs', ...
        'indicesOfOptimallyMappedRGCs', ...
        'theMRGCMosaicOptimallyMappedVisualRFmaps');

    optimallyMappedVisualRFmaps = cell(1, numel(indicesOfOptimallyMappedRGCs));

    dx = spatialSupportDegs(2)-spatialSupportDegs(1);
    for iCell = 1:numel(theMRGCMosaicOptimallyMappedVisualRFmaps)
        optimallyMappedVisualRFmaps{iCell} = struct(...
            'theRFmap', theMRGCMosaicOptimallyMappedVisualRFmaps{iCell}, ...
            'spatialSupportDegsX', spatialSupportDegs+stimPositionDegs(1)-dx, ...
            'spatialSupportDegsY', spatialSupportDegs+stimPositionDegs(2)-dx);
    end

    % Save all the optimally mapped visual RF maps
    fprintf('Saving optimally mapped subspace RF maps to %s\n', optimallyMappedRFmapsFileName);
    save(optimallyMappedRFmapsFileName, ...
        'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs', '-v7.3');

end

function theRFmaps = computeRFs(indicesOfOptimallyMappedRGCs, ...
    theRFmappingLinearResponses, spatialModulationPatterns, spatialSupportDegs)

    nStim = size(theRFmappingLinearResponses,1);
    cellsNum = size(theRFmappingLinearResponses,2);
    pixelsNum = size(spatialModulationPatterns,2);

    m = max(abs(theRFmappingLinearResponses),[],1);
    cellsWithNonZeroResponse = find(m > 0);
    theRFmaps = cell(cellsNum, 1);

    parfor iCell = 1:numel(indicesOfOptimallyMappedRGCs)
        fprintf('Computing visual RF by accumulating m-sequence frames for the %d of %d optimally mapped RGC... \n', iCell, numel(indicesOfOptimallyMappedRGCs));
        theRGCindex = indicesOfOptimallyMappedRGCs(iCell);

        if (ismember(theRGCindex, cellsWithNonZeroResponse))
            theRFmap = zeros(pixelsNum, pixelsNum, 'single');
            allResponses = squeeze(theRFmappingLinearResponses(:,theRGCindex));

            for iStim = 1:nStim
                theRFmap = theRFmap + single(squeeze(spatialModulationPatterns(iStim,:,:)) * allResponses(iStim));
            end
            theRFmaps{iCell} = theRFmap;
        end

    end
end


function visualizeAllOptimallyMappedRFmapLocations(optimallyMappedRFmapsFileName,...
    theMRGCMosaic, visualizedRGCindex, stimulusChromaticity)

    % Load all the optimally mapped visual RF maps
    fprintf('Loading optimally mapped subspace RF maps from %s\n', optimallyMappedRFmapsFileName);
    load(optimallyMappedRFmapsFileName, 'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs');


end

