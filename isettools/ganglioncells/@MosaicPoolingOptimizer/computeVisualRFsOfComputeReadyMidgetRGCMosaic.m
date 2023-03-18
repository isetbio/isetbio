function computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, ...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses)

   
    % Compute responses of the input cone mosaic to the subspace RF mapping
    % stimuli

    if (reComputeInputConeMosaicSubspaceRFmappingResponses)
        MosaicPoolingOptimizer.generateInputConeMosaicSubspaceRFmappingLinearResponses(...
         theComputeReadyMRGCmosaic, coneMosaicResponsesFileName);
    end
    fprintf('\nLoading cone mosaic subspace responses and Hartley spatial modulation patterns ...');
    % Load the previously computed responses
    load(coneMosaicResponsesFileName, ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceResponses');
    HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);
    fprintf('Done loading !\n');
    

    % Compute mRGC mosaic responses to theConeMosaicSubspaceResponses
    % Compute RF maps of cells in the mRGC mosaic

    [HartleyStimNum, nCones] = size(theConeMosaicSubspaceResponses);
    
    if (reComputeMRGCMosaicSubspaceRFmappingResponses)
        fprintf('MRGC mosaic subspace RF maps and responses will be saved to %s \n', mRGCMosaicResponsesFileName);
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
             theConeMosaicResponse = squeeze(theConeMosaicSubspaceResponses(iStim,:));
             theConeMosaicResponse = reshape(theConeMosaicResponse, [nTrials nTimePoints nCones]);
             % Compute !
             [theMRGCMosaicResponse, theMRGCresponseTemporalSupportSeconds] = ...
                    theComputeReadyMRGCmosaic.compute(theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicSubspaceRFmappingLinearResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(1, :,:)));
         end
    
    
        fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses to %s ...', mRGCMosaicResponsesFileName);
        save(mRGCMosaicResponsesFileName, ...
            'theMRGCMosaicSubspaceRFmappingLinearResponses', ...
            'spatialSupportDegs', 'lIndices', 'mIndices', ...
            '-v7.3');
    end

    % Load theMRGCMosaicSubspaceRFmappingLinearResponses
    load(mRGCMosaicResponsesFileName, ...
            'theMRGCMosaicSubspaceRFmappingLinearResponses', ...
            'spatialSupportDegs', 'lIndices', 'mIndices');

    % Compute RF maps of cones in the MRGC mosaic
    theMRGCMosaicVisualRFmaps = computeRFs(...
        theComputeReadyMRGCmosaic, ...
        theMRGCMosaicSubspaceRFmappingLinearResponses, ...
        HartleySpatialModulationPatterns, ...
        spatialSupportDegs, lIndices, mIndices);

    fprintf('\nSaving computed visual RFs to %s ...', mRGCMosaicResponsesFileName);
    save(mRGCMosaicResponsesFileName, ...
        'theMRGCMosaicVisualRFmaps', '-append', '-v7.3');
    fprintf('Done saving! \n')
end

function theRFmaps = computeRFs(theComputeReadyMRGCmosaic, ...
    theSubspaceRFmappingLinearResponses, ...
    HartleySpatialModulationPatterns, ...
    spatialSupportDegs, lIndices, mIndices)

    nStim = size(theSubspaceRFmappingLinearResponses,1);
    cellsNum = size(theSubspaceRFmappingLinearResponses,2);
    pixelsNum = size(HartleySpatialModulationPatterns,2);

    fprintf('Computing visual RFs for all RGCs in the mosaic ... \n')
    theRFmaps = cell(cellsNum, 1);

    hFig = figure(1);
    set(hFig, 'Position', [10 10 2048 1024], 'Color', [1 1 1]);
    ax1 = subplot(1,2,1);
    ax2 = subplot(1,2,2);

    for iCell = 1:cellsNum

        HartleyMapSize = sqrt(nStim);
        theHartleyTuningMap = zeros(HartleyMapSize, HartleyMapSize);
        for iStim = 1:nStim
            theHartleyTuningMap(lIndices(iStim), mIndices(iStim)) = theSubspaceRFmappingLinearResponses(iStim,iCell);
        end
        theHartleyTuningMap = theHartleyTuningMap / max(abs(theHartleyTuningMap(:)));
        imagesc(ax1,theHartleyTuningMap);
        set(ax, 'CLim', [-1 1]);
        axis(ax1, 'image')
        colormap(ax1,brewermap(1024, '*RdBu'));

        theRFmap = zeros(pixelsNum, pixelsNum, 'single');
        for iStim = 1:nStim
            fprintf('Updating RF of cell %d based on responses to Hartley pattern %d of %d\n', iCell, iStim, nStim);
            r = theSubspaceRFmappingLinearResponses(iStim,iCell);
            theRFmap = theRFmap + ...
                    single(squeeze(HartleySpatialModulationPatterns(iStim,:,:)) * r);
        end
        theRFmap = theRFmap / max(abs(theRFmap(:)));
       
        imagesc(ax2,theRFmap);
        set(ax2, 'CLim', [-1 1]);
        axis(ax2, 'image')
        title(ax2,sprintf('RF %d of %d (%2.1f,%2.1f) degs', iCell, cellsNum, ...
            theComputeReadyMRGCmosaic.rgcRFpositionsDegs(iCell,1), ...
            theComputeReadyMRGCmosaic.rgcRFpositionsDegs(iCell,2)));
        colormap(ax2,brewermap(1024, '*RdBu'));
        drawnow
        pause
        theRFmaps{iCell} = theRFmap;
    end

end


