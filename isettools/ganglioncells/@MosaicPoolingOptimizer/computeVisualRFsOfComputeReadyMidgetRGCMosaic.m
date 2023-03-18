function computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, ...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses)

   
    % Compute responses of the input cone mosaic to the subspace RF mapping
    % stimuli

    if (reComputeInputConeMosaicSubspaceRFmappingResponses)
        MosaicPoolingOptimizer.generateInputConeMosaicSubspaceRFmappingLinearResponses(...
         theComputeReadyMRGCmosaic, coneMosaicResponsesFileName);
    else
        fprintf('\nLoading cone mosaic subspace responses and Hartley spatial modulation patterns ...');
        % Load the previously computed responses
        load(coneMosaicResponsesFileName, ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceResponses');
        HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);
        fprintf('Done loading !\n');
    end

    % Compute mRGC mosaic responses to theConeMosaicSubspaceResponses
    % Compute RF maps of cells in the mRGC mosaic

    [HartleyStimNum, nCones] = size(theConeMosaicSubspaceResponses);
    
    fprintf('MRGC mosaic subspace RF maps and responses will be saved to %s \n', mRGCMosaicResponsesFileName);
    fprintf('Computing visual subspace RF mapping responses for all RGCs in the mosaic ... \n');
    nTimePoints = 1;
    nTrials = 1;
    theConeMosaicResponseTemporalSupportSeconds = [0];
    
    theMRGCMosaicSubspaceRFmappingLinearResponses = zeros(...
                 HartleyStimNum, theComputeReadyMRGCmosaic.rgcsNum, 'single');
    parfor iStim = 1:HartleyStimNum
         theConeMosaicResponse = squeeze(theConeMosaicSubspaceResponses(iStim,:));
         theConeMosaicResponse = reshape(theConeMosaicResponse, [nTrials nTimePoints nCones]);
         % Compute !
         [theMRGCMosaicResponse, theMRGCresponseTemporalSupportSeconds] = ...
                theComputeReadyMRGCmosaic.compute(theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds);
         theMRGCMosaicSubspaceRFmappingLinearResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(1, :,:)));
     end


    % Compute RF maps of cones in the MRGC mosaic
    theMRGCMosaicVisualRFmaps = computeRFs(theComputeReadyMRGCmosaic, theMRGCMosaicSubspaceRFmappingLinearResponses, HartleySpatialModulationPatterns);

    fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses and computed visual RFs to %s ...', mRGCMosaicResponsesFileName);
    save(mRGCMosaicResponsesFileName, ...
        'theMRGCMosaicSubspaceRFmappingLinearResponses', ...
        'spatialSupportDegs', 'lIndices', 'mIndices', ...
        'theMRGCMosaicVisualRFmaps', '-v7.3');
    fprintf('Done saving! \n')
end

function theRFmaps = computeRFs(theComputeReadyMRGCmosaic, theSubspaceRFmappingLinearResponses, HartleySpatialModulationPatterns)
    nStim = size(theSubspaceRFmappingLinearResponses,1);
    cellsNum = size(theSubspaceRFmappingLinearResponses,2);
    pixelsNum = size(HartleySpatialModulationPatterns,2);

    fprintf('Computing visual RFs for all RGCs in the mosaic ... \n')
    theRFmaps = cell(cellsNum, 1);

     hFig = figure(1);
     set(hFig, 'Position', [10 10 1024 1024], 'Color', [1 1 1]);
     ax = subplot('Position', [0.05 0.05 0.94 0.94]);
    for iCell = 1:cellsNum
        theRFmap = zeros(pixelsNum, pixelsNum, 'single');
        for iStim = 1:nStim
        fprintf('Updating RFs based on responses to Hartley pattern %d of %d\n', iStim, nStim);
            r = theSubspaceRFmappingLinearResponses(iStim,iCell);
            theRFmap = theRFmap + ...
                    single(squeeze(HartleySpatialModulationPatterns(iStim,:,:)) * r);
        end
        theRFmap = theRFmap / max(abs(theRFmap(:)));
       
        imagesc(ax,theRFmap, [-1 1]);
        axis(ax, 'image')
        title(ax,sprintf('RF %d of %d (%2.1f,%2.1f) degs', iCell, cellsNum, ...
            theComputeReadyMRGCmosaic.rgcRFpositionsDegs(iCell,1), ...
            theComputeReadyMRGCmosaic.rgcRFpositionsDegs(iCell,2)));
        colormap(ax,brewermap(1024, '*RdBu'));
        drawnow
        pause
        theRFmaps{iCell} = theRFmap;
    end

end


