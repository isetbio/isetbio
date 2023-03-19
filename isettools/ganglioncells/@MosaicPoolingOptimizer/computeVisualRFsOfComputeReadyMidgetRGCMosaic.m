function computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, theOptics, ...
            maxSFcyclesPerDegree, stimSizeDegs, stimXYpositionDegs, ...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, varargin)

   
    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.parse(varargin{:});
    parPoolSize = p.Results.parPoolSize;

    % Encode examined spatial position
    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f.mat', stimXYpositionDegs(1), stimXYpositionDegs(2));
    coneMosaicResponsesFileName = strrep(coneMosaicResponsesFileName, '.mat', positionPostFix);
    mRGCMosaicResponsesFileName = strrep(mRGCMosaicResponsesFileName, '.mat', positionPostFix);

    % Compute responses of the input cone mosaic to the subspace RF mapping stimuli
    if (reComputeInputConeMosaicSubspaceRFmappingResponses)
        MosaicPoolingOptimizer.generateInputConeMosaicSubspaceRFmappingLinearResponses(...
            theComputeReadyMRGCmosaic, theOptics, ...
            maxSFcyclesPerDegree, stimSizeDegs, stimXYpositionDegs, ...
            coneMosaicResponsesFileName, ...
            'parPoolSize', parPoolSize);
    end
    fprintf('\nLoading cone mosaic subspace responses and Hartley spatial modulation patterns ...');
    % Load the previously computed responses
    load(coneMosaicResponsesFileName, ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceLinearResponses', 'theConeMosaicSubspaceEnergyResponses');
    HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);
    fprintf('Done loading !\n');
    

    % Compute mRGC mosaic responses to theConeMosaicSubspaceResponses
    % Compute RF maps of cells in the mRGC mosaic
    [HartleyStimNum, nCones] = size(theConeMosaicSubspaceLinearResponses);
    
    if (reComputeMRGCMosaicSubspaceRFmappingResponses)
        fprintf('MRGC mosaic subspace RF maps and responses will be saved to %s \n', mRGCMosaicResponsesFileName);
        fprintf('Computing visual subspace RF mapping responses for all RGCs in the mosaic ... \n');
        nTimePoints = 1;
        nTrials = 1;
        theConeMosaicResponseTemporalSupportSeconds = [0];
        
        theMRGCMosaicSubspaceRFmappingLinearResponses = zeros(...
                     HartleyStimNum, theComputeReadyMRGCmosaic.rgcsNum, 'single');
        theMRGCMosaicSubspaceRFmappingEnergyResponses = theMRGCMosaicSubspaceRFmappingLinearResponses;

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


             theConeMosaicResponse = squeeze(theConeMosaicSubspaceEnergyResponses(iStim,:));
             theConeMosaicResponse = reshape(theConeMosaicResponse, [nTrials nTimePoints nCones]);

             % Compute !
             [theMRGCMosaicResponse, theMRGCresponseTemporalSupportSeconds] = ...
                    theComputeReadyMRGCmosaic.compute(theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicSubspaceRFmappingEnergyResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(1, :,:)));

         end
    
    
        fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses to %s ...', mRGCMosaicResponsesFileName);
        save(mRGCMosaicResponsesFileName, ...
            'theMRGCMosaicSubspaceRFmappingLinearResponses', ...
            'theMRGCMosaicSubspaceRFmappingEnergyResponses', ...
            'spatialSupportDegs', 'lIndices', 'mIndices', ...
            '-v7.3');
    end

    % Load theMRGCMosaicSubspaceRFmappingLinearResponses
    load(mRGCMosaicResponsesFileName, ...
            'theMRGCMosaicSubspaceRFmappingLinearResponses', ...
            'theMRGCMosaicSubspaceRFmappingEnergyResponses',...
            'spatialSupportDegs', 'lIndices', 'mIndices');

    % Compute RF maps of cones in the MRGC mosaic
    theMRGCMosaicVisualRFmaps = computeRFs(...
        theComputeReadyMRGCmosaic, ...
        theMRGCMosaicSubspaceRFmappingLinearResponses, ...
        theMRGCMosaicSubspaceRFmappingEnergyResponses, ...
        HartleySpatialModulationPatterns, ...
        spatialSupportDegs, lIndices, mIndices);

    fprintf('\nSaving computed visual RFs to %s ...', mRGCMosaicResponsesFileName);
    save(mRGCMosaicResponsesFileName, ...
        'theMRGCMosaicVisualRFmaps', '-append');
    fprintf('Done saving! \n')
end

function theRFmaps = computeRFs(theComputeReadyMRGCmosaic, ...
    theSubspaceRFmappingLinearResponses, ...
    theSubspaceRFmappingEnergyResponses, ...
    HartleySpatialModulationPatterns, ...
    spatialSupportDegs, lIndices, mIndices)

    nStim = size(theSubspaceRFmappingLinearResponses,1);
    omega = (sqrt(nStim)-1)/2;
    cellsNum = size(theSubspaceRFmappingLinearResponses,2);
    pixelsNum = size(HartleySpatialModulationPatterns,2);

    fprintf('Computing visual RFs for all RGCs in the mosaic ... \n')
    theRFmaps = cell(cellsNum, 1);

    hFig = figure(1);
    set(hFig, 'Position', [10 10 2048 1024], 'Color', [1 1 1]);
    ax1 = subplot(1,2,1);
    ax2 = subplot(1,2,2);

    m = max(abs(theSubspaceRFmappingLinearResponses),[],1);
    cellsWithNonZeroResponse = find(m > 0);

    for idx = 1:numel(cellsWithNonZeroResponse)

        iCell = cellsWithNonZeroResponse(idx);
        HartleyMapSize = 2*omega+1;
        theHartleyTuningMap = zeros(HartleyMapSize, HartleyMapSize);
        for iStim = 1:nStim
            theHartleyTuningMap(lIndices(iStim)+omega+1, mIndices(iStim)+omega+1) = theSubspaceRFmappingEnergyResponses(iStim,iCell);
        end

        theHartleyTuningMap = theHartleyTuningMap / max(theHartleyTuningMap(:));
        
        imagesc(ax1,-omega:1:omega, -omega:1:omega, abs(theHartleyTuningMap));
        set(ax1, 'CLim', [0 1]);
        axis(ax1, 'image')
        colormap(ax1,brewermap(1024, '*greys'));

        theRFmap = zeros(pixelsNum, pixelsNum, 'single');
        for iStim = 1:nStim
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
  
        theRFmaps{iCell} = theRFmap;
    end

end


