function computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, opticsParams, ...
            maxSFcyclesPerDegree, stimSizeDegs, posIncrementDegs, ...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, varargin)

   
    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('visualizedResponses', false, @islogical);
    p.parse(varargin{:});
    parPoolSize = p.Results.parPoolSize;
    visualizedResponses = p.Results.visualizedResponses;

    [X,Y] = generateSamplingGrid(...
        theComputeReadyMRGCmosaic.inputConeMosaic.sizeDegs, ...
        theComputeReadyMRGCmosaic.inputConeMosaic.eccentricityDegs, ...
        posIncrementDegs);

    stimXYpositionGridDegs = [X(:) Y(:)];

    for theGridNodeIndex = 1:size(stimXYpositionGridDegs,1)
        % stimulus position within the mRGC mosaic
        gridNodeXYpositionDegs = stimXYpositionGridDegs(theGridNodeIndex,:);
    
        % Generate native optics
        opticsParamsAtThisPosition = opticsParams;
        opticsParamsAtThisPosition.positionDegs = gridNodeXYpositionDegs;
        theComputeReadyMRGCmosaic.generateNativeOptics(opticsParams);

        % Retrieve the native optics
        theOptics = theComputeReadyMRGCmosaic.theNativeOptics;

        if (reComputeInputConeMosaicSubspaceRFmappingResponses || ...
            reComputeMRGCMosaicSubspaceRFmappingResponses || ...
            reComputeRFs)

            computeRFsForGridPosition(...
                theComputeReadyMRGCmosaic, theOptics, ...
                maxSFcyclesPerDegree, stimSizeDegs, gridNodeXYpositionDegs, ...
                parPoolSize, ...
                coneMosaicResponsesFileName, ...
                mRGCMosaicResponsesFileName, ...
                reComputeInputConeMosaicSubspaceRFmappingResponses, ...
                reComputeMRGCMosaicSubspaceRFmappingResponses, ...
                visualizedResponses);

        else
            visualizeRFsForGridPosition(theComputeReadyMRGCmosaic, stimXYpositionGridDegs, theGridNodeIndex, mRGCMosaicResponsesFileName);
        end
    end

end

function visualizeRFsForGridPosition(theComputeReadyMRGCmosaic, stimXYpositionGridDegs, theGridNodeIndex, mRGCMosaicResponsesFileName)

    % Encode examined spatial position
    gridNodeXYpositionDegs = stimXYpositionGridDegs(theGridNodeIndex,:);
    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f.mat', gridNodeXYpositionDegs(1), gridNodeXYpositionDegs(2));
    mRGCMosaicResponsesFileName = strrep(mRGCMosaicResponsesFileName, '.mat', positionPostFix);

    load(mRGCMosaicResponsesFileName, ...
        'spatialSupportDegs', 'lIndices', 'mIndices', ...
        'theMRGCMosaicVisualRFmaps');

    HartleyMapSize = sqrt(numel(lIndices));

    nStim = size(theMRGCMosaicSubspaceRFmappingEnergyResponses,1);
    omega = (HartleyMapSize-1)/2;

    cellsNum = numel(theMRGCMosaicVisualRFmaps);
    for iCell = 1:cellsNum

        if isempty(theMRGCMosaicVisualRFmaps{iCell})
            continue;
        end

        % Find the grid position that this cell is closest to
        d = sqrt(sum((bsxfun(@minus, stimXYpositionGridDegs, theComputeReadyMRGCmosaic.rgcRFpositionsDegs(iCell,:))).^2,2));
        [~, theClosestGridNodeIndex] = min(d);

        if (theClosestGridNodeIndex ~= theGridNodeIndex)
            continue;
        end

        hFig = figure(1);clf;
        set(hFig, 'Position', [10 10 1800 500], 'Color', [1 1 1]);
        ax1 = subplot(1,3,1);
        ax2 = subplot(1,3,2);
        ax3 = subplot(1,3,3);

        theRFmap = theMRGCMosaicVisualRFmaps{iCell};
        imagesc(ax2, spatialSupportDegs, spatialSupportDegs, theRFmap);
        
        set(ax2, 'CLim', 0.1*[-1 1]);
        axis(ax2, 'image')
        title(ax2,sprintf('RF at (%2.1f,%2.1f), grid node position: (%2.1f, %2.1f)', ...
                theComputeReadyMRGCmosaic.rgcRFpositionsDegs(iCell,1), ...
                theComputeReadyMRGCmosaic.rgcRFpositionsDegs(iCell,2), ...
                gridNodeXYpositionDegs(1), gridNodeXYpositionDegs(2)));
        colormap(ax2,brewermap(1024, '*RdBu'));
        
        theRFprofileX = squeeze(sum(theRFmap,1));
        theRFprofileY = squeeze(sum(theRFmap,2));
        maxProfile = max([max(abs(theRFprofileX(:))) max(abs(theRFprofileY(:)))]);

        plot(ax3, spatialSupportDegs, theRFprofileX/maxProfile, 'r-', 'LineWidth', 1.5);
        hold(ax3, 'on');
        plot(ax3, spatialSupportDegs, theRFprofileY/maxProfile, 'b-', 'LineWidth', 1.5);
        set(ax3, 'YLim', [-1 1]);
        pause
        drawnow

    end
end


function [X,Y] = generateSamplingGrid(inputConeMosaicSizeDegs, inputConeMosaicEccDegs, posIncrementDegs)

    k = round(0.5*inputConeMosaicSizeDegs(1)/posIncrementDegs)-1;
    xCoords = inputConeMosaicEccDegs(1) + ...
            (-(k+1):1:k)*posIncrementDegs + posIncrementDegs*0.5;

    k = round(0.5*inputConeMosaicSizeDegs(2)/posIncrementDegs)-1;
    yCoords = inputConeMosaicEccDegs(2) + ...
            (-(k+1):1:k)*posIncrementDegs + posIncrementDegs*0.5;

    [X,Y] = meshgrid(xCoords, yCoords);
    X = X(:); Y = Y(:);
end


function computeRFsForGridPosition( ...
            theComputeReadyMRGCmosaic, theOptics, ...
            maxSFcyclesPerDegree, stimSizeDegs, gridNodeXYpositionDegs, ...
            parPoolSize,...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            visualizedResponses)

    % Encode examined spatial position
    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f.mat', gridNodeXYpositionDegs(1), gridNodeXYpositionDegs(2));
    coneMosaicResponsesFileName = strrep(coneMosaicResponsesFileName, '.mat', positionPostFix);
    mRGCMosaicResponsesFileName = strrep(mRGCMosaicResponsesFileName, '.mat', positionPostFix);

    % Compute responses of the input cone mosaic to the subspace RF mapping stimuli
    if (reComputeInputConeMosaicSubspaceRFmappingResponses)
        MosaicPoolingOptimizer.generateInputConeMosaicSubspaceRFmappingLinearResponses(...
            theComputeReadyMRGCmosaic, theOptics, ...
            maxSFcyclesPerDegree, stimSizeDegs, gridNodeXYpositionDegs, ...
            coneMosaicResponsesFileName, ...
            'visualizedResponses', visualizedResponses, ...
            'parPoolSize', parPoolSize);
    end
    fprintf('\nLoading cone mosaic subspace responses and Hartley spatial modulation patterns ...');
    % Load the previously computed responses
    load(coneMosaicResponsesFileName, ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceLinearResponses');
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
    
    
        fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses to %s ...', mRGCMosaicResponsesFileName);
        save(mRGCMosaicResponsesFileName, ...
            'theMRGCMosaicSubspaceRFmappingLinearResponses', ...
            'spatialSupportDegs', 'lIndices', 'mIndices', ...
            '-v7.3');
    end

    % Load theMRGCMosaicSubspaceRFmappingLinearResponses
    load(mRGCMosaicResponsesFileName, ...
            'theMRGCMosaicSubspaceRFmappingLinearResponses');

    % Compute RF maps of cones in the MRGC mosaic
    theMRGCMosaicVisualRFmaps = computeRFs(...
        theMRGCMosaicSubspaceRFmappingLinearResponses, ...
        HartleySpatialModulationPatterns);

    fprintf('\nSaving computed visual RFs to %s ...', mRGCMosaicResponsesFileName);
    save(mRGCMosaicResponsesFileName, ...
        'theMRGCMosaicVisualRFmaps', '-append');
    fprintf('Done saving! \n')
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


