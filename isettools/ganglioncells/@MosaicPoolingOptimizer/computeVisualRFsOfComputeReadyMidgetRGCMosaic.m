function computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            computeReadyMosaicFilename, ...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName)

    % Load the compute-ready MRGC mosaic
    load(computeReadyMosaicFilename, 'theComputeReadyMRGCmosaic');

    % Compute responses of the input cone mosaic to the subspace RF mapping
    % stimuli
    MosaicPoolingOptimizer.generateInputConeMosaicSubspaceRFmappingResponses(...
        theComputeReadyMRGCmosaic, coneMosaicResponsesFileName);

    % Load computed responses
    load(coneMosaicResponsesFileName, ...
                'HartleySpatialModulationPatterns', ...
                'theConeMosaicSubspaceResponses');

    % Compute RF maps of cones in the cone mosaic
    theRFmaps = computeRFs(theConeMosaicSubspaceResponses, HartleySpatialModulationPatterns);

    % Compute mRGC mosaic responses to theConeMosaicSubspaceResponses
    % Compute RF maps of cells in the mRGC mosaic
end

function theRFmaps = computeRFs(theSubspaceRFmappingResponses, HartleySpatialModulationPatterns)
    nStim = size(theSubspaceRFmappingResponses,1);
    cellsNum = size(theSubspaceRFmappingResponses,2);
    pixelsNum = size(HartleySpatialModulationPatterns,2);

    theRFmaps = zeros(cellsNum, pixelsNum, pixelsNum);
    for iFrame = 1:nStim
        parfor iCell = 1:cellsNum
            r = theSubspaceRFmappingResponses(iFrame,iCell);
            theRFmaps(iCell,:,:) = theRFmaps(iCell,:,:) + ...
                    HartleySpatialModulationPatterns(iFrame,:,:) * r;
        end
    end

end


