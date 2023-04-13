function performComputeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic(mosaicParams, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('maxSFcyclesPerDegree', 30, @isscalar);
    p.addParameter('stimSizeDegs', 1.0, @isscalar);
    p.addParameter('posIncrementDegs', 0.5, @isscalar);
    p.addParameter('gridNodeIndicesToCompute', [], @(x)(isempty(x)||(isnumeric(x))));
    p.addParameter('reComputeInputConeMosaicSubspaceRFmappingResponses', false, @islogical);
    p.addParameter('reComputeMRGCMosaicSubspaceRFmappingResponses', false, @islogical);
    p.addParameter('reComputeRFs', false, @islogical);
    p.addParameter('onlyVisualizeOptimallyMappedRFmaps', false, @islogical);
    p.parse(varargin{:});

    % Subspace RF mapping params:
    % max SF to explore
    maxSFcyclesPerDegree = p.Results.maxSFcyclesPerDegree; 

    % stimulus patch size
    stimSizeDegs = p.Results.stimSizeDegs;    
    
    % simulus delta position
    posIncrementDegs = p.Results.posIncrementDegs;

    % Which grid nodes to compute RFs
    gridNodeIndicesToCompute = p.Results.gridNodeIndicesToCompute;

    reComputeInputConeMosaicSubspaceRFmappingResponses = p.Results.reComputeInputConeMosaicSubspaceRFmappingResponses;
    reComputeMRGCMosaicSubspaceRFmappingResponses = p.Results.reComputeMRGCMosaicSubspaceRFmappingResponses;
    reComputeRFs = p.Results.reComputeRFs;
    onlyVisualizeOptimallyMappedRFmaps = p.Results.onlyVisualizeOptimallyMappedRFmaps;

    parpoolSize = [];

    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Ask the user which H1 cell index to use for optimizing the RF
    % surround pooling model
    retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);

    % Generate the filename of the compute-ready mRGCMosaic to generate
    [computeReadyMosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    % Load the compute-ready MRGC mosaic
    load(fullfile(resourcesDirectory, computeReadyMosaicFileName), 'theComputeReadyMRGCmosaic');

    % Generate filename for the coneMosaic subspace responses
    [coneMosaicSubspaceResponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('coneMosaicSubspaceResponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);

    % Generate filename for the mRGCMosaic subspace responses
    [mRGCMosaicSubspaceResponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSubspaceResponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);

    % Optimally generated RF maps filename
    optimallyMappedSubspaceRFmapsFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', '_optimallyMappedRFs.mat');

    MosaicPoolingOptimizer.computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, opticsParams, ...
            maxSFcyclesPerDegree, stimSizeDegs, posIncrementDegs, ...
            fullfile(resourcesDirectory, coneMosaicSubspaceResponsesFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSubspaceResponsesFileName), ...
            fullfile(resourcesDirectory, optimallyMappedSubspaceRFmapsFileName), ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            onlyVisualizeOptimallyMappedRFmaps, ...
            'parPoolSize', parpoolSize);
end
