function performComputeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic(mosaicParams, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('maxSFcyclesPerDegree', 30, @isscalar);
    p.addParameter('stimSizeDegs', 1.0, @isscalar);
    p.addParameter('stimPositionDegs', [], @(x)(isempty(x)||(numel(x) == 2)));
    p.addParameter('reComputeInputConeMosaicSubspaceRFmappingResponses', false, @islogical);
    p.addParameter('reComputeMRGCMosaicSubspaceRFmappingResponses', false, @islogical);
    p.addParameter('reComputeRFs', false, @islogical);
    p.addParameter('visualizeOptimallyMappedRFmapLocations', false, @islogical);
    p.parse(varargin{:});

    % Subspace RF mapping params:
    % max SF to explore
    maxSFcyclesPerDegree = p.Results.maxSFcyclesPerDegree; 

    % stimulus patch size
    stimSizeDegs = p.Results.stimSizeDegs;    
    
    % simulus position
    stimPositionDegs = p.Results.stimPositionDegs;


    % Which grid nodes to compute RFs
    reComputeInputConeMosaicSubspaceRFmappingResponses = p.Results.reComputeInputConeMosaicSubspaceRFmappingResponses;
    reComputeMRGCMosaicSubspaceRFmappingResponses = p.Results.reComputeMRGCMosaicSubspaceRFmappingResponses;
    reComputeRFs = p.Results.reComputeRFs;
    visualizeOptimallyMappedRFmapLocations = p.Results.visualizeOptimallyMappedRFmapLocations;

    parpoolSize = [];

    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    fprintf('\n---> Select the optics that were used to compute the compute-ready mosaic\n');
    opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);


    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    fprintf('\n---> Select the optics to use for subspace RF mapping\n');
    [opticsParamsForSubSpaceMapping, opticsToEmploy] = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);


    % Ask the user which H1 cell index to use for optimizing the RF
    % surround pooling model
    retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);

    % Generate the filename of the compute-ready mRGCMosaic to generate
    [computeReadyMosaicFileName, computeReadyMosaicResourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    % Load the compute-ready MRGC mosaic
    load(fullfile(computeReadyMosaicResourcesDirectory, computeReadyMosaicFileName), 'theComputeReadyMRGCmosaic');

    % Generate and set the optics
    theComputeReadyMRGCmosaic.setTheOptics(opticsParamsForSubSpaceMapping);

    % Visualize the generated optics
    MosaicPoolingOptimizer.visualizeVlambdaWeightedPSF(theComputeReadyMRGCmosaic, opticsParamsForSubSpaceMapping);

    % Generate filename for the coneMosaic subspace responses
    [coneMosaicSubspaceResponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('coneMosaicSubspaceResponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParamsForSubSpaceMapping);

    % Generate filename for the mRGCMosaic subspace responses
    [mRGCMosaicSubspaceResponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSubspaceResponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParamsForSubSpaceMapping);

    % Ask the user what stimulus chromaticity to use
    % and update coneMosaicSubspaceResponsesFileName
    [stimulusChromaticity, coneMosaicSubspaceResponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
        coneMosaicSubspaceResponsesFileName, 'coneMosaicSubspaceResponses');
        
    % Update mRGCMosaicSubspaceResponsesFileName
    [~, mRGCMosaicSubspaceResponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
        mRGCMosaicSubspaceResponsesFileName, 'mRGCMosaicSubspaceResponses', ...
        'doNotQueryUserInsteadEmployThisStimulusChromaticity', stimulusChromaticity);

    % Optimally generated RF maps filename
    optimallyMappedSubspaceRFmapsFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', '_optimallyMappedRFs.mat');

    MosaicPoolingOptimizer.computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            maxSFcyclesPerDegree, stimSizeDegs, stimPositionDegs, ...
            fullfile(resourcesDirectory, coneMosaicSubspaceResponsesFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSubspaceResponsesFileName), ...
            fullfile(resourcesDirectory, optimallyMappedSubspaceRFmapsFileName), ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            visualizeOptimallyMappedRFmapLocations, ...
            'parPoolSize', parpoolSize);
end
