function performComputeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic(mosaicParams, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('employConeFundamentalsDerivedFromInputConeMosaicAtStimPosition', false, @islogical);
    p.addParameter('stimSizeDegs', 1.0, @(x)(isscalar(x)||numel(x)==2));
    p.addParameter('stimPositionDegs', [], @(x)(isempty(x)||(numel(x) == 2)));
    p.addParameter('maxSFLimit', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('rfMappingPixelMagnificationFactor', 1, @(x)(isscalar(x)&&(x>=1)));
    p.addParameter('reComputeInputConeMosaicSubspaceRFmappingResponses', false, @islogical);
    p.addParameter('reComputeMRGCMosaicSubspaceRFmappingResponses', false, @islogical);
    p.addParameter('reComputeRFs', false, @islogical);
    p.addParameter('visualizeOptimallyMappedRFmapLocations', false, @islogical);
    p.parse(varargin{:});


    % Whether to emply cone fundamentals appropriate for the position
    % within the cone mosaic that we are stimulating and the current optics (lens, really)
    employConeFundamentalsDerivedFromInputConeMosaicAtStimPosition = p.Results.employConeFundamentalsDerivedFromInputConeMosaicAtStimPosition;

    % stimulus patch size
    stimSizeDegs = p.Results.stimSizeDegs;    
    
    % simulus position
    stimPositionDegs = p.Results.stimPositionDegs;

    % stimulus max spatial frequency (can be set to lower that optimal retinal sf)
    maxSFLimit = p.Results.maxSFLimit;

    % pixel magnification factor >= 1
    rfMappingPixelMagnificationFactor = p.Results.rfMappingPixelMagnificationFactor;

    % Which grid nodes to compute RFs
    reComputeInputConeMosaicSubspaceRFmappingResponses = p.Results.reComputeInputConeMosaicSubspaceRFmappingResponses;
    reComputeMRGCMosaicSubspaceRFmappingResponses = p.Results.reComputeMRGCMosaicSubspaceRFmappingResponses;
    reComputeRFs = p.Results.reComputeRFs;
    visualizeOptimallyMappedRFmapLocations = p.Results.visualizeOptimallyMappedRFmapLocations;

    % Set the parpoolsize to [] to do 
    parpoolSize = [];

    % Set the parpoolize to 1 to do one stimulus frame at a time so we can 
    % visualize cone mosaic responses to each frame
    %parpoolSize = 1;

    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    fprintf('\n---> Select the optics that were used to compute the compute-ready mosaic\n');
    opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);


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

   
    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
    fprintf('\n---> Select the optics to use for subspace RF mapping\n');
    [opticsParamsForSubSpaceMapping, opticsToEmploy] = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

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
    [stimulusChromaticity, mRGCMosaicSubspaceResponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
        mRGCMosaicSubspaceResponsesFileName, 'mRGCMosaicSubspaceResponses');
    

    % Update the coneMosaicSumspaceResponsesFileName
    [~, coneMosaicSubspaceResponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
            coneMosaicSubspaceResponsesFileName, 'coneMosaicSubspaceResponses', ...
            'doNotQueryUserInsteadEmployThisStimulusChromaticity', stimulusChromaticity);


    % Optimally generated RF maps filename
    optimallyMappedSubspaceRFmapsFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', '_optimallyMappedRFs.mat');

    MosaicPoolingOptimizer.computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, stimulusChromaticity, ...
            employConeFundamentalsDerivedFromInputConeMosaicAtStimPosition, ...
            maxSFLimit, rfMappingPixelMagnificationFactor, ...
            fullfile(resourcesDirectory, coneMosaicSubspaceResponsesFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSubspaceResponsesFileName), ...
            fullfile(resourcesDirectory, optimallyMappedSubspaceRFmapsFileName), ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            visualizeOptimallyMappedRFmapLocations, ...
            'parPoolSize', parpoolSize);
end
