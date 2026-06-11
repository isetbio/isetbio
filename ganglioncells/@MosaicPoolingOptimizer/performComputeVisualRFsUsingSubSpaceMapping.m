function performComputeVisualRFsUsingSubSpaceMapping(mosaicParams, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('stimSizeDegs', 1.0, @(x)(isscalar(x)||numel(x)==2));
    p.addParameter('stimPositionDegs', [], @(x)(isempty(x)||(numel(x) == 2)));
    p.addParameter('maxSFLimit', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('maxSFToBeAnalyzed', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('rfMappingPixelMagnificationFactor', 1, @(x)(isscalar(x)&&(x>=1)));
    p.addParameter('reComputeInputConeMosaicResponses', false, @islogical);
    p.addParameter('reComputeMRGCMosaicResponses', false, @islogical);
    p.addParameter('reComputeRFs', false, @islogical);
    p.addParameter('visualizeOptimallyMappedRFmapLocations', false, @islogical);
    p.addParameter('visualizedRGCindex', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('msequencePixelSizeDegs', 0.01, @isscalar);
    p.parse(varargin{:});

    % stimulus patch size
    stimSizeDegs = p.Results.stimSizeDegs;    
    
    % simulus position
    stimPositionDegs = p.Results.stimPositionDegs;

    % stimulus max spatial frequency (can be set to lower that optimal retinal sf)
    maxSFLimit = p.Results.maxSFLimit;

    % max spatial frequency to include in the RFmap
    maxSFToBeAnalyzed = p.Results.maxSFToBeAnalyzed;

    % pixel magnification factor >= 1
    rfMappingPixelMagnificationFactor = p.Results.rfMappingPixelMagnificationFactor;

    % Pixel size for the simulation of m-sequence RF map
    msequencePixelSizeDegs = p.Results.msequencePixelSizeDegs;

    % What to compute
    reComputeInputConeMosaicResponses = p.Results.reComputeInputConeMosaicResponses;
    reComputeMRGCMosaicResponses = p.Results.reComputeMRGCMosaicResponses;
    reComputeRFs = p.Results.reComputeRFs;
    visualizeOptimallyMappedRFmapLocations = p.Results.visualizeOptimallyMappedRFmapLocations;
    
    visualizedRGCindex = p.Results.visualizedRGCindex;

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

   
    % Ask the user which optics to use
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
    [stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, mRGCMosaicSubspaceResponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
        mRGCMosaicSubspaceResponsesFileName, 'mRGCMosaicSubspaceResponses');
    

    % Update the coneMosaicSumspaceResponsesFileName 
    [~, ~, coneMosaicSubspaceResponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
            coneMosaicSubspaceResponsesFileName, 'coneMosaicSubspaceResponses', ...
            'doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals', struct(...
                    'stimulusChromaticity', stimulusChromaticity,...
                    'coneFundamentalsOptimizedForStimPosition', coneFundamentalsOptimizedForStimPosition));


    if (isempty(maxSFToBeAnalyzed)) || (maxSFToBeAnalyzed>maxSFLimit)
        maxSFToBeAnalyzed = maxSFLimit;
    end

    % Optimally generated RF maps filename
    optimallyMappedSubspaceRFmapsFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', '_optimallyMappedRFs.mat');

    if (~isempty(maxSFToBeAnalyzed))
        optimallyMappedSubspaceRFmapsFileName = strrep(optimallyMappedSubspaceRFmapsFileName, '.mat', sprintf('_LimitedTo%2.2fCPD.mat', maxSFToBeAnalyzed));
    end

    MosaicPoolingOptimizer.computeVisualRFsUsingSubSpaceMapping(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, ...
            maxSFLimit, maxSFToBeAnalyzed, rfMappingPixelMagnificationFactor, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            fullfile(resourcesDirectory, coneMosaicSubspaceResponsesFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSubspaceResponsesFileName), ...
            fullfile(resourcesDirectory, optimallyMappedSubspaceRFmapsFileName), ...
            reComputeInputConeMosaicResponses, ...
            reComputeMRGCMosaicResponses, ...
            reComputeRFs, ...
            visualizeOptimallyMappedRFmapLocations, ...
            'msequencePixelSizeDegs', msequencePixelSizeDegs, ...
            'visualizedRGCindex', visualizedRGCindex, ...
            'parPoolSize', parpoolSize);
end
