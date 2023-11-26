function performComputeVisualRFsUsingMSequenceMapping(mosaicParams, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('stimSizeDegs', 1.0, @(x)(isscalar(x)||numel(x)==2));
    p.addParameter('stimPositionDegs', [], @(x)(isempty(x)||(numel(x) == 2)));
    p.addParameter('rfPixelsAcross', 16, @(x)(isscalar(x)&&(x>0)));
    p.addParameter('mSequenceBitLength', 10, @(x)(isscalar(x)&&(x>0)));
    p.addParameter('ternaryInsteadOfBinaryMsequence', false, @islogical);
    p.addParameter('reComputeInputConeMosaicResponses', false, @islogical);
    p.addParameter('reComputeMRGCMosaicResponses', false, @islogical);
    p.addParameter('reComputeRFs', false, @islogical);
    p.addParameter('visualizeOptimallyMappedRFmapLocations', false, @islogical);
    p.parse(varargin{:});


    % stimulus patch size
    stimSizeDegs = p.Results.stimSizeDegs;    
    
    % simulus position
    stimPositionDegs = p.Results.stimPositionDegs;

    % pixels (x,y)
    rfPixelsAcross = p.Results.rfPixelsAcross;

    % Bit length of m-sequence
    mSequenceBitLength = p.Results.mSequenceBitLength;

    % Type of msequence
    ternaryInsteadOfBinaryMsequence = p.Results.ternaryInsteadOfBinaryMsequence;

    % What to compute
    reComputeInputConeMosaicResponses = p.Results.reComputeInputConeMosaicResponses;
    reComputeMRGCMosaicResponses = p.Results.reComputeMRGCMosaicResponses;
    reComputeRFs = p.Results.reComputeRFs;
    visualizeOptimallyMappedRFmapLocations = p.Results.visualizeOptimallyMappedRFmapLocations;

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
    fprintf('\n---> Select the optics to use for m-sequence RF mapping\n');
    [opticsParamsForMSequenceMapping, opticsToEmploy] = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Generate and set the optics
    theComputeReadyMRGCmosaic.setTheOptics(opticsParamsForMSequenceMapping);

    % Visualize the generated optics
    MosaicPoolingOptimizer.visualizeVlambdaWeightedPSF(theComputeReadyMRGCmosaic, opticsParamsForMSequenceMapping);

    % Generate filename for the coneMosaic subspace responses
    [coneMosaicMSequenceResponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('coneMosaicMSequenceResponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParamsForMSequenceMapping);

    % Generate filename for the mRGCMosaic subspace responses
    [mRGCMosaicMSequenceResponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicMSequenceResponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParamsForMSequenceMapping);

    % Ask the user what stimulus chromaticity to use
    % and update mRGCMosaicMSequenceResponsesFileName
    [stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, mRGCMosaicMSequenceResponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
        mRGCMosaicMSequenceResponsesFileName, 'mRGCMosaicMSequenceResponses');
    

    % Update the coneMosaicMSequenceResponsesFileName
    [~, ~, coneMosaicMSequenceResponsesFileName] = ...
        MosaicPoolingOptimizer.chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
            coneMosaicMSequenceResponsesFileName, 'coneMosaicMSequenceResponses', ...
            'doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals', struct(...
                    'stimulusChromaticity', stimulusChromaticity,...
                    'coneFundamentalsOptimizedForStimPosition', coneFundamentalsOptimizedForStimPosition));


    % Encode m-sequence type in responses filename
    if (ternaryInsteadOfBinaryMsequence)
        coneMosaicMSequenceResponsesFileName = strrep(coneMosaicMSequenceResponsesFileName, 'MSequence', sprintf('Ternary%2.0fBitMSequence', mSequenceBitLength));
        mRGCMosaicMSequenceResponsesFileName = strrep(mRGCMosaicMSequenceResponsesFileName, 'MSequence', sprintf('Ternary%2.0fBitMSequence', mSequenceBitLength));
    else
        coneMosaicMSequenceResponsesFileName = strrep(coneMosaicMSequenceResponsesFileName, 'MSequence', sprintf('Binary%2.0fBitMSequence', mSequenceBitLength));
        mRGCMosaicMSequenceResponsesFileName = strrep(mRGCMosaicMSequenceResponsesFileName, 'MSequence', sprintf('Binary%2.0fBitMSequence', mSequenceBitLength));
    end

    % Encode rfPixelsAcross
    coneMosaicMSequenceResponsesFileName = strrep(coneMosaicMSequenceResponsesFileName, 'MSequence', sprintf('MSequence%dx%dPixels', rfPixelsAcross, rfPixelsAcross));
    mRGCMosaicMSequenceResponsesFileName = strrep(mRGCMosaicMSequenceResponsesFileName, 'MSequence', sprintf('MSequence%dx%dPixels', rfPixelsAcross, rfPixelsAcross));


    % Optimally generated RF maps filename
    optimallyMappedMSequenceRFmapsFileName = strrep(mRGCMosaicMSequenceResponsesFileName, '.mat', '_optimallyMappedRFs.mat');


    MosaicPoolingOptimizer.computeVisualRFsUsingMSequenceMapping(...
        theComputeReadyMRGCmosaic, opticsToEmploy, ...
        stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
        ternaryInsteadOfBinaryMsequence, mSequenceBitLength, ...
        stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
        fullfile(resourcesDirectory, coneMosaicMSequenceResponsesFileName), ...
        fullfile(resourcesDirectory, mRGCMosaicMSequenceResponsesFileName), ...
        fullfile(resourcesDirectory, optimallyMappedMSequenceRFmapsFileName), ...
        reComputeInputConeMosaicResponses, ...
        reComputeMRGCMosaicResponses, ...
        reComputeRFs, ...
        visualizeOptimallyMappedRFmapLocations);

end
