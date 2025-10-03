function [stimParams, thePresentationDisplay] = setupSubspaceRFmappingExperiment(wavelengthSupport, ...
        stimSizeDegs, optimalRetinalPixelSizeDegs, maxSFLimit, stimulusChromaticity)

    % At least 7 samples / period
    minNumberOfPixelsPerSpatialPeriod = 7;
    maxSFretinal = 1/(minNumberOfPixelsPerSpatialPeriod*sqrt(2.0)*optimalRetinalPixelSizeDegs);

    % Maybe cap the maxSF to reduce the # of stimuli
    if (~isempty(maxSFLimit)) && (maxSFLimit < maxSFretinal)
        maxSF = maxSFLimit;
        fprintf('Optimal SF: %2.1f c/deg. Capping max SF to %2.1f c/deg.\n', maxSFretinal, maxSFLimit);
    else
        maxSF = maxSFretinal;
        fprintf('Optimal SF: %2.1f c/deg. Max SF set to this value.\n', maxSFretinal);
    end

    pixelSizeDegs = 1/(minNumberOfPixelsPerSpatialPeriod*maxSF*sqrt(2.0));


    % Generate a presentation display with a desired resolution
    stimulusPixelsNum = round(max(stimSizeDegs(:))/pixelSizeDegs);
    
    omega = round(maxSF * max(stimSizeDegs));
    nStim = (2*omega+1)^2;

    a = single(1);
    s = whos('a');
    HartleyMemoryRequirementGBytes = (nStim * stimulusPixelsNum * stimulusPixelsNum * s.bytes)/1024/1024/1024;
    fprintf('\nTo probe RFs with spatial frequencies up to %2.1f c/deg\nusing a patch size of %2.1f degs,\n%d (%d x %d) Hartley patterns (retinal res:%2.3f arc min, optimal retinal res: %2.3f arc min) will be employed.\nThis requires %2.1f GBytes or RAM\n', ...
        maxSF, stimSizeDegs(1), nStim, stimulusPixelsNum, stimulusPixelsNum, pixelSizeDegs*60, optimalRetinalPixelSizeDegs*60, HartleyMemoryRequirementGBytes);

    % Generate a presentation display with a desired resolution
    viewingDistanceMeters = 4;

    thePresentationDisplay = rfMappingStimulusGenerator.presentationDisplay(...
            wavelengthSupport, pixelSizeDegs, ...
            viewingDistanceMeters);

    % Cone contrasts and overall contrast for desired stimulus chromaticity
    [coneContrasts, contrast] = MosaicPoolingOptimizer.contrastForChromaticity(stimulusChromaticity);

    % Stim params for the RF mapping
    stimParams = struct(...
            'backgroundLuminanceCdM2', 50.0, ...
            'backgroundChromaticity', [0.301 0.301], ...
            'coneContrasts', coneContrasts, ...
            'contrast', contrast, ...
            'omega', omega, ...
            'maxSF', maxSF, ...
            'pixelSizeDegs', pixelSizeDegs, ...
            'stimSizeDegs', max(stimSizeDegs), ...
            'wavelengthSupport', displayGet(thePresentationDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(thePresentationDisplay, 'viewing distance') ...
            );

end
