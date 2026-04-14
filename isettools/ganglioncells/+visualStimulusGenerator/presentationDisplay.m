%
%
%
function [theDisplay, backgroundChromaticity, backgroundLuminanceCdM2] = presentationDisplay(...
    wavelengthSupport, desiredPixelSizeDegs, viewingDistanceMeters, varargin)

    p = inputParser;
    p.addParameter('bitDepth', 20, @isscalar);
    p.addParameter('meanLuminanceCdPerM2', 50, @isscalar);
    p.addParameter('luminanceHeadroom', 0.1, @isscalar);
    p.addParameter('displayType', '', @ischar);
    p.addParameter('adjustBackgroundChromaticityToEqualizeLandMconeExcitations', false, @islogical);
    p.addParameter('backgroundChromaticity', [], @(x)(isempty(x)||numel(x)==2));
    p.addParameter('backgroundLuminanceCdM2', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('coneFundamentalsToEmploy', []);
    p.parse(varargin{:});


    displayParams = generateConventionalxyYDisplayDefaultParams;

    if (~isempty(p.Results.displayType))
        displayParams.whichDisplay = p.Results.displayType;
    end
    displayParams.viewingDistanceMeters = viewingDistanceMeters;
    displayParams.spectralSupport = wavelengthSupport;
    displayParams.meanLuminanceCdPerM2 = p.Results.meanLuminanceCdPerM2;
    displayParams.luminanceHeadroom = p.Results.luminanceHeadroom;

    adjustBackgroundChromaticityToEqualizeLandMconeExcitations = p.Results.adjustBackgroundChromaticityToEqualizeLandMconeExcitations;
    backgroundChromaticity = p.Results.backgroundChromaticity;
    backgroundLuminanceCdM2 = p.Results.backgroundLuminanceCdM2;

    coneFundamentalsToEmploy = p.Results.coneFundamentalsToEmploy;

    % Generate display
    theDisplay = generateConventionalxyYDisplay(displayParams);

    % Linear LUT
    bitDepth = p.Results.bitDepth;
    N = 2^bitDepth;
    gTable = repmat(linspace(0, 1, N), 3, 1)';
    theDisplay = displaySet(theDisplay, 'gTable', gTable);

    % Correct the dpi so we end up with the desired pixel size (in visual degrees)
    
    % 1. Get current pixel size in degs
    pixelSizeDegs = displayGet(theDisplay, 'degperpixel');
    scaleFactorToMatchDesiredPixelSizeDegs = pixelSizeDegs/desiredPixelSizeDegs;
    
    % 3. original dots per inch
    dpiOriginal = displayGet(theDisplay, 'dpi');
    dpiDesired = dpiOriginal * scaleFactorToMatchDesiredPixelSizeDegs;
    
    % 4. Set desired dots per inch
    theDisplay = displaySet(theDisplay, 'dpi', dpiDesired);

    % 5. Adjust background chromaticity to enable equal L and M cone excitations
    if (~isempty(adjustBackgroundChromaticityToEqualizeLandMconeExcitations)) && ...
       (~isempty(backgroundChromaticity)) && (~isempty(backgroundLuminanceCdM2)) 
        beforexyY = backgroundChromaticity;
        beforexyY(3) = backgroundLuminanceCdM2;
        [backgroundChromaticity, backgroundLuminanceCdM2] = ...
            RGCMosaicConstructor.helper.simulateExperiment.updateBackgroundToAchieveEqualLandMconeActivation(...
                    theDisplay, backgroundChromaticity, backgroundLuminanceCdM2, ...
                    'coneFundamentalsToEmploy', coneFundamentalsToEmploy);

        fprintf(2,'Adjusted background chromaticity from (%2.2f, %2.2f, lum = %2.1f cd/m2) to (%2.2f,%2.2f, lum = %2.1f cd/m2) to achieve equal L and M cone excitations\n', ...
            beforexyY(1), beforexyY(2), beforexyY(3), backgroundChromaticity(1), backgroundChromaticity(2), backgroundLuminanceCdM2)
    end
end