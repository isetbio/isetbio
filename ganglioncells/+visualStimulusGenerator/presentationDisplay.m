%
%
%
function [theDisplay, backgroundChromaticity, backgroundLuminanceCdM2] = presentationDisplay(...
    wavelengthSupport, desiredPixelSizeDegs, viewingDistanceMeters, varargin)

    p = inputParser;
    p.addParameter('bitDepth', 20, @(x)(isscalar(x) && x > 0));
    p.addParameter('meanLuminanceCdPerM2', 50, @(x)(isscalar(x) && x > 0));
    p.addParameter('luminanceHeadroom', 0.1, @(x)(isscalar(x) && x >= 0 && x < 1));
    p.addParameter('displayType', '', @ischar);
    p.addParameter('adjustBackgroundChromaticityToEqualizeLandMconeExcitations', false, @islogical);
    p.addParameter('backgroundChromaticity', [], @(x)(isempty(x)||numel(x)==2));
    p.addParameter('backgroundLuminanceCdM2', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('coneFundamentalsToEmploy', []);
    p.parse(varargin{:});


    displayType = p.Results.displayType;
    if (isempty(displayType))
        displayType = 'LCD-Apple';
    end

    adjustBackgroundChromaticityToEqualizeLandMconeExcitations = p.Results.adjustBackgroundChromaticityToEqualizeLandMconeExcitations;
    backgroundChromaticity = p.Results.backgroundChromaticity;
    backgroundLuminanceCdM2 = p.Results.backgroundLuminanceCdM2;

    coneFundamentalsToEmploy = p.Results.coneFundamentalsToEmploy;

    % Generate display
    theDisplay = displayCreate(displayType, 'wave', wavelengthSupport);
    theDisplay = displaySet(theDisplay, 'viewing distance', viewingDistanceMeters);

    % Scale the white point so the requested mean luminance leaves the
    % requested fraction of display range above it.
    targetPeakLuminanceCdM2 = p.Results.meanLuminanceCdPerM2/(1-p.Results.luminanceHeadroom);
    theDisplay = localSetPeakLuminance(theDisplay, targetPeakLuminanceCdM2);

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
    if (adjustBackgroundChromaticityToEqualizeLandMconeExcitations) && ...
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

function theDisplay = localSetPeakLuminance(theDisplay, targetPeakLuminanceCdM2)

    currentPeakLuminanceCdM2 = displayGet(theDisplay, 'peak luminance');
    assert(currentPeakLuminanceCdM2 > 0, ...
        'Display %s has non-positive peak luminance.', displayGet(theDisplay, 'name'));

    spdScaleFactor = targetPeakLuminanceCdM2/currentPeakLuminanceCdM2;
    theDisplay = displaySet(theDisplay, 'spd', spdScaleFactor*displayGet(theDisplay, 'spd'));

end
