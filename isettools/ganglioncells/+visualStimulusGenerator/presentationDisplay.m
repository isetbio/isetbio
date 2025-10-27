%
%
%
function theDisplay = presentationDisplay(...
    wavelengthSupport, desiredPixelSizeDegs, viewingDistanceMeters, varargin)

    p = inputParser;
    p.addParameter('bitDepth', 20, @isscalar);
    p.addParameter('meanLuminanceCdPerM2', 50, @isscalar);
    p.addParameter('luminanceHeadroom', 0.1, @isscalar);
    p.parse(varargin{:});


    displayParams = generateConventionalxyYDisplayDefaultParams;

    displayParams.viewingDistanceMeters = viewingDistanceMeters;
    displayParams.spectralSupport = wavelengthSupport;
    displayParams.meanLuminanceCdPerM2 = p.Results.meanLuminanceCdPerM2;
    displayParams.luminanceHeadroom = p.Results.luminanceHeadroom;
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
end