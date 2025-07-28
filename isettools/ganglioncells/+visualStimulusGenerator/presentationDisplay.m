function theDisplay = presentationDisplay(wavelengthSupport, desiredPixelSizeDegs, viewingDistanceMeters)

    theDisplay = displayCreate('LCD-Apple', ...
        'wave', wavelengthSupport, ...
        'viewing distance',viewingDistanceMeters);
    
    % Linear, 12-bit LUT
    bitDepth = 12;
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