function presentationDisplay = generateCRTDisplay()
    % Load a CRT display, d
    load('CRT-MODEL.mat', 'd');
    presentationDisplay = d;
    clear 'd';

    % Linear, 12-bit LUT
    bitDepth = 12;
    N = 2^bitDepth;
    gTable = repmat(linspace(0, 1, N), 3, 1)';

    presentationDisplay = displaySet(presentationDisplay, 'gTable', gTable);
end