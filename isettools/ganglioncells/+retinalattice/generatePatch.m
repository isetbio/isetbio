function generatePatch(fovDegs, neuronType, whichEye, patchSaveFileName)
    % Parse input
    retinalattice.parseInput(neuronType, whichEye);
    
    % Configure algorithm params
    params = retinalattice.configure(neuronType);
    
    % Start timing
    tStart = tic;
    
    % Generate initial RF positions and downsample according to the density
    [rfPositionsMicrons, params.radius] = retinalattice.initialize(fovDegs, whichEye, params, tStart);
    
    % Compute RF spacing look-up tables
    [tabulatedRFspacingMicrons, tabulatedEccMicrons] = retinalattice.rfSpacingLookUpTables(...
        rfPositionsMicrons, whichEye, params.rfSpacingExactFunction, params.eccentricityLookUpTableSamplesNum);
    
    % Iteratively smooth the grid of RFs
    dataOut = retinalattice.smoothGrid(rfPositionsMicrons, tabulatedEccMicrons, tabulatedRFspacingMicrons, ...
            params, tStart, patchSaveFileName);
    
    % Save lattice generation data
    save(patchSaveFileName, 'dataOut', 'params', 'fovDegs', 'neuronType', 'whichEye', '-v7.3');
    fprintf('History saved  in %s\n', patchSaveFileName);
    
end


