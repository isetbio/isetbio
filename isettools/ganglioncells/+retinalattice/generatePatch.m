function generatePatch(fovDegs, neuronType, whichEye)

    % Validate input
    validateInput(fovDegs, neuronType, whichEye);
    
    % Configure algorithm params
    params = retinalattice.configure(fovDegs, neuronType, whichEye);
    
    % Start timing
    tStart = tic;
    
    % Generate initial RF positions and downsample according to the density
    [rfPositionsMicrons, params.radius] = retinalattice.initialize(fovDegs, whichEye, params, tStart);
    
    % Compute RF spacing look-up tables
    [tabulatedRFspacingMicrons, tabulatedEccMicrons] = retinalattice.compute.rfSpacingLookUpTables(...
        rfPositionsMicrons, whichEye, params.rfSpacingExactFunction, params.eccentricityLookUpTableSamplesNum);
    
    % Iteratively smooth the grid of RFs
    dataOut = retinalattice.compute.smoothGrid(rfPositionsMicrons, tabulatedEccMicrons, tabulatedRFspacingMicrons, ...
            params, tStart);
    
    % Save lattice generation data
    save(fullfile(params.latticeGalleryDir ,params.patchSaveFileName), 'dataOut', 'params', 'fovDegs', 'neuronType', 'whichEye', '-v7.3');
    fprintf('Lattice generation finished with status: ''%s''.\n', dataOut.terminationReason);
    fprintf('Lattice data for %d nodes and history saved  in ''%s'' (gallery dir: ''%s'').\n', ...
        size(dataOut.rfPositions,1), params.patchSaveFileName, params.latticeGalleryDir);
    
end


function validateInput(fovDegs, neuronType, whichEye)
    validNeuronTypes = retinalattice.validvalues.neuronTypes;
    validEyes = retinalattice.validvalues.eyes; 
    assert(ismember(neuronType, validNeuronTypes), sprintf('Unknown neuron type: ''%s''.', neuronType));
    assert(ismember(whichEye, validEyes), sprintf('Unknown eye: ''%s''.', whichEye));
    assert(isscalar(fovDegs), sprintf('fovDegs must be a scalar'));
end