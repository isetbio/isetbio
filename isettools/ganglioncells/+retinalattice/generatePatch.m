function bestQualityRFpositions = generatePatch(fovDegs, neuronType, whichEye, exportHistoryToFile, visualizeConvergence, useParfor, maxIterations)
   
    % Validate input
    validateInput(fovDegs, neuronType, whichEye);
    
    % Configure algorithm params
    params = retinalattice.configure(fovDegs, neuronType, whichEye);
    
    if (~isempty(maxIterations))
        params.maxIterations = maxIterations;
    end
    
    % Start timing
    tStart = tic;
    
    % Generate initial RF positions and downsample according to the density
    [rfPositionsMicrons, params.radius] = retinalattice.initialize(fovDegs, whichEye, params, useParfor, tStart);
    
    % Compute RF spacing look-up tables
    [tabulatedRFspacingMicrons, tabulatedEccMicrons] = retinalattice.compute.rfSpacingLookUpTables(...
        rfPositionsMicrons, whichEye, useParfor, params.rfSpacingExactFunction, params.eccentricityLookUpTableSamplesNum);
    
    % Iteratively smooth the grid of RFs
    dataOut = retinalattice.compute.smoothGrid(rfPositionsMicrons, tabulatedEccMicrons, tabulatedRFspacingMicrons, ...
            params, visualizeConvergence, tStart);
        
    % Return positions corresponding to max lattice quality
    [~, idx] = max(dataOut.qualityHistory);
    bestQualityRFpositions = double(squeeze(dataOut.rfPositionsHistory(idx,:,:)));
    
    % Report back
    fprintf('Lattice generation finished with status: ''%s''.\n', dataOut.terminationReason);
    
    % Save lattice generation data
    if (exportHistoryToFile)
        save(fullfile(params.latticeGalleryDir ,params.patchSaveFileName), 'dataOut', 'params', 'fovDegs', 'neuronType', 'whichEye', '-v7.3');
        fprintf('Lattice data for %d nodes and history saved  in ''%s'' (gallery dir: ''%s'').\n', ...
            size(dataOut.rfPositions,1), params.patchSaveFileName, params.latticeGalleryDir);
    end
    
end


function validateInput(fovDegs, neuronType, whichEye)
    validNeuronTypes = retinalattice.validvalues.neuronTypes;
    validEyes = retinalattice.validvalues.eyes; 
    assert(ismember(neuronType, validNeuronTypes), sprintf('Unknown neuron type: ''%s''.', neuronType));
    assert(ismember(whichEye, validEyes), sprintf('Unknown eye: ''%s''.', whichEye));
    assert(isscalar(fovDegs), sprintf('fovDegs must be a scalar'));
end