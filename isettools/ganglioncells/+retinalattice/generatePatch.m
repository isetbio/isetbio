function bestQualityRFpositions = generatePatch(fovDegs, neuronType, whichEye, exportHistoryToFile, visualizeConvergence, useParfor, maxIterations, varargin)
   
    p = inputParser;
    p.addParameter('randomSeed', [],  @(x)(isempty(x) || isscalar(x)));
    p.addParameter('customDegsToMMsConversionFunction', [], @(x) (isempty(x) || isa(x,'function_handle')));
    p.addParameter('customRFspacingFunction', [], @(x) (isempty(x) || isa(x,'function_handle')));
    p.addParameter('customMinRFspacing', [], @(x) (isempty(x) || isscalar(x)));
    p.parse(varargin{:});
    
    randomSeed = p.Results.randomSeed;
    customDegsToMMsConversionFunction = p.Results.customDegsToMMsConversionFunction;
    customRFspacingFunction = p.Results.customRFspacingFunction;
    customMinRFspacing = p.Results.customMinRFspacing;
    
    % Validate input
    validateInput(fovDegs, neuronType, whichEye);
    
    % Configure algorithm params
    params = retinalattice.configure(fovDegs, neuronType, whichEye);
    
    if (~isempty(customRFspacingFunction))
        params.rfSpacingExactFunction = customRFspacingFunction;
    end
    
    if (~isempty(customMinRFspacing))
        params.lambdaMinMicrons = customMinRFspacing;
    end
    
    if (~isempty(maxIterations))
        params.maxIterations = maxIterations;
    end
    
    if (~isempty(randomSeed))
        params.rng = randomSeed;
    end
    
    % Start timing
    tStart = tic;
    
    % Generate initial RF positions and downsample according to the density
    [rfPositionsMicrons, params.radius] = retinalattice.initialize(fovDegs, whichEye, params, useParfor, tStart, ...
        'customDegsToMMsConversionFunction', customDegsToMMsConversionFunction);
    
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
    fprintf('Lattice generation finished with status: ''%s'' with quality %2.4f at iteration %d.\n', dataOut.terminationReason, dataOut.qualityHistory(idx), idx);
    
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