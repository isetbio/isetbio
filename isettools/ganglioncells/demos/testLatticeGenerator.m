function testLatticeGenerator
    fovDegs = 45;
    whichEye = 'right eye';
    neuronType = 'cones'; %'midget ganglion cells';

    generateNewPatch = ~true;
    
    if (generateNewPatch)
        neuronType = 'cones';
        retinalattice.generatePatch(fovDegs, neuronType, whichEye);
        neuronType = 'midget ganglion cells';
        retinalattice.generatePatch(fovDegs, neuronType, whichEye);
    else
        retinalattice.savePositionsAtIteration(fovDegs, neuronType, whichEye);
        %retinalattice.inspectPatch(fovDegs, neuronType, whichEye);
    end
    
end

