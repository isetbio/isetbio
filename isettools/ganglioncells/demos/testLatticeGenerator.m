function testLatticeGenerator
    fovDegs = 2;
    whichEye = 'right eye';
    neuronType = 'midget ganglion cells';

    generateNewPatch = ~true;
    
    if (generateNewPatch)
        neuronType = 'cones';
        retinalattice.generatePatch(fovDegs, neuronType, whichEye);
        neuronType = 'midget ganglion cells';
        retinalattice.generatePatch(fovDegs, neuronType, whichEye);
    else
        retinalattice.inspectPatch(fovDegs, neuronType, whichEye);
    end
    
end

