% Function to either generate a lattice or save the lattice at some point
% in its progression.
function testLatticeGenerator
    fovDegs = 58;
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

