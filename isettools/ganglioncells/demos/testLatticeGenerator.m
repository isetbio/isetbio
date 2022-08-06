% Function to either generate a lattice or save the lattice at some point
% in its progression.
function testLatticeGenerator
    fovDegs = 60;
    exportHistoryToFile = true;
    visualizeConvergence = true; 
    useParfor = true;
    maxIterations = 8000;
    neuronType = 'midget ganglion cells'; % Select from { 'cones', 'midget ganglion cells'};

    generateNewPatch = ~true;
    
    if (generateNewPatch)
        whichEye = 'right eye';
        retinalattice.generatePatch(fovDegs, neuronType, whichEye, ...
            exportHistoryToFile, visualizeConvergence, useParfor, maxIterations, ...
                'eccentricityLookUpTableSamplesNum', 256);

        whichEye = 'left eye';
        retinalattice.generatePatch(fovDegs, neuronType, whichEye, ...
            exportHistoryToFile, visualizeConvergence, useParfor, maxIterations, ...
             'eccentricityLookUpTableSamplesNum', 256);

        %neuronType = 'midget ganglion cells';
        %retinalattice.generatePatch(fovDegs, neuronType, whichEye);
    else
        whichEye = 'left eye';
        retinalattice.savePositionsAtIteration(fovDegs, neuronType, whichEye);
        %retinalattice.inspectPatch(fovDegs, neuronType, whichEye);
    end
    
end

