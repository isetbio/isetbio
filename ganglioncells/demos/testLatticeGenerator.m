% Function to either generate a lattice or save the lattice at some point
% in its progression.
function testLatticeGenerator
    fovDegs = 64;
    exportHistoryToFile = true;
    visualizeConvergence = true; 
    useParfor = true;
    maxIterations = 8000;
    neuronType = 'cones'; % midget ganglion cells'; % Select from { 'cones', 'midget ganglion cells'};

    generateNewPatch = true;
    generateLeftEyeMosaic = true;
    generateRightEyeMosaic = ~true;

    if (generateNewPatch)
        if (generateRightEyeMosaic)
            whichEye = 'right eye';
            retinalattice.generatePatch(fovDegs, neuronType, whichEye, ...
                exportHistoryToFile, visualizeConvergence, useParfor, maxIterations, ...
                    'eccentricityLookUpTableSamplesNum', 256);
        end

        if (generateLeftEyeMosaic)
            whichEye = 'left eye';
            retinalattice.generatePatch(fovDegs, neuronType, whichEye, ...
                exportHistoryToFile, visualizeConvergence, useParfor, maxIterations, ...
                 'eccentricityLookUpTableSamplesNum', 256);
        end

        %neuronType = 'midget ganglion cells';
        %retinalattice.generatePatch(fovDegs, neuronType, whichEye);
    else
        if (generateLeftEyeMosaic)
            whichEye = 'left eye';
            retinalattice.savePositionsAtIteration(fovDegs, neuronType, whichEye);
            %retinalattice.inspectPatch(fovDegs, neuronType, whichEye);
        end

        if (generateRightEyeMosaic)
            whichEye = 'right eye';
            retinalattice.savePositionsAtIteration(fovDegs, neuronType, whichEye);
            %retinalattice.inspectPatch(fovDegs, neuronType, whichEye);
        end

    end
    
end

