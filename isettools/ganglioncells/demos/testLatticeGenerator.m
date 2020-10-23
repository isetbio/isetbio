function testLatticeGenerator
    fovDegs = 2;
    neuronType = 'midget ganglion cells';
    whichEye = 'right eye';
    
    patchSaveFileName = 'coneLeftEyePatch';
    generateNewPatch = ~false;
    
    if (generateNewPatch)
        retinalattice.generatePatch(fovDegs, neuronType, whichEye, patchSaveFileName);
    else
        retinalattice.inspectPatch(patchSaveFileName);
    end
    
end

