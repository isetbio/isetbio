function testLatticeGenerator
    fovDegs = 2;
    neuronType = 'midget ganglion cells';
    whichEye = 'right eye';
    patchSaveFileName = sprintf('%s_%s_%1.0fdeg_mosaic_progress', ...
        strrep(whichEye, ' ', '_'), strrep(neuronType, ' ', '_'), fovDegs);
    
    generateNewPatch = ~false;
    
    if (generateNewPatch)
        retinalattice.generatePatch(fovDegs, neuronType, whichEye, patchSaveFileName);
    else
        retinalattice.inspectPatch(patchSaveFileName);
    end
    
end

