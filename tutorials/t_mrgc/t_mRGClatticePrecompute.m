function t_mRGClatticePrecompute

    obj.eccentricityDegs = [0 0];
    obj.sizeDegs = 45*[1 1];
    obj.useParfor = true;
    
    maxIterations = 500;
    
    % Regenerate lattice whose FOV is large enough to encopass the desired size at the desired eccentricity
    fovDegs = sqrt(sum(obj.eccentricityDegs.^2,2)) + max(obj.sizeDegs)*1.3;

    obj.whichEye = 'left eye';
    generatePatch(obj, maxIterations, fovDegs);
    
    obj.whichEye = 'right eye';
    generatePatch(obj, maxIterations, fovDegs);
    
end

function obj = generatePatch(obj, maxIterations, fovDegs)
    exportHistoryToFile = true;
    visualizeConvergence = false;
    
    retinalattice.generatePatch(fovDegs, ...
        'midget ganglion cells', obj.whichEye, exportHistoryToFile, visualizeConvergence, obj.useParfor, maxIterations);
    
end
