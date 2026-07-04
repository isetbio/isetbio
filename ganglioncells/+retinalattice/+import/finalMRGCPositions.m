function [rfPositionsMicrons, rfPositionsDegs] = finalMRGCPositions(...
    sourceLatticeSizeDegs, eccMicrons, sizeMicrons, whichEye, ...
    MMsToDegsConversionFunction)

    % Load mRGC RF positions
    p = retinalattice.configure(sourceLatticeSizeDegs, 'midget ganglion cells', whichEye);
    
    load(fullfile(p.latticeGalleryDir, p.patchFinalPositionsSaveFileName), 'rfPositions');
    
    % Reverse the polarity
    rfPositions = -rfPositions;
    
    rfPositionsMicrons = double(retinalattice.compute.croppedPositions(rfPositions, eccMicrons, sizeMicrons));
    
    % Convert positions to degs
    rfPositionsDegs = MMsToDegsConversionFunction(rfPositionsMicrons*1e-3);
end

