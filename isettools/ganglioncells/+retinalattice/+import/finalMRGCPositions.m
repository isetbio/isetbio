function [rfPositionsMicrons, rfPositionsDegs] = finalMRGCPositions(sourceLatticeSizeDegs, eccDegs, sizeDegs, whichEye)

    % Convert degs to retinal microns
    eccMicrons  = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs);
    sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(sizeDegs, sqrt(sum(eccDegs.^2,2)));
    
    % Load mRGC RF positions
    p = retinalattice.configure(sourceLatticeSizeDegs, 'midget ganglion cells', whichEye);
    
    load(fullfile(p.latticeGalleryDir, p.patchFinalPositionsSaveFileName), 'rfPositions');
    rfPositionsMicrons = double(retinalattice.compute.croppedPositions(rfPositions, eccMicrons, sizeMicrons));
    
    % Convert positions to degs
    rfPositionsDegs = RGCmodels.Watson.convert.rhoMMsToDegs(rfPositionsMicrons*1e-3);
end

