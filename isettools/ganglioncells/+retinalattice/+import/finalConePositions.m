function [rfPositionsMicrons, rfPositionsDegs] = finalConePositions(sourceLatticeSizeDegs, eccDegs, sizeDegs, whichEye)

    % Convert degs to retinal microns
    eccMicrons  = 1000*RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs);
    sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(sizeDegs, sqrt(sum(eccDegs.^2,2)));

    % Load final cone positions
    p = retinalattice.configure(sourceLatticeSizeDegs, 'cones', whichEye);
    load(fullfile(p.latticeGalleryDir, p.patchFinalPositionsSaveFileName), 'rfPositions');
    rfPositionsMicrons = double(retinalattice.compute.croppedPositions(rfPositions, eccMicrons, sizeMicrons));
    
     % Convert positions to degs
    rfPositionsDegs = RGCmodels.Watson.convert.rhoMMsToDegs(rfPositionsMicrons*1e-3);
end

