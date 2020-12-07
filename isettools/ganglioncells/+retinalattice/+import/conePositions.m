function [rfPositionsMicrons, rfPositionsDegs] = conePositions(sourceLatticeSizeDegs, eccDegs, sizeDegs, whichEye)

    % Convert degs to retinal microns
    eccMicrons  = 1000*RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs);
    sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(sizeDegs, sqrt(sum(eccDegs.^2,2)));

    % Load cone positions
    p = retinalattice.configure(sourceLatticeSizeDegs, 'cones', sprintf('%s eye',whichEye));
    load(fullfile(p.latticeGalleryDir, p.patchSaveFileName), 'dataOut');
    rfPositionsMicrons = dataOut.rfPositions;
    rfPositionsMicrons = retinalattice.compute.croppedPositions(rfPositionsMicrons, eccMicrons, sizeMicrons);
    
     % Convert positions to degs
    rfPositionsDegs = RGCmodels.Watson.convert.rhoMMsToDegs(rfPositionsMicrons*1e-3);
end

