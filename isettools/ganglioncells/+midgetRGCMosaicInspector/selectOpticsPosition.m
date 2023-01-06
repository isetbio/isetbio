function  opticsPositionDegs = selectOpticsPosition(theMidgetRGCmosaic)

    opticsPositionDegs = [];
    while (numel(opticsPositionDegs) ~= 2)
        fprintf('\nSelect the desired optics position. For reference, in this midgetRGCMosaic, RTVFmodels were fitted at the following positions: ')
        for idx = 1:numel(theMidgetRGCmosaic.theConesNumPooledByTheRFcenterGrid)
            fprintf('\n[%+2.2f %+2.2f] (for %d-cone RF centers)', ...
                theMidgetRGCmosaic.theSamplingPositionGrid(idx,1), ...
                theMidgetRGCmosaic.theSamplingPositionGrid(idx,2), ...
                theMidgetRGCmosaic.theConesNumPooledByTheRFcenterGrid(idx));
        end

        opticsPositionDegs = input('\nEnter the optics position to use for computing retinal stimulus images ([x y]): ');
    end
end
