function  opticsPositionDegs = selectOpticsPosition(theMidgetRGCmosaic)

    opticsPositionDegs = [];
    while (numel(opticsPositionDegs) ~= 2)
        fprintf('\nSelect the desired optics position. For reference, in this midgetRGCMosaic, RTVFmodels were fitted at the following positions: ')
        centerConesNumExamined = sort(unique(theMidgetRGCmosaic.theConesNumPooledByTheRFcenterGrid), 'ascend');
        for iCone = 1:numel(centerConesNumExamined)
            iidx = find(theMidgetRGCmosaic.theConesNumPooledByTheRFcenterGrid == centerConesNumExamined(iCone));
            for idx = 1:numel(iidx)
                iObj = iidx(idx);
                fprintf('\n[%+2.2f %+2.2f] (for %d-cone RF centers)', ...
                    theMidgetRGCmosaic.theSamplingPositionGrid(iObj,1), ...
                    theMidgetRGCmosaic.theSamplingPositionGrid(iObj,2), ...
                    theMidgetRGCmosaic.theConesNumPooledByTheRFcenterGrid(iObj));
            end
        end

        opticsPositionDegs = input('\nEnter the optics position to use for computing retinal stimulus images ([x y]): ');
    end
end
