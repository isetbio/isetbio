function mappedSizeDegs = mSequenceRFmapSizeDegsForMosaic(theMRGCMosaic, mappedPositionDegs, targetNumberOfMappedCells)

    mosaicEccDegs = theMRGCMosaic.eccentricityDegs;
    mosaicRadialEccDegs = sqrt(sum(mosaicEccDegs.^2));

    if (isempty(targetNumberOfMappedCells))
            % Some reasonable sizes
            if (mosaicRadialEccDegs < 1)
                mappedSizeDegs = 0.2*[1 1];
            elseif (mosaicRadialEccDegs <= 2)
                mappedSizeDegs = 0.3*[1 1];
            elseif (mosaicRadialEccDegs <= 3)
                mappedSizeDegs = 0.5*[1 1];
            elseif (mosaicRadialEccDegs <= 4)
                mappedSizeDegs = 0.6*[1 1];
            elseif (mosaicRadialEccDegs <= 8)
                mappedSizeDegs = 1.0*[1 1];
            elseif (mosaicRadialEccDegs <= 12)
                mappedSizeDegs = 2.0*[1 1];
            elseif (mosaicRadialEccDegs <= 16)
                mappedSizeDegs = 3.0 *[1 1];
            else
                mappedSizeDegs = 3.0*[1 1]; % More than 4 with AO optics not feasible in my 192 GB RAM M3 MacStudio
            end
    else
            theRGCIndices = theMRGCMosaic.indicesOfRGCsWithinROI(mappedPositionDegs, 0.1*theMRGCMosaic.sizeDegs);
            
            if (isempty(theRGCIndices))
                error('No RGCs exist within the specified area (%f,%f) from xyPos: %2.2f,%2.2f\nfor mosaic at ecc: %f,%f and size: %f,%f\n', ...
                    0.1*theMRGCMosaic.sizeDegs(1), 0.1*theMRGCMosaic.sizeDegs(2), ...
                    mappedPositionDegs(1), mappedPositionDegs(2), ...
                    theMRGCMosaic.eccentricityDegs(1), theMRGCMosaic.eccentricityDegs(2), ...
                    theMRGCMosaic.sizeDegs(1), theMRGCMosaic.sizeDegs(2));
            end
            theTargetRGCindex = theRGCIndices(1);
            maxNeighborsNum = targetNumberOfMappedCells(1);
            nearbyRGCindices = theMRGCMosaic.neighboringRGCsToTargetRGCs(theTargetRGCindex, maxNeighborsNum);
            theXYpositionsOfTheTargetNumberOfMappedCells = theMRGCMosaic.rgcRFpositionsDegs(nearbyRGCindices,:);
            minXY = min(theXYpositionsOfTheTargetNumberOfMappedCells,[], 1);
            maxXY = max(theXYpositionsOfTheTargetNumberOfMappedCells,[], 1);
            mappedSizeDegs = max(maxXY-minXY)*[1 1];
    end
end