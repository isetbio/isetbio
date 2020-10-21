function sizeMicrons = sizeVisualDegsToSizeRetinalMicrons(sizeVisualDegs, eccDegs)
% Convert a retinal size in visual degrees at given ecc (also in degs) to its
% equivalent linear size in retinal microns
    
    leftCoordDegs  = (eccDegs - sizeVisualDegs/2);
    rightCoordDegs = (eccDegs + sizeVisualDegs/2);
    sizeMicrons =  1e3 * ( RGCmodels.Watson.convert.rhoDegsToMMs(rightCoordDegs) ...
                          -RGCmodels.Watson.convert.rhoDegsToMMs(leftCoordDegs));
          
end

