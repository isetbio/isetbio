function testMRGCmosaic
    % Instantiate a [0.3 x 0.3] deg wide midget RGC mosaic positioned 
    % at an eccentricity of [-1 0] degs , in the right eye
    theMidgetRGCmosaic = mRGCmosaic([-1 0], [0.3 0.3], 'right');
    
    % Visualize mRGC positions together with imported ecc-varying cone positions
    % and together with cone positions in the equivalent employed reg-hex mosaic
    theMidgetRGCmosaic.visualizeInputPositions();
    
    % Visualize connections of cones to mRGC RF centers
    theMidgetRGCmosaic.visualizeConeMosaicTesselation('degrees');
    
    theMidgetRGCmosaic.visualizeSynthesizedParams();
    
    theMidgetRGCmosaic.visualizeConeWeights();
end

