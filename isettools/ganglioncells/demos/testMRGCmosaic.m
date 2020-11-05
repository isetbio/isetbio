function testMRGCmosaic
    % Instantiate a [0.4 x 0.4] deg wide midget RGC mosaic positioned 
    % at an eccentricity of [-2 0] degs 
    theMidgetRGCmosaic = mRGCmosaic([-2 0], [0.2 0.2]);
    
    % Visualize the mosaic.
    theMidgetRGCmosaic.visualizeConeMosaicTesselation('degrees');
    
    theMidgetRGCmosaic.visualizeSynthesizedParams();
    
    theMidgetRGCmosaic.visualizeConeWeights();
end

