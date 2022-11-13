function generateMidgetRGCMosaics()
    
    eccDegs = [0 0];
    sizeDegs = [0.5 0.5];


    theMidgetRGCmosaic = midgetRGCMosaic(...
            'sourceLatticeSizeDegs', 60, ...
            'eccentricityDegs', eccDegs, ...
            'sizeDegs', sizeDegs ...
            );
    theMidgetRGCmosaic.visualizeRetinalRFs();



    save(sprintf('mRGCmosaic_eccDegs_%2.2f.mat', eccDegs(1)), ...
        'theMidgetRGCmosaic', '-v7.3');


end

function generate