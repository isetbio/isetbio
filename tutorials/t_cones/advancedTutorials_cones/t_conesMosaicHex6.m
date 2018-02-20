function t_conesMosaicHex6
%%t_conesMosaicHex6

    %% Initialize
    ieInit; clear; close all;

    % Freeze random number generator
    rng('default'); rng(219347);

    customLambda = 1.7; 
    resamplingFactor = 7;
    spatiallyVaryingConeDensity = false;
    
    for rotationDegs = 0:2:30
        theHexMosaic = coneMosaicHex(resamplingFactor, ...
            'eccBasedConeDensity', spatiallyVaryingConeDensity, ...
            'customLambda', customLambda, ...
            'rotationDegs', rotationDegs, ...
            'size', [20 20] ...
            );
        rng(1);
        theHexMosaic.reassignConeIdentities();
        theHexMosaic.visualizeGrid('visualizedConeAperture', 'lightCollectingArea');
    end
    
end

