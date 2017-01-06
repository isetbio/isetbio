function t_coneMosaicHex6

    %% Initialize
    ieInit; clear; close all;

    % Freeze random number generator
    rng('default'); rng(219347);

    customLambda = 1.7; 
    resamplingFactor = 9;
    spatiallyVaryingConeDensity = false;
    
    for rotationDegs = 0:0.25:30
        theHexMosaic = coneMosaicHex(...
            resamplingFactor, spatiallyVaryingConeDensity, customLambda, ...
            'rotationDegs', rotationDegs, ...
            'size', [16 16]);
        theHexMosaic.visualizeGrid('visualizedConeAperture', 'lightCollectingArea');
    end
    
end

