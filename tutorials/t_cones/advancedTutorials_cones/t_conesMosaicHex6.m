function t_conesMosaicHex6
% How to generate a hexagonal mosaic with a eccentricity based cone density
%
% Syntax:
%   t_conesMosaicHex6
% Description:
%    This tutorial covers how to generate hexagonal mosaic with
%    cone density based on eccentricity
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: Change plot style to incorporate subplots to display changes
%      in the cone mosaic.
%

    %% Initialize
    ieInit;
    clear;
    close all;

    % Freeze random number generator
    rng('default');
    rng(219347);

    customLambda = 1.7;
    resamplingFactor = 7;
    spatiallyVaryingConeDensity = false;
    
    for rotationDegs = 0:2:30
    % Mosaic Parameters
    %    resamplingFactor    - Numeric. Controls hex mosaic grid accuracy.
    %                          Factor of 7 listed here.
    %    eccBasedConeDensity - Boolean. Whether to have an eccentricity
    %                          based, spatially-varying density.
    %    rotationDegs          - Numeric. The degrees of rotation.
    %    size                  - Vector. The mosaic size.
        theHexMosaic = coneMosaicHex(resamplingFactor, ...
            'eccBasedConeDensity', spatiallyVaryingConeDensity, ...
            'customLambda', customLambda, ...
            'rotationDegs', rotationDegs, ...
            'size', [20 20]);
        rng(1);
        theHexMosaic.reassignConeIdentities();
        theHexMosaic.visualizeGrid(...
            'visualizedConeAperture', 'lightCollectingArea');
    end
    
end
