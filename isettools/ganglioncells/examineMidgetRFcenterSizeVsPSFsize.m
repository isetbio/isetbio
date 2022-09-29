function examineMidgetRFcenterSizeVsPSFsize()

    % Eccentricity in temporal retina
    radialEcc = 12;
    retinaQuadrant = 'temporal meridian';
    whichEye = 'right eye';
    [horizontalEcc, verticalEcc] = cMosaic.eccentricitiesForRetinaMeridianInEye(...
            radialEcc, retinaQuadrant, whichEye);

    eccDegs = [horizontalEcc verticalEcc];
    opticsDataBase = 'Artal2012';
    opticsSubjectRankOrder = 10;

    [theMidgetRFcenterRcDegs, thePSFRcDegs] = analyzeMidgetRFcenterAndPSF(eccDegs, ...
        opticsDataBase, opticsSubjectRankOrder, whichEye);
end

function [theMidgetRFcenterRcDegs, thePSFRcDegs] = analyzeMidgetRFcenterAndPSF(eccDegs, ...
        opticsDataBase, opticsSubjectRankOrder, whichEye)

    sizeDegs = 0.5+(sqrt(sum(eccDegs.^2,2))+1)*0.4*[0.5 0.25];
    

    theConeMosaic = generateConeMosaic(eccDegs, sizeDegs, whichEye);
    theConeMosaic.visualize();

    theMidgetRFcenterRcDegs = 1;
    thePSFRcDegs = 1;
end

  
function theConeMosaic = generateConeMosaic(eccDegs, sizeDegs, whichEye)
    % Set cone aperture modifiers
    % Use a Gaussian cone aperture with
    % sigma equal to 0.204 x inner segment diameter (cone diameter)
    sigmaGaussian = 0.204;  % From McMahon et al, 2000
    coneApertureModifiers = struct(...
            'smoothLocalVariations', true, ...
            'sigma',  sigmaGaussian, ...
            'shape', 'Gaussian');

    sourceLatticeSizeDegs = 60;
    customDegsToMMsConversionFunction = @(x)RGCmodels.Watson.convert.rhoDegsToMMs(x);
    customMMsToDegsConversionFunction = @(x)RGCmodels.Watson.convert.rhoMMsToDegs(x);

    theConeMosaic = cMosaic(...
       'sourceLatticeSizeDegs', sourceLatticeSizeDegs, ...
       'eccentricityDegs', eccDegs, ...
       'sizeDegs', sizeDegs, ...
       'whichEye', whichEye, ...
       'coneDensities', [0.6 0.3 0.1], ...
       'overlappingConeFractionForElimination', 0.5, ...
       'rodIntrusionAdjustedConeAperture', true, ...
       'coneApertureModifiers', coneApertureModifiers, ...
       'customDegsToMMsConversionFunction', customDegsToMMsConversionFunction, ...
       'customMMsToDegsConversionFunction', customMMsToDegsConversionFunction);
end



