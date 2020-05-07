function [coneDensityPerMM2, maxDensity] = coneDensityFunctionFull(rfPositions, whichEye)
    eccentricitiesInMicrons = sqrt(sum(rfPositions .^ 2, 2));
    eccentricitiesInMeters = eccentricitiesInMicrons * 1e-6;
    angles = atan2(rfPositions(:, 2), rfPositions(:, 1)) / pi * 180;
    [~,~,coneDensityPerMM2] = coneSizeReadData('eccentricity', eccentricitiesInMeters, ...
        'angle', angles, 'whichEye', whichEye);  
    [~,~,maxDensity] = coneSizeReadData('eccentricity', 0, 'angle', 0, 'whichEye', whichEye); 
end