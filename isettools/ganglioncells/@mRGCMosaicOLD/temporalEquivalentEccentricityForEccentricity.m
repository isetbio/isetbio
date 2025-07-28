function temporalEquivalentEccDegs = temporalEquivalentEccentricityForEccentricity(obj, eccDegs)

    if (size(eccDegs,2) ~= 2)
        error('eccDegs must be an N x 2 matrix of (x,y) eccentricities.');
    end

    % Find which RGCs are in the nasal reridian
    horizontalEcc = eccDegs(:,1);
    switch (obj.whichEye)
        case 'left eye'
            nasalMeridianRGCindices = find(horizontalEcc < 0);

        case 'right eye'
            nasalMeridianRGCindices = find(horizontalEcc > 0);
    end % switch

    switch (obj.temporalEquivantEccentricityFactor)
        case 'WatanabeRodieckBased'
          % Multiply by 0.61 - See Watanabe & Rodieck 1989
          scalingFactor = 0.61;
        case 'ISETBioMosaicsBased'
           % Multiply by 0.70 - to match our asymmetry
           scalingFactor = 0.70;
    end % switch

    % Apply scaling factor to RGCs in the nasal meridian
    horizontalEcc(nasalMeridianRGCindices) = ...
        scalingFactor * horizontalEcc(nasalMeridianRGCindices);

    % Put it back to eccDegs
    temporalEquivalentEccDegs = eccDegs;
    temporalEquivalentEccDegs(:,1) = horizontalEcc;
 
    debugConversion = false;
    if (debugConversion)
        figure();
        subplot(1,2,1)
        plot(eccDegs(:,1), temporalEquivalentEccDegs(:,1), 'r.');
        axis 'square'
        set(gca, 'XLim', [-1 1], 'YLim', [-1 1])
        subplot(1,2,2)
        plot(eccDegs(:,2), temporalEquivalentEccDegs(:,2), 'b.');
        axis 'square'
        set(gca, 'XLim', [-1 1], 'YLim', [-1 1]);
    end

end