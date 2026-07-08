function contourData = subregionOutlineContourFromSpacing(...
    rgcPos, rgcRFradius,  ...
    xSupport, ySupport, spatialSupportSamples)

    % Compute spatial support
    xSep = rgcRFradius;
    if (isnan(xSep))
        contourData = [];
        return;
    end

    if (isempty(xSupport))
        xx = rgcPos(:,1);
        xSupport = linspace(min(xx)-xSep,max(xx)+xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = rgcPos(:,2);
        ySupport = linspace(min(yy)-xSep,max(yy)+xSep,spatialSupportSamples);
    end

    [X,Y] = meshgrid(xSupport, ySupport);
    spatialSupportXY(:,1) = xSupport(:);
    spatialSupportXY(:,2) = ySupport(:);
    
    % 2*sigma = radius
    rSigma = rgcRFradius/3;
    % Compute aperture2D x weight
    XX = X-rgcPos(1);
    YY = Y-rgcPos(2);
    theAperture2D = exp(-0.5*(XX/rSigma).^2) .* exp(-0.5*(YY/rSigma).^2);
    RF = theAperture2D;

    zLevels(1) = exp(-2.5);
    zLevels(2) = exp(-2.4);
    s = mRGCMosaic.contourDataFromDensityMap(spatialSupportXY, RF, zLevels);
    contourData = s{1};
end