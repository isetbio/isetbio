function [contourData, RF] = subregionOutlineContourFromPooledCones(...
    conePos, coneRc, poolingWeights, ...
    xSupport, ySupport, spatialSupportSamples)

    % Compute spatial support
    xSep = max(coneRc)*2*sqrt(numel(poolingWeights));
    if (isempty(xSupport))
        xx = conePos(:,1);
        xSupport = linspace(min(xx)-xSep,max(xx)+xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = conePos(:,2);
        ySupport = linspace(min(yy)-xSep,max(yy)+xSep,spatialSupportSamples);
    end

    [X,Y] = meshgrid(xSupport, ySupport);
    spatialSupportXY(:,1) = xSupport(:);
    spatialSupportXY(:,2) = ySupport(:);

    RF = zeros(size(X));
    for iCone = 1:numel(poolingWeights)
        % Characteristic radius of the input RF
        rC = coneRc(iCone);
        % Compute aperture2D x weight
        XX = X-conePos(iCone,1);
        YY = Y-conePos(iCone,2);
        theAperture2D = poolingWeights(iCone) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
        % Accumulate 2D apertures
        RF = RF + theAperture2D;
    end

    zLevels(1) = 0.02*min(poolingWeights);
    zLevels(2) = max(poolingWeights);
    s = mRGCMosaic.contourDataFromDensityMap(spatialSupportXY, RF, zLevels);
    contourData = s{1};
end
