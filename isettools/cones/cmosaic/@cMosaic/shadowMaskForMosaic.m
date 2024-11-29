function shadowMaskForMosaic(obj, shadowMaskMatrix)

    if (isempty(shadowMaskMatrix))
        obj.shadowMask = [];
        return;
    end

    interpolationMethod = 'nearest';
    extrapolationMethod = 'nearest';

    % Compute gridded interpolant
    F = scatteredInterpolant(shadowMaskMatrix(:,1),shadowMaskMatrix(:,2),shadowMaskMatrix(:,3),interpolationMethod, extrapolationMethod);

    % Interpolate at cone positions
    shadowMask = F([obj.coneRFpositionsDegs(:,1) obj.coneRFpositionsDegs(:,2)]);

    obj.shadowMask = reshape(shadowMask, [1 obj.conesNum]);
end