% Method for suggesting an optimal pixel size for the input scene
% depending on the aperture size of the mosaic's cones and whether
% the computation is done using ecc-dependent blur mode
function minScenePixelSizeDegs = suggestedScenePixelSizeDegs(obj, eccVaryingConeBlur)

    % Compute cone aperture diameters based on their local spacing
    coneApertureDiametersDegs = obj.coneRFspacingsDegs * obj.coneApertureToDiameterRatio;
    
    if (eccVaryingConeBlur)
        minConeApertureDiameterDegs = min(coneApertureDiametersDegs);
    else
        minConeApertureDiameterDegs = median(coneApertureDiametersDegs);
    end
    
    % We want at least a 3x3 aperture kernel
    minScenePixelSizeDegs = minConeApertureDiameterDegs/3.0;
end
