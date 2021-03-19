% Since this correction is applied on the optical image we need to compute
% correction factors at each emPosition because the relative position
% between the optical image and the mosaic varies during the emPath
function macularPigmentDensityBoostFactors = computeMPBoostFactors(obj, oiPositionsDegs, emPositionDegs, oiWave, oiSize, oiResMicrons)

    if (cachedMacularPigmentDensityBoostFactorsIsValid(obj, oiWave, oiSize, oiResMicrons, emPositionDegs))
        fprintf('No change in oiWave, oiSize, oiResMicrons, and same emPosition, so will use cached macularPigmentDensityBoostFactors.\n');
        macularPigmentDensityBoostFactors = obj.cachedMacularPigmentDensityBoostFactors.macularPigmentDensityBoostFactors; 
        return;
    end
    
    fprintf('Computing new MP density boost factors\n');
    % Adjust oiPositionsDegs to take into account the current eye position
    oiPositionsDegs = bsxfun(@minus, oiPositionsDegs, emPositionDegs);
    
    % Compute ecc-based MP optical densities
    eccBasedMacularPigmentDensities = obj.macular.eccDensity([], 'eccDegs2', sum(oiPositionsDegs.^2,2));

    % And corresponding transmittances. For a 2x2 deg oi, this step takes 2.7 seconds
    eccBasedMacularPigmentTransmittances = 10.^(-eccBasedMacularPigmentDensities * obj.macular.unitDensity');
    
    % Boost factors. For a 2x2 deg oi, this step takes 1.0 seconds
    macularPigmentDensityBoostFactors = eccBasedMacularPigmentTransmittances ./ repmat(obj.macular.transmittance', [size(eccBasedMacularPigmentTransmittances,1) 1]);

    % Save to cache
    obj.cachedMacularPigmentDensityBoostFactors.eccVaryingMacularPigmentDensity = obj.eccVaryingMacularPigmentDensity;
    obj.cachedMacularPigmentDensityBoostFactors.macularPigmentDensityBoostFactors = macularPigmentDensityBoostFactors;
    obj.cachedMacularPigmentDensityBoostFactors.oiResMicrons = oiResMicrons;
    obj.cachedMacularPigmentDensityBoostFactors.oiSize = oiSize;
    obj.cachedMacularPigmentDensityBoostFactors.oiWave = oiWave;
    obj.cachedMacularPigmentDensityBoostFactors.emPositionDegs = emPositionDegs;
end


function isValid =  cachedMacularPigmentDensityBoostFactorsIsValid(obj, oiWave, oiSize, oiResMicrons, emPositionDegs)
    if (isempty(obj.cachedMacularPigmentDensityBoostFactors))
        isValid = false;
    else
        isValid = ...
            (oiResMicrons == obj.cachedMacularPigmentDensityBoostFactors.oiResMicrons) && ...
            (oiSize(1) == obj.cachedMacularPigmentDensityBoostFactors.oiSize(1)) && ...
            (oiSize(2) == obj.cachedMacularPigmentDensityBoostFactors.oiSize(2)) && ...
            (numel(oiWave) == numel(obj.cachedMacularPigmentDensityBoostFactors.oiWave)) && ...
            all(oiWave == obj.cachedMacularPigmentDensityBoostFactors.oiWave) && ...
            (emPositionDegs(1) == obj.cachedMacularPigmentDensityBoostFactors.emPositionDegs(1)) && ...
            (emPositionDegs(2) == obj.cachedMacularPigmentDensityBoostFactors.emPositionDegs(2)) && ...
            (obj.eccVaryingMacularPigmentDensity == obj.cachedMacularPigmentDensityBoostFactors.eccVaryingMacularPigmentDensity);
    end
end
