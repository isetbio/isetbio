% Since this correction is applied on the optical image we need to compute
% correction factors at each emPosition because the relative position
% between the optical image and the mosaic varies during the emPath
function macularPigmentDensityBoostFactors = computeMPBoostFactors(obj, oiPositionsDegs, emPositionDegs, oiWave, oiSize, oiResMicrons)

    if (cachedMacularPigmentDensityBoostFactorsIsValid(obj, oiWave, oiSize, oiResMicrons, emPositionDegs))
        macularPigmentDensityBoostFactors = obj.cachedMacularPigmentDensityBoostFactors.macularPigmentDensityBoostFactors; 
        return;
    end
    
    % Adjust oiPositionsDegs to take into account the current eye position
    oiPositionsDegs = bsxfun(@minus, oiPositionsDegs, emPositionDegs);

    % Compute the macular pigment density boost factors for all oi pixel positions
    macularPigmentDensityBoostFactors = cMosaic.macularPigmentBoostFactors(obj.macular, oiPositionsDegs);

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
        unchangedOIRes  = (oiResMicrons == obj.cachedMacularPigmentDensityBoostFactors.oiResMicrons);
        unchangedOISize = all(oiSize == obj.cachedMacularPigmentDensityBoostFactors.oiSize);
        unchangedOIWave = (numel(oiWave) == numel(obj.cachedMacularPigmentDensityBoostFactors.oiWave)) && ...
                           all(oiWave == obj.cachedMacularPigmentDensityBoostFactors.oiWave);
        unchangedEMPos = all(emPositionDegs == obj.cachedMacularPigmentDensityBoostFactors.emPositionDegs);
        unchangedMPDensityFlag = (obj.eccVaryingMacularPigmentDensity == obj.cachedMacularPigmentDensityBoostFactors.eccVaryingMacularPigmentDensity);
        isValid = unchangedOIRes && unchangedOISize && unchangedOIWave && unchangedEMPos && unchangedMPDensityFlag;
    
        beVerbose = false;
        if (beVerbose)
            if (isValid)
                fprintf('No change in oiWave, oiSize, oiResMicrons, and same emPosition, so will use cached macularPigmentDensityBoostFactors.\n');
            else
                fprintf('Recomputing macularPigmentDensityBoostFactors ');
                if (~unchangedOIRes)
                    fprintf('Because of a difference in the oiRes\n');
                end
                if (~unchangedOISize)
                    fprintf('Because of a difference in the oiSize\n');
                end
                if (~unchangedOIWave)
                    fprintf('Because of a difference in the oiWave\n');
                end
                if (~unchangedEMPos)
                    fprintf('Because of a difference in the eye position\n');
                end
                if (~unchangedMPDensityFlag)
                    fprintf('Because of a change in the eccVaryingMacularPigmentDensity Flag\n');
                end
            end
        end
    end
end
