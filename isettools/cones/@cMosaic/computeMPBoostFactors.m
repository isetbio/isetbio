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
    
    % Compute ecc-based MP optical densities
    eccBasedMacularPigmentDensities = obj.macular.eccDensity([], 'eccDegs2', sum(oiPositionsDegs.^2,2));

    
    computeMethod = 2;
    
    if (computeMethod == 1)
        %Total: 59.1
        eccBasedMacularPigmentTransmittances = 10.^(-eccBasedMacularPigmentDensities * obj.macular.unitDensity'); % 43.5
        macularPigmentDensityBoostFactors = eccBasedMacularPigmentTransmittances ./ repmat(obj.macular.transmittance', [size(eccBasedMacularPigmentTransmittances,1) 1]); % 11.1
    elseif (computeMethod == 2)
        % Total: 57.1
        eccBasedMacularPigmentTransmittances = 10.^(-eccBasedMacularPigmentDensities * obj.macular.unitDensity'); % 43.4
        macularPigmentDensityBoostFactors = bsxfun(@rdivide, eccBasedMacularPigmentTransmittances, obj.macular.transmittance'); %9.3
    else
        % Total: 62
        macularUnitDensity = obj.macular.unitDensity';
        a = eccBasedMacularPigmentDensities * macularUnitDensity; % 8
        b = repmat(obj.macular.density * macularUnitDensity,  [size(eccBasedMacularPigmentDensities,1) 1]); %6
        diff = b-a; % 8
        macularPigmentDensityBoostFactors = 10 .^ diff; % 30
    end
    
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
