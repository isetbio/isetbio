function macularPigmentDensityBoostFactors = macularPigmentBoostFactors(theMacular, positionsDegs)

    % Compute ecc-based MP optical densities
    eccBasedMacularPigmentDensities = theMacular.eccDensity([], 'eccDegs2', sum(positionsDegs.^2,2));

    computeMethod = 2;
    
    if (computeMethod == 1)
        %Total: 59.1
        eccBasedMacularPigmentTransmittances = 10.^(-eccBasedMacularPigmentDensities * theMacular.unitDensity'); % 43.5
        macularPigmentDensityBoostFactors = eccBasedMacularPigmentTransmittances ./ repmat(theMacular.transmittance', [size(eccBasedMacularPigmentTransmittances,1) 1]); % 11.1
    elseif (computeMethod == 2)
        % Total: 57.1
        eccBasedMacularPigmentTransmittances = 10.^(-eccBasedMacularPigmentDensities * theMacular.unitDensity'); % 43.4
        macularPigmentDensityBoostFactors = bsxfun(@rdivide, eccBasedMacularPigmentTransmittances, theMacular.transmittance'); %9.3
    else
        % Total: 62
        macularUnitDensity = theMacular.unitDensity';
        a = eccBasedMacularPigmentDensities * macularUnitDensity; % 8
        b = repmat(theMacular.density * macularUnitDensity,  [size(eccBasedMacularPigmentDensities,1) 1]); %6
        diff = b-a; % 8
        macularPigmentDensityBoostFactors = 10 .^ diff; % 30
    end

end