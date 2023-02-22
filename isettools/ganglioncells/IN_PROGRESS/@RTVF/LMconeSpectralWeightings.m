function [L,M, wavelengthSupport] = LMconeSpectralWeightings(theConeMosaic, theTargetEccDegs)

    if (theConeMosaic.eccVaryingMacularPigmentDensity)
        % Each cone has a different boost factor, based on its eccentricity
        macularPigmentBoostFactors = cMosaic.macularPigmentBoostFactors(theConeMosaic.macular, theTargetEccDegs);
    else
        % Single factor based on the eccentricity of the mosaic
        macularPigmentBoostFactors = cMosaic.macularPigmentBoostFactors(theConeMosaic.macular, theConeMosaic.eccentricityDegs);
    end
    boostVector = macularPigmentBoostFactors';
    
    % Quantal efficiencies with the foveal MP density and no lens
    coneQuantalEfficienciesFovealMP = theConeMosaic.qe;

    % Quantal efficiencies with the MP density at the target ecc and no lens
    coneQuantalEfficienciesEccBasedMP = diag(boostVector) * coneQuantalEfficienciesFovealMP;

    % Quantal efficiencies with the MP density at the target ecc and the lens
    theLens = Lens('wave',theConeMosaic.wave);
    lensTransmittance = theLens.transmittance;
    coneQuantalEfficienciesFinal = diag(lensTransmittance)*coneQuantalEfficienciesEccBasedMP;

    % Return the L- and M-cone quantal efficiencies
    L = coneQuantalEfficienciesFinal(:,1);
    M = coneQuantalEfficienciesFinal(:,2);
    wavelengthSupport = theConeMosaic.wave;
end