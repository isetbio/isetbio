function computeOuterSegmentLengthEccVariationAttenuationFactors(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('useParfor', true, @islogical);
    p.parse(varargin{:});
    useParfor = p.Results.useParfor;
    
    % Compute radial eccDegs of all cones
    coneEccentricityDegs = sqrt((obj.coneRFpositionsDegs(:,1)).^2 + (obj.coneRFpositionsDegs(:,2)).^2);
    
    % Compute outer segment length for each cone based on its eccentricity
    [osLengthMicrons, osLengthMicronsFoveal] = cMosaic.outerSegmentLengthFromEccentricity(coneEccentricityDegs);

    % Compute relative to the fovea os length
    relativeWithRespectToFoveaOuterSegmentLength = osLengthMicrons / osLengthMicronsFoveal;
    
    % Instantiate a cPhotoPigment object
    p = cPhotoPigment;
    
    % Retrieve the axial optical density of L-cones
    fovealLconeAxialOpticalDensity = p.opticalDensity(1);
    % Retrieve the absorbance spectrum for L-cones
    LconeAbsorbanceSpectrum = p.absorbance(:,1);
    % Compute the foveal absorptance
    fovealLconeAbsorptanceSpectrum = 1 - 10 .^ (-LconeAbsorbanceSpectrum * fovealLconeAxialOpticalDensity);
    
    % Compute the foveal  axial outer segment length assuming that
    % all cone types have same axial length
    % Specific density of foveal L-cones is 0.013 +/- 0.002 per micron for the L-cones (Bowmaker et al., 1978)
    LConeSpecificDensity = 0.013;
    fovealConeOuterSegmentAxialLength = fovealLconeAxialOpticalDensity / LConeSpecificDensity;
    
    % Estimate how axial outer segment lengths changes with the eccentricity of cones in the mosaic
    eccDependentOuterSegmentAxialLengths = fovealConeOuterSegmentAxialLength * relativeWithRespectToFoveaOuterSegmentLength;
    
    % Estimate optical densities of all cones, asuming specific density & absorbance do not change with eccentricity
    eccDependentAxialOpticalDensities = eccDependentOuterSegmentAxialLengths * LConeSpecificDensity;
    
    conesNum = size(obj.coneRFpositionsDegs,1);
    attenuationFactors = ones(1, conesNum);
    if (useParfor)
        parfor coneIndex = 1:conesNum
            eccDependentLconeAbsorptanceSpectrum = 1 - 10 .^ (-LconeAbsorbanceSpectrum * eccDependentAxialOpticalDensities(coneIndex));
            photonAbsorptionChange = eccDependentLconeAbsorptanceSpectrum ./ fovealLconeAbsorptanceSpectrum;
            % take the mean sensitivity change across all wavelengths
            attenuationFactors(coneIndex) = mean(photonAbsorptionChange);
        end
    else
        for coneIndex = 1:conesNum
            eccDependentLconeAbsorptanceSpectrum = 1 - 10 .^ (-LconeAbsorbanceSpectrum * eccDependentAxialOpticalDensities(coneIndex));
            photonAbsorptionChange = eccDependentLconeAbsorptanceSpectrum ./ fovealLconeAbsorptanceSpectrum;
            % take the mean sensitivity change across all wavelengths
            attenuationFactors(coneIndex) = mean(photonAbsorptionChange);
        end
    end
    
    obj.outerSegmentLengthEccVariationAttenuationFactors = attenuationFactors;
end