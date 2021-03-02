function computeOuterSegmentLengthEccVariationAttenuationFactors(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('useParfor', true, @islogical);
    p.parse(varargin{:});
    useParfor = p.Results.useParfor;
    
    % Compute radial eccDegs of all cones
    coneEccentricityDegs = sqrt((obj.coneRFpositionsDegs(:,1)).^2 + (obj.coneRFpositionsDegs(:,2)).^2);
    
    % Compute outer segment length for each cone based on its eccentricity
    relativeWithRespectToFoveaOuterSegmentLength = outerSegmentLengthFromEccentricity(coneEccentricityDegs);
    
    
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

function relativeWithRespectToFoveaOuterSegmentLength = outerSegmentLengthFromEccentricity(eccDegs)
    % scanned data (eccentricity, osLength in microns) from Figure 1 (right panel)
    % Banks, Sekuler and Anderson (1991). Peripher spatial vision: limits
    % imposed by optics, photoreceptors and receptor pooling
    s = [ ...
        0.00  47.81;
        1.82  26.16;
        4.86  21.2;
        9.86  21.20;
        19.78  21.2;
        39.90  13.22;
    ];
  scannedData.eccDegsRaw = s(:,1);
  scannedData.lengthMicronsRaw = s(:,2);
  interpolationMethod = 'pchip';
  osLengthMicrons = interp1(scannedData.eccDegsRaw, scannedData.lengthMicronsRaw, eccDegs, interpolationMethod);
  osLengthAtZeroEcc = s(1,2);
  relativeWithRespectToFoveaOuterSegmentLength = osLengthMicrons / osLengthAtZeroEcc;
end

