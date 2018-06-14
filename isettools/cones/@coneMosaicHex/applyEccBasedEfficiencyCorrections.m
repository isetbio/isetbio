function absorptionsCorrected = applyEccBasedEfficiencyCorrections(obj, absorptions)
% Static method for computing absorptions with ecc-based efficiency corrections
%
% Syntax:
%   absorptionsCorrected  = APPLYECCBASEDEFFICIENCYCORRECTIONS(aConeMosaicHexObject, absorptions)
%
% Description:
%    This method is called by @coneMosaic's computeSingleFrame method
%    when the mosaic has hexagonal packing AND its eccBasedConeDensity flag
%    is set to true, i.e, when the mosaic has an eccentricity-based cone
%    density. The purpose of this method is to compute corrected absorptions
%    by accounting for each cone's variation in outer segment length and 
%    inner segment aperture variation with eccentricity.
%
% Inputs:
%    obj          - A coneMosaicHex object
%    absorptions  - The absorptions integrated over default cone aperture (zero eccentricity)
%
% Outputs:
%    absorptionsCorrected  - The eccentricity based corrected cone absorptions
%
% Optional key/value pairs:
%    None' 
%
% See Also:
%    coneMosaic.computeSingleFrame
%

% History:
%    06/14/18  NPC, ISETBIO Team    Wrote it

    fprintf('Applying ecc-based correction in absorptions\n');
    [rows,cols,coneTypes] = size(absorptions);
    
    % Compute cone eccentricities in meters
    coneXYEccentricities = obj.coneLocs / obj.resamplingFactor;
    coneEccentricitiesInMeters = (sqrt(sum(coneXYEccentricities.^2,2)))';
    coneEccentricitiesInDegs = coneEccentricitiesInMeters*1e6/obj.micronsPerDegree;
    
    % Compute cone angles
    coneAnglesInDegrees = atan2(squeeze(coneXYEccentricities(:,2)), squeeze(coneXYEccentricities(:,1))) / pi * 180;
    absorptionsCorrected  = absorptions*0;
    
    % Compute cone aperture for each cone based on its 2D location
    [~, apertureMeters, ~] = coneSizeReadData(...
        'eccentricity',coneEccentricitiesInMeters,...
        'angle',coneAnglesInDegrees);
    [~,apertureMetersAtZeroEcc, ~] = coneSizeReadData(...
        'eccentricity', 0.0,...
        'angle', 0.0);
    
    % Compute outer segment length for each cone based on its eccentricity
    [osLengthMicrons, osLengthAtZeroEcc] = outerSegmentLengthFromEccentricity(coneEccentricitiesInDegs);

    plotEccDependentChanges = false;
    if (plotEccDependentChanges)
        visualizeOSlengthAndISaperture(coneEccentricitiesInDegs, ...
            apertureMeters*1e6, osLengthMicrons);
    end

    % Apply corrections separately for each cone type
    for coneTypeIndex = 1:coneTypes
        absorptionsForTargetCone = reshape(...
            squeeze(absorptions(:,:,coneTypeIndex)), ...
            [1, rows*cols]);

        absorptionsCorrected(:,:, coneTypeIndex) = ...
            reshape(computeAbsorptionsUsingEccentricity(...
                            absorptionsForTargetCone, ...
                            apertureMeters/apertureMetersAtZeroEcc, ...
                            osLengthMicrons/osLengthAtZeroEcc, ...
                            coneTypeIndex), ...
                    [rows cols 1]);
    
    end
end

function absorptionsCorrected = computeAbsorptionsUsingEccentricity(absorptions, apertureChangeWithRespectToZeroEcc, osChangeWithRespectToZeroEcc, coneTypeIndex)
    
    % Sensitivity change due to change in aperture diameter (A = r^2);
    sensitivityChangeDueToInnerSegmentDiam = apertureChangeWithRespectToZeroEcc.^2;
    
    % Sensitivity change due to change in outer segment length
    sensitivityChangeDueToOuterSegmentLength = zeros(size(osChangeWithRespectToZeroEcc));
    
    % Get the axial optical density for the target cone type
    p = photoPigment;
    axialOpticalDensity = p.opticalDensity(coneTypeIndex);
    
    % Get the absorbance & absorptance for the target cone type at zero eccentricity
    fovealAbsorbance = squeeze(p.absorbance(:,coneTypeIndex));
    fovealAbsorptance = 1 - 10 .^ (-fovealAbsorbance * axialOpticalDensity);
    
    % Compute the absorbance & absorptance for the target cone type for
    % cone at all eccentricities
    for coneIndex = 1:numel(osChangeWithRespectToZeroEcc)
        peripheralAbsorbance = fovealAbsorbance * osChangeWithRespectToZeroEcc(coneIndex);
        peripheralAbsorptance = 1 - 10 .^ (-peripheralAbsorbance * axialOpticalDensity);
        photonAbsorptionChange = peripheralAbsorptance./fovealAbsorptance;
        % take the mean sensitivity change across all wavelengths
        sensitivityChangeDueToOuterSegmentLength(coneIndex) = mean(photonAbsorptionChange);
    end
    
    absorptionsCorrected = absorptions .* sensitivityChangeDueToInnerSegmentDiam .* sensitivityChangeDueToOuterSegmentLength;
end

function [osLengthMicrons, osLengthAtZeroEcc] = outerSegmentLengthFromEccentricity(eccDegs)
    % scanned data (eccentricity, osLength) from Figure 1 (right panel)
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
end

function visualizeOSlengthAndISaperture(coneEccentricitiesInDegs, apertureMicrons, osLengthMicrons)
    figure(); clf;
    subplot(1,2,1)
    plot(coneEccentricitiesInDegs, apertureMicrons, 'k.');
    ylabel('inner segment diam (microns)');
    xlabel('eccentricity (degs)');
    
    subplot(1,2,2)
    plot(coneEccentricitiesInDegs, osLengthMicrons, 'k.');
    ylabel('outer segment length (microns)');
    xlabel('eccentricity (degs)');
end

