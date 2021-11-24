function [correctionFactors, outerSegmentLengthAttenationFactors, innerSegmentDiameterBoostFactors] = computeConeEfficiencyCorrectionFactors(obj, triggerFunctionName, varargin)
% Static method for computing ecc-based absorption correction factors
%
% Syntax:
%   correctionFactors = COMPUTECONEEFFICIENCYCORRECTIONFACTORS(aConeMosaicHexObject, triggerFunctionName, varargin)
%
% Description:
%    This method is called by @coneMosaic's computeSingleFrame method
%    when the mosaic has hexagonal packing AND its eccBasedConeDensity flag
%    is set to true, AND its eccBasedConeQuantalEfficiency is set to true.
%    The purpose of this method is to compute correction factors for the
%    absorptions by accounting for each cone's variation in outer segment 
%    length and inner segment aperture variation with eccentricity.
%
% Inputs:
%    obj          - A coneMosaicHex object
%    triggerFunctionName - The function calling this method (debug purpose)
%
% Outputs:
%    correctionFactors  - The computed correction factors
%
% Optional key/value pairs:
%    'beVerbose'                - Boolean, Whether or not to display the
%                                 calling function name (used for debugging)
%                                 Default is False.
%    'plotEccDependentChanges'  - Boolean, Whether or not to plot how the 
%                                 aperture size and outer segment length 
%                                 vary with eccentricity.
%                                 Default is False.
%
% See Also:
%    coneMosaic.compute
%

% History:
%    06/14/18  NPC, ISETBIO Team    Wrote it

    p = inputParser;
    p.addRequired('obj', @(x)(isa(x,'coneMosaic')));
    p.addRequired('triggerFunctionName', @(x)(ischar(x)));
    p.addParameter('beVerbose', false, @islogical);
    p.addParameter('plotEccDependentChanges', false, @islogical);
    p.parse(obj, triggerFunctionName, varargin{:});
    
    if (p.Results.beVerbose)
        fprintf('>>> Computing ecc-based correction factors in cone quantal efficiency.\n\tTriggerred by ''%s'' method\n', triggerFunctionName);
    end

    if (p.Results.plotEccDependentChanges)
        figure(); clf;
        colors = [1 0 0; 0 1 0; 0 0 1];
    end
    
    correctionFactors = zeros(obj.rows, obj.cols);
    innerSegmentDiameterBoostFactors = correctionFactors;
    outerSegmentLengthAttenationFactors = correctionFactors;
    
    for coneTypeIndex = 1:3
        coneIndices = find(obj.pattern == coneTypeIndex+1);
        % Compute cone eccentricities in meters
        coneXYEccentricities = obj.coneLocs(coneIndices,:) / obj.resamplingFactor;
        coneEccentricitiesInMeters = (sqrt(sum(coneXYEccentricities.^2,2)))';
        coneEccentricitiesInDegs = coneEccentricitiesInMeters*1e6/obj.micronsPerDegree;
        
        % Compute cone angles
        coneAnglesInDegrees = atan2(squeeze(coneXYEccentricities(:,2)), squeeze(coneXYEccentricities(:,1))) / pi * 180;
        
        % Compute cone aperture for each cone based on its 2D location
        [~, apertureMeters, ~] = coneSizeReadData(...
            'eccentricity',coneEccentricitiesInMeters,...
            'angle',coneAnglesInDegrees);
        [~,apertureMetersAtZeroEcc, ~] = coneSizeReadData(...
            'eccentricity', 0.0,...
            'angle', 0.0);
        
        % Compute outer segment length for each cone based on its eccentricity
        [osLengthMicrons, osLengthAtZeroEcc] = outerSegmentLengthFromEccentricity(coneEccentricitiesInDegs);
        
        [correctionFactors(coneIndices), innerSegmentDiameterBoostFactors(coneIndices), outerSegmentLengthAttenationFactors(coneIndices)] = ...
            computeAbsorptionCorrectionFactors(...
              apertureMeters/apertureMetersAtZeroEcc, ...
              osLengthMicrons/osLengthAtZeroEcc, coneTypeIndex);
        
        if (p.Results.plotEccDependentChanges)
            visualizeOSlengthAndISaperture(coneEccentricitiesInDegs, ...
                apertureMeters*1e6, osLengthMicrons, squeeze(colors(coneTypeIndex,:)));
        end
    end  
end

function [correctionFactors, sensitivityChangeDueToInnerSegmentDiam, sensitivityChangeDueToOuterSegmentLength]= computeAbsorptionCorrectionFactors(...
    apertureChangeWithRespectToZeroEcc, osChangeWithRespectToZeroEcc, coneTypeIndex)
    
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
    
    correctionFactors = sensitivityChangeDueToInnerSegmentDiam .* sensitivityChangeDueToOuterSegmentLength;
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

function visualizeOSlengthAndISaperture(coneEccentricitiesInDegs, apertureMicrons, osLengthMicrons, color)
    
    subplot(1,2,1)
    plot(coneEccentricitiesInDegs, apertureMicrons, '.', 'Color', color);
    hold on
    ylabel('inner segment diam (microns)');
    xlabel('eccentricity (degs)');
    
    subplot(1,2,2)
    plot(coneEccentricitiesInDegs, osLengthMicrons, '.', 'Color', color);
    hold on;
    ylabel('outer segment length (microns)');
    xlabel('eccentricity (degs)');
end

