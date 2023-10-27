function [osLengthMicrons, osLengthMicronsFoveal] = outerSegmentLengthFromEccentricity(eccDegs)
% Scanned data (eccentricity, osLength) from Figure 1 (right panel)
% Banks, Sekuler and Anderson (1991). "Peripheral spatial vision: limits
% imposed by optics, photoreceptors and receptor pooling".
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
  osLengthMicronsFoveal = s(1,2);
end