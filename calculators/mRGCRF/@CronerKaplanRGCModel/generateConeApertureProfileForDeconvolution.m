function coneApertureProfile = generateConeApertureProfileForDeconvolution(thePSFsupportDegs, coneAperturesDegs)
    
    % Spatial support for the cone aperture
    spatialSampleSize = thePSFsupportDegs(2)-thePSFsupportDegs(1);
    theConeApertureSupportDegs = (spatialSampleSize:spatialSampleSize:0.6*max(coneAperturesDegs));
    theConeApertureSupportDegs = [-fliplr(theConeApertureSupportDegs) 0 theConeApertureSupportDegs];
    [Xcone, Ycone] = meshgrid(theConeApertureSupportDegs, theConeApertureSupportDegs);
    
    % Cone aperture profile
    apertureShape = 'flattopGaussian';
    coneCharacteristicRadiusDegs = 0.5*mean(coneAperturesDegs)/3;
    if (strcmp(apertureShape, 'disk'))
        r = sqrt(Xcone.^2+Ycone.^2);
        coneApertureProfile = zeros(size(Xcone));
        coneApertureProfile(r<=0.5*mean(coneAperturesDegs)) = 1;
    elseif (strcmp(apertureShape, 'flattopGaussian'))
        coneApertureProfile = (exp(-(Xcone/coneCharacteristicRadiusDegs).^2) .* exp(-(Ycone/coneCharacteristicRadiusDegs).^2));
        flatTopThreshold = exp(-5)*exp(-5);
        coneApertureProfile(coneApertureProfile<flatTopThreshold) = 0;
        coneApertureProfile = coneApertureProfile .^ 0.01;
    else  
        coneApertureProfile = exp(-(Xcone/coneCharacteristicRadiusDegs).^2) .* exp(-(Ycone/coneCharacteristicRadiusDegs).^2);
    end
end