function [coneApertureProfile, theConeApertureSupportDegs] = generateConeApertureProfileForDeconvolution(thePSFsupportDegs, coneAperturesDegs)
    
    % Spatial support for the cone aperture
    spatialSampleSize = thePSFsupportDegs(2)-thePSFsupportDegs(1);
    theConeApertureSupportDegs = (spatialSampleSize:spatialSampleSize:1.0*max(coneAperturesDegs));
    theConeApertureSupportDegs = [-fliplr(theConeApertureSupportDegs) 0 theConeApertureSupportDegs];
    [Xcone, Ycone] = meshgrid(theConeApertureSupportDegs, theConeApertureSupportDegs);
    
    % Cone aperture profile
    apertureShape = 'disk';
    coneRadiusDegs = 0.5*mean(coneAperturesDegs);
    
    if (strcmp(apertureShape, 'disk'))
        r = sqrt(Xcone.^2+Ycone.^2);
        coneApertureProfile = zeros(size(Xcone));
        coneApertureProfile(r<=coneRadiusDegs) = 1;
        
    elseif (strcmp(apertureShape, 'flattopGaussian'))
        nExp = 5;
        coneApertureProfile = (exp(-(abs(Xcone/coneRadiusDegs)).^nExp) .* ...
                               exp(-(abs(Ycone/coneRadiusDegs)).^nExp));
                           
    elseif (strcmp(apertureShape, 'Gaussian'))
        nExp = 2;
        coneRadiusSigma = coneRadiusDegs/3;
        coneApertureProfile = (exp(-(abs(Xcone/coneRadiusSigma)).^nExp) .* ...
                               exp(-(abs(Ycone/coneRadiusSigma)).^nExp));
    end
    
end