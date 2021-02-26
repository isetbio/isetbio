function apertureKernel = generateApertureKernel(coneApertureDiameterMicrons, oiResMicrons)
    % Compute convolution kernel size in pixels
    apertureSamples = ceil(coneApertureDiameterMicrons / oiResMicrons);
    
    % Make sure it is odd
    if (mod(apertureSamples, 2) == 0)
        apertureSamples = apertureSamples + 1;
    end
    
    %fprintf('Aperture kernel is %d x %d pixels\n', apertureSamples, apertureSamples);
    %fprintf('Aperture diameter: %2.2f microns, optical image resolution: %2.2f microns\n', coneApertureDiameterMicrons , oiResMicrons);
    
    % Generate aperture kernel: pillbox with unit volume
    apertureKernel = zeros(apertureSamples, apertureSamples);
    [xx,yy] = sample2space(1:apertureSamples, 1:apertureSamples, oiResMicrons, oiResMicrons);
    [apertureRows, apertureCols] = meshgrid(xx,yy);
    apertureRadii = sqrt(apertureRows .^ 2 + apertureCols .^ 2);
    apertureKernel(apertureRadii <= 0.5*coneApertureDiameterMicrons) = 1;
    apertureKernel = apertureKernel ./ sum(apertureKernel(:));
end