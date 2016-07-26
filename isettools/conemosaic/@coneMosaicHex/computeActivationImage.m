function [activationImage, activationImageLMScone, sampledHexMosaicXaxis, sampledHexMosaicYaxis] = computeActivationImage(obj, activation)
% Compute activation images for the hex mosaic (all cones +  LMS submosaics)
%
% NPC, ISETBIO TEAM, 2015

    sampledHexMosaicXaxis = squeeze(obj.patternSupport(1,:,1));
    sampledHexMosaicYaxis = squeeze(obj.patternSupport(:,1,2));
    
    dx = sampledHexMosaicXaxis(2)-sampledHexMosaicXaxis(1);
    interpolationF = 2; dx = dx / interpolationF;
    sampledHexMosaicXaxis = sampledHexMosaicXaxis(1):dx:sampledHexMosaicXaxis(end);
    sampledHexMosaicYaxis = sampledHexMosaicYaxis(1):dx:sampledHexMosaicYaxis(end);
    
    % Generate cone aperture kernel - flat top Gaussian
    nSamples = ceil(0.5*obj.pigment.width/dx);
    apertureSupport = (-nSamples:nSamples)*dx;
    apertureSigma = obj.pigment.width/6;
    [x,y] = meshgrid(apertureSupport, apertureSupport);
    apertureKernel = exp(-0.5*(x/apertureSigma).^2) .* exp(-0.5*(y/apertureSigma).^2);
    apertureKernel = apertureKernel / max(apertureKernel(:));
    apertureKernel = apertureKernel .^ 0.7;
    apertureKernel(apertureKernel>0.15) = 0.15;
    apertureKernel = apertureKernel / max(apertureKernel(:));
    
    visualizeConeApertureKernel = false;
    if (visualizeConeApertureKernel)
        figure(222); clf;
        imagesc(apertureKernel);
        set(gca, 'CLim', [0 1]);
        axis 'image'
        colormap(gray)
        drawnow;
        pause
    end
    
    zeroFrame = zeros(interpolationF*obj.rows, interpolationF*obj.cols);
    activationImageLMScone = zeros(interpolationF*obj.rows, interpolationF*obj.cols, 3);
    
    for coneID = 2:4
        ix = find(obj.pattern == coneID);
        for k = 1:numel(ix)
            [r,c] = ind2sub(size(obj.pattern),ix(k));
            frame = zeroFrame;
            frame((r-1)*interpolationF + interpolationF/2,(c-1)*interpolationF + interpolationF/2) = activation(r,c);
            frame = conv2(frame, apertureKernel, 'same');
            activationImageLMScone(:,:,coneID-1) = activationImageLMScone(:,:,coneID-1) + frame;
        end
    end
    activationImage = sum(activationImageLMScone, 3);
end
