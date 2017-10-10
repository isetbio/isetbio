function [activationImage, activationImageLMScone, sampledHexMosaicXaxis, sampledHexMosaicYaxis] = computeActivationDensityMap(obj, activation)
% Compute activation images for the hex mosaic (all cones +  LMS submosaics)
%
% The variable names could be improved (xxxx).
%
% We need to store the activationImage, particularly when it is a
% movie.  Then we can play it multiple times.
%
% Maybe it should be stored in the coneMosaic object, or at least
% returned by the plot routine for the movie.
%
% Good start.  Other 
% NPC, ISETBIO TEAM, 2015

    sampledHexMosaicXaxis = squeeze(obj.patternSupport(1,:,1)) + obj.center(1);
    sampledHexMosaicYaxis = squeeze(obj.patternSupport(:,1,2)) + obj.center(2);
    
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
    apertureKernel = apertureKernel .^ 1.0;
    apertureKernel(apertureKernel<exp(-0.5*(3)^2)) = 0;
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
    
    % [Row,Col, 3]
    activationImageLMScone = zeros(interpolationF*obj.rows, interpolationF*obj.cols, 3);

    % First is blank, 2,3,4 are L,M,S
    % Not sure what the algorithm is doing here.
    % But it is very slow when we have a movie (e.g., 100 temporal samples).
    % Must rethink.
    for coneID = 2:4
        % Calculating individual maps, and then summing at the end
        fprintf('Cone ID %d\n',coneID);
        ix = find(obj.pattern == coneID);  % Positions for this cone class
        frame = zeroFrame;
        [r,c] = ind2sub(size(obj.pattern),ix);
        for k = 1:numel(ix)  % For each cone position in this cone class
            yy = (r(k)-1)*interpolationF + interpolationF/2 + round(apertureSupport/dx);
            xx = (c(k)-1)*interpolationF + interpolationF/2 + round(apertureSupport/dx);
            xx = xx(xx>0 & xx <= size(frame,2));
            yy = yy(yy>0 & yy <= size(frame,1));
            [xxx,yyy] = meshgrid(xx,yy);
            xxx = xxx(:);
            yyy = yyy(:);
            yyyy = yyy-((r(k)-1)*interpolationF + interpolationF/2)+1 + nSamples;
            xxxx = xxx-((c(k)-1)*interpolationF + interpolationF/2)+1 + nSamples;
            frame(yyy, xxx) = frame(yyy,xxx) + activation(r(k),c(k))*apertureKernel(yyyy,xxxx);
        end
        % frame has been interpolated to the whole image size.
        activationImageLMScone(:,:,coneID-1) = frame;
    end
    fprintf('Combining cone maps\n');
    activationImage = sum(activationImageLMScone, 3);
end
