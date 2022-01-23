function [activationImage, activationImageLMScone, ...
    sampledHexMosaicXaxis, sampledHexMosaicYaxis] = ...
    computeActivationDensityMap(obj, activationFrame)
% Compute activation images for the hex mosaic (all cones + LMS submosaics)
%
% Syntax:
%   [activationImage, activationImageLMScone, sampledHexMosaicXaxis, ...
%       sampledHexMosaicYaxis] = computeActivationDensityMap(obj, 
%       activation)
%
% Description:
%    Compute the activation images for the hex mosaic (all of the cones +
%    LMS sub-mosaics)
%
%    We need to store the activationImage, particularly when it is a
%    movie.  Then we can play it multiple times.
%
%    Maybe it should be stored in the coneMosaic object, or at least
%    returned by the plot routine for the movie.
%
% Inputs:
%    obj                    - The cone hex mosaic object
%    activation             - The activation to map
%
% Outputs:
%    activationImage        - The created activation image
%    activationImageLMScone - The LMS cones for the activation image
%    sampledHexMosaicXaxis  - The X-axis for the hex mosaic
%    sampledHexMosaicYaxis  - The Y-axis for the hex mosaic
%
% Optional key/value pairs:
%    None.
%
% Notes:
%

% History:
%    xx/xx/15  NPC  ISETBIO TEAM, 2015
%    02/16/18  jnm  Formatting
%    4/5/18    NPC  Sped up algorithm
%    5/26/20   NPC  Fixed iRow issue, which was causing the mosaic plotting Y-coord flip 
%                   (rows grow top -> bottom, whereas Y-coords grow bottom -> top)
    
    sampledHexMosaicXaxis = squeeze(obj.patternSupport(1, :, 1)) + ...
        obj.center(1);
    sampledHexMosaicYaxis = squeeze(obj.patternSupport(:, 1, 2)) + ...
        obj.center(2);

    dxOriginal = sampledHexMosaicXaxis(2) - sampledHexMosaicXaxis(1);
    interpolationF = 4;
    dx = dxOriginal / interpolationF;
    sampledHexMosaicXaxis = ...
        sampledHexMosaicXaxis(1):dx:(sampledHexMosaicXaxis(end)+dxOriginal+dx);
    sampledHexMosaicYaxis = ...
        sampledHexMosaicYaxis(1):dx:(sampledHexMosaicYaxis(end)+dxOriginal+dx);

    % Generate cone aperture kernel - flat top Gaussian
    nSamples = ceil(0.5 * obj.pigment.width / dx);
    apertureSupport = (-nSamples:nSamples) * dx;
    apertureSigma = obj.pigment.width / 6;
    [x, y] = meshgrid(apertureSupport, apertureSupport);
    apertureKernel = exp(-0.5 * (x / apertureSigma) .^ 2) .* exp(-0.5 * ...
        (y / apertureSigma) .^ 2);
    apertureKernel = apertureKernel / max(apertureKernel(:));
    apertureKernel = apertureKernel .^ 1.0;
    apertureKernel(apertureKernel < exp(-0.5 * 3 ^ 2)) = 0;
    apertureKernel(apertureKernel > 0.15) = 0.15;
    apertureKernel = apertureKernel / max(apertureKernel(:));

    visualizeConeApertureKernel = ~true;
    if (visualizeConeApertureKernel)
        figure(222);
        clf;
        imagesc(apertureKernel);
        set(gca, 'CLim', [0 1]);
        axis 'image'
        colormap(gray)
        drawnow;
    end

    activeConeIndices= find(obj.pattern > 1);

    activationImageLMScone = zeros(...
        numel(sampledHexMosaicYaxis), ...
        numel(sampledHexMosaicXaxis), 3);
    
    for coneID = 2:4
        activationImageSingleCone = zeros(...
        numel(sampledHexMosaicYaxis), ...
        numel(sampledHexMosaicXaxis));
        % Calculating individual maps, and then summing at the end
        ix = find(obj.pattern(activeConeIndices) == coneID);  % Positions for this cone class
        [r, c] = ind2sub(size(obj.pattern), activeConeIndices(ix));
        for k = 1:numel(r)
            iRow = (r(k)-1)*interpolationF+round(interpolationF/2);
            iCol = (c(k)-1)*interpolationF+round(interpolationF/2);
            activationImageSingleCone(end-iRow,iCol) = activationFrame(ix(k));
        end
         activationImageLMScone(:,:, coneID-1) = conv2(activationImageSingleCone, apertureKernel, 'same');
    end

    activationImage = sum(activationImageLMScone, 3);
end
