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

    sampledHexMosaicXaxis = squeeze(obj.patternSupport(1, :, 1)) + ...
        obj.center(1);
    sampledHexMosaicYaxis = squeeze(obj.patternSupport(:, 1, 2)) + ...
        obj.center(2);

    dxOriginal = sampledHexMosaicXaxis(2) - sampledHexMosaicXaxis(1);
    interpolationF = 4;
    dx = dxOriginal / interpolationF;
    sampledHexMosaicXaxis = ...
        sampledHexMosaicXaxis(1):dx:(sampledHexMosaicXaxis(end)+dxOriginal);
    sampledHexMosaicYaxis = ...
        sampledHexMosaicYaxis(1):dx:(sampledHexMosaicYaxis(end)+dxOriginal);

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
            rr = (r(k)-1)*interpolationF+round(interpolationF/2);
            cc = (c(k)-1)*interpolationF+round(interpolationF/2);
            activationImageSingleCone(rr,cc) = activationFrame(ix(k));
        end
         activationImageLMScone(:,:, coneID-1) = conv2(activationImageSingleCone, apertureKernel, 'same');
    end

    activationImage = sum(activationImageLMScone, 3);
end



%     size(obj.patternSupport)
%     obj.visualizeGrid();
%     pause
%     zeroFrame = zeros(interpolationF * obj.rows, ...
%         interpolationF * obj.cols);
% 
%     [size(zeroFrame) obj.rows obj.cols]
%     
%     % [Row, Col, 3]
%     activationImageLMScone = zeros(interpolationF * obj.rows, ...
%         interpolationF * obj.cols, 3);
% 
%     % First is blank, 2, 3, 4 are L, M, S
%     % Not sure what the algorithm is doing here.
%     % But it is very slow when we have a movie (e.g., 100 temporal samples)
%     % Must rethink. (TODO)
%     for coneID = 2:4
%         % Calculating individual maps, and then summing at the end
%         fprintf('Cone ID %d\n', coneID);
%         ix = find(obj.pattern == coneID);  % Positions for this cone class
%         frame = zeroFrame;
%         [r, c] = ind2sub(size(obj.pattern), ix);
%         [numel(r) numel(c)]
%         for k = 1:numel(ix)  % For each cone position in this cone class
%             yy = (r(k) - 1) * interpolationF + interpolationF / 2 + ...
%                 round(apertureSupport / dx);
%             xx = (c(k) - 1) * interpolationF + interpolationF / 2 + ...
%                 round(apertureSupport / dx);
%             xx = xx(xx > 0 & xx <= size(frame, 2));
%             yy = yy(yy > 0 & yy <= size(frame, 1));
%             [xxx, yyy] = meshgrid(xx, yy);
%             xxx = xxx(:);
%             yyy = yyy(:);
%             [numel(r) numel(c) k]
%             yyyy = yyy - ((r(k) - 1) * interpolationF + ...
%                 interpolationF / 2) + 1 + nSamples;
%             xxxx = xxx - ((c(k) - 1) * interpolationF + ...
%                 interpolationF / 2) + 1 + nSamples;
%             frame(yyy, xxx) = frame(yyy, xxx) + ...
%                 activation(r(k), c(k)) * apertureKernel(yyyy, xxxx);
%         end
%         % frame has been interpolated to the whole image size.
%         activationImageLMScone(:, :, coneID - 1) = frame;
%     end
%     fprintf('Combining cone maps\n');
%     activationImage = sum(activationImageLMScone, 3);
% end
