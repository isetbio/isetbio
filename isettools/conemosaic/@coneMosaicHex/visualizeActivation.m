function visualizeActivation(obj, activation)
% Visualize activation images for the hex mosaic (all cones +  LMS submosaics)
%
% NPC, ISETBIO TEAM, 2015

    % Compute activation image maps
    [activationImage, activationImageLMScone, sampledHexMosaicXaxis, sampledHexMosaicYaxis] = obj.computeActivationImage(activation);
    
    figure();
    imagesc(sampledHexMosaicXaxis, sampledHexMosaicYaxis, activationMap);
    axis 'image'; axis 'xy';
    colormap(gray(1024));
    drawnow
end
