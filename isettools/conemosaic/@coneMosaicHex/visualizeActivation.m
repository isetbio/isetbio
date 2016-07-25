function visualizeActivation(obj, activation)

    [activationImage, sampledHexMosaicXaxis, sampledHexMosaicYaxis] = obj.computeActivationImage(activation);
    figure();
    imagesc(sampledHexMosaicXaxis, sampledHexMosaicYaxis, activationMap);
    axis 'image'; axis 'xy';
    colormap(gray(1024));
    drawnow
end
