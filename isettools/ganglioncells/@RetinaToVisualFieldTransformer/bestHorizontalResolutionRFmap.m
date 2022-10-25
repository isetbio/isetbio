function [theRotatedRF, rotationDegs] = bestHorizontalResolutionRFmap(theRF, rotationDegs)
    % Determine rotation that maximizes horizontal resolution
    if (isempty(rotationDegs))
        binaryImage = theRF;
        m1 = min(binaryImage(:));
        m2 = max(binaryImage(:));
        binaryImage = imbinarize((binaryImage - m1)/(m2-m1));
        s = regionprops('table', binaryImage,'Orientation');
        rotationDegs = s.Orientation(1);
        fprintf(2,'RF center horizontal resolution maximizing rotation: %2.1f degs\n', rotationDegs);
    end

    % Rotate RF so as to maximize horizontal resolution
    theRotatedRF = imrotate(theRF, 90-rotationDegs, 'bilinear', 'crop');

    debugRotation = false;
    if (debugRotation)
        figure(999); clf;
        subplot(1,2,1);
        imagesc(theRF);
        axis 'image';
    
        subplot(1,2,2);
        imagesc(theRotatedRF)
        axis 'image';
    
        colormap(gray(1024));
        pause
    end

end