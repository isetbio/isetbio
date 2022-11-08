function [theRotatedRF, rotationDegs] = bestHorizontalResolutionRFmap(theRF, rotationDegs)

    if (isempty(rotationDegs))  
        % Sample the rotation axis (theRF rotation) every 1 degree
        theta = 0:1:179;

        % Do the Radon transform
        [R,xOffset] = radon(theRF,theta);

        % Search only the central offsets
        rowIndices = find(abs(xOffset) < 0.2*max(abs(xOffset)));
        R = R(rowIndices,:);

        % Find peak of Radon transform 
        [~,idx] = max(R(:));
        [~, col] = ind2sub(size(R), idx);

        % Peak rotation is theta(col), so rotate by its negative
        rotationDegs = -theta(col);

        % Debug the radon transform information
        debugRadon = false;
        if (debugRadon)
            theRotatedRF = imrotate(theRF, rotationDegs, 'bilinear', 'crop');
            theRotatedRF2 = imrotate(theRF, -rotationDegs, 'bilinear', 'crop');

            figure(1999)
            subplot(2,2,1);
            imagesc(theRF); axis 'image'
        
            subplot(2,2,2);
            xOffset = xOffset(rowIndices);
            imshow(R,[],'Xdata',theta,'Ydata',xOffset,'InitialMagnification','fit')
            set(gca, 'XTick', theta)
            xlabel('\theta (degrees)')
            ylabel('offset (pixels from center)')
            colormap(gca,hot), colorbar
        
            subplot(2,2,3);
            imagesc(theRotatedRF); axis 'image'
        
            subplot(2,2,4);
            imagesc(theRotatedRF2); axis 'image'
            drawnow
        end
   end

   theRotatedRF = imrotate(theRF, rotationDegs, 'bilinear', 'crop');
end