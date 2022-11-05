function [theRotatedRF, rotationDegs] = bestHorizontalResolutionRFmap(theRF, rotationDegs)

    if (isempty(rotationDegs))  
        theta = 0:1:179;
        [R,xOffset] = radon(theRF,theta);

        rowIndices = find(abs(xOffset) < 0.2*max(abs(xOffset)));
        xOffset = xOffset(rowIndices);
        R = R(rowIndices,:);

        [~,idx] = max(R(:));
        [row, col] = ind2sub(size(R), idx)
        rotationDegs = -theta(col);

        debugRadon = false;
        if (debugRadon)
            theRotatedRF = imrotate(theRF, rotationDegs, 'bilinear', 'crop');
            theRotatedRF2 = imrotate(theRF, -rotationDegs, 'bilinear', 'crop');

            figure(1999)
            subplot(2,2,1);
            imagesc(theRF); axis 'image'
        
            subplot(2,2,2)
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