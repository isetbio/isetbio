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
        [row, col] = ind2sub(size(R), idx);
        
        % Peak rotation is theta(col), so rotate by its negative
        rotationDegs = -theta(col);

        % Debug the radon transform
        debugRadonTransformAnalysis = true;
    else
        debugRadonTransformAnalysis = false;
    end

    % The rotatedRF so that the elongation is along the y-axis
    theRotatedRF = imrotate(theRF, rotationDegs, 'bilinear', 'crop');

    % Debug the radon transform information
    if (debugRadonTransformAnalysis)
        hFig = figure(1999);
        set(hFig, 'Position', [10 10 1500 530]);
        subplot(1,3,1);
        imagesc(theRF); axis 'image'
    
        subplot(1,3,2);
        xOffset = xOffset(rowIndices);
        imshow(R,[],'Xdata',theta,'Ydata',xOffset,'InitialMagnification','fit');
        hold on;
        plot(theta(col), xOffset(row), 'kx');
        set(gca, 'XTick', theta, 'XTickLabel', sprintf('%2.0f\n',theta));
        xlabel('\theta (degrees)')
        ylabel('offset (pixels from center)')
        colormap(gca,hot), colorbar
    
        subplot(1,3,3);
        imagesc(theRotatedRF); axis 'image'
        title(sprintf('best orientation = %2.0f degs', rotationDegs));

        drawnow;
    end
   
end