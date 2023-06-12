function [centeredOTF,  translationVector, centeredPSF, xGridMinutes,yGridMinutes] = otfWithZeroCenteredPSF(OTF, PSF, translationVector,  xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showTranslation)
    if (isempty(translationVector))
        % Compute center of mass
        centerOfMass = computeCenterOfMass(PSF, 'useMaxAsCenterOfMass', true);
        centerPosition = floor(size(PSF,1)/2) + 1 * [1 1];
        
        if (showTranslation)
            figure(200);clf
            imagesc(1:size(PSF,2), 1:size(PSF,1), PSF);
            hold on;
            plot(centerPosition(1)*[1 1], [1 size(PSF,1)], 'k-');
            plot([1 size(PSF,2)], centerPosition(2)*[1 1], 'k-');
            plot(centerOfMass(1), centerOfMass(2), 'ro', 'MarkerFaceColor', [1 0.5 0.5]);
            plot([centerPosition(1) centerOfMass(1)], [centerPosition(2)  centerOfMass(2)], 'r-');
            hold off;
            axis 'square'
            colormap(gray(1024));
            drawnow;
        end
        
        % Compute translation vector to bring center of mass at 0,0
        translationVector = (centerPosition-centerOfMass);
    end
    centeredOTF = shiftInFourierPlane(OTF, translationVector);

    [xGridMinutes,yGridMinutes,centeredPSF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(centeredOTF));
    
    if (abs(sum(centeredPSF(:))-1) > 1000*eps(1))
        fprintf('centeredPSF min = %1.9f, max = %1.9f, sum = %1.9f\n', min(centeredPSF(:)), max(centeredPSF(:)), sum(centeredPSF(:)));
        error('PSF volume does not equal 1\n');
    end
end

function centerOfMass = computeCenterOfMass(PSF, varargin)
    p = inputParser;
    p.addParameter('useMaxAsCenterOfMass', true, @islogical);
    p.parse(varargin{:});
    
    if (p.Results.useMaxAsCenterOfMass)
        [~,idx] = max(PSF(:));
        [i,j] = ind2sub(size(PSF), idx);
        centerOfMass = [j i];
    else
        binaryImage = true(size(PSF));
        measurements = regionprops(bwlabel(binaryImage), PSF, 'WeightedCentroid');
        centerOfMass = measurements.WeightedCentroid;
    end
end