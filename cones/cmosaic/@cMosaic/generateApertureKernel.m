function apertureKernel = generateApertureKernel(obj, blurApertureDiameterMicrons, oiResMicrons, lowOpticalImageResolutionWarning)

    % Compute convolution kernel size in pixels - make it big here, we'll
    % trim at the end
    apertureSamples = ceil(blurApertureDiameterMicrons*2 / oiResMicrons);
    
    if (apertureSamples < 3)
        if (lowOpticalImageResolutionWarning)
            fprintf(2, 'Warning: The optical image resolution (%2.2f microns) is too low relative to the cone aperture (%2.2f microns) for accurate computation of cone blur.\n\tSkipping blur by cone aperture.\n\tConsider increasing the # of pixels in the stimulus or decrease the stimulus FOV.\n', oiResMicrons, blurApertureDiameterMicrons);
        end
        % Return delta function aperture
        apertureKernel = 1;
        return;
    end
    
    % Make sure it is odd
    if (mod(apertureSamples, 2) == 0)
        apertureSamples = apertureSamples + 1;
    end
    
    %fprintf('Aperture kernel is %d x %d pixels\n', apertureSamples, apertureSamples);
    %fprintf('Aperture diameter: %2.2f microns, optical image resolution: %2.2f microns\n', blurApertureDiameterMicrons , oiResMicrons);
    
    % Generate aperture kernel
    [xx,yy] = sample2space(1:apertureSamples, 1:apertureSamples, oiResMicrons, oiResMicrons);
    [apertureRows, apertureCols] = meshgrid(xx,yy);
        
    if (isfield(obj.coneApertureModifiers, 'shape')) 
        %fprintf('Aperture shape is %s\n', obj.coneApertureModifiers.apertureShape)
        switch (obj.coneApertureModifiers.shape)
            case 'Gaussian'
                % Gaussian with sigma = obj.coneApertureModifiers.sigma
                blurApertureSigma = obj.coneApertureModifiers.sigma * blurApertureDiameterMicrons;
                [apertureKernel, theoreticalAreaMetersSquared, theoreticalAreaMicronsSquared, actualAreaMicronsSquared] = ...
                    gaussianKernel(apertureRows, apertureCols, blurApertureSigma);
            case 'Pillbox'
                % Pillbox with unit volume
                [apertureKernel, theoreticalAreaMetersSquared, theoreticalAreaMicronsSquared, actualAreaMicronsSquared] = pillboxKernel(apertureRows, apertureCols, blurApertureDiameterMicrons);
            otherwise
                error('Do not know how to generate a ''%s'' aperture.', obj.coneApertureModifiers.apertureShape)
        end
    else
        % Default aperture: pillbox with unit volume
        [apertureKernel, theoreticalAreaMetersSquared, theoreticalAreaMicronsSquared, actualAreaMicronsSquared] = pillboxKernel(apertureRows, apertureCols, blurApertureDiameterMicrons);
    end
    
    % Keep pixels up to 1% of the peak amplitude
    minAmplitudeKept = 0.01 * max(apertureKernel(:));
    idx = find(apertureKernel>= minAmplitudeKept);
    [keptRows, keptCols] = ind2sub(size(apertureKernel), idx);
    rr = [min(keptRows) max(keptRows)];
    cc = [min(keptCols) max(keptCols)];
        
    % Crop it symmetrically
    midPoint = (apertureSamples-1)/2+1;
  
    if ((rr(2)-rr(1)) >= (cc(2)-cc(1)))
        rr = (midPoint - (midPoint-rr(1))):1:(midPoint + (midPoint-rr(1)));
    else
        rr = (midPoint - (midPoint-cc(1))):1:(midPoint + (midPoint-cc(1)));
    end

    debugAperture = ~true;
    if (debugAperture)
        
        figure(2); clf;
        subplot(1,2,1)
        xAxis = apertureCols(:,1);
        yAxis = apertureRows(1,:);
        imagesc(xAxis, yAxis,apertureKernel);
        set(gca, 'XTick', [xAxis(1) 0 xAxis(end)], 'YTick',  [yAxis(1) 0 yAxis(end)]);
        axis 'image';
        colormap(brewermap(1024, '*greys'));

        apertureKernel = apertureKernel(rr,rr);
        apertureRows = apertureRows(rr,rr);
        apertureCols = apertureCols(cc,cc);

        subplot(1,2,2)
        xAxis = apertureCols(:,1);
        yAxis = apertureRows(1,:);
        imagesc(xAxis, yAxis,apertureKernel);
        set(gca, 'XTick', [xAxis(1) 0 xAxis(end)], 'YTick',  [yAxis(1) 0 yAxis(end)]);
        axis 'image';
        colormap(brewermap(1024, '*greys'));
        
        
        m = (size(apertureKernel,1)-1)/2+1;
        dx = 0*(xx(2)-xx(1))/2;
        figure(12345);
        clf;
        subplot(1,2,1);
        imagesc(xx,yy,apertureKernel);
        hold on;
        stem(xx-dx, min(yy) + (max(yy)-min(yy))*0.95*apertureKernel(m,:)/max(apertureKernel(:)), ...
            'filled', 'Basevalue', min(yy), 'Color', 'k', 'LineWidth', 1.5);
        hold off
        set(gca, 'XTick', -5:1:5, 'YTick', -5:1:5);
        axis 'image'; axis 'xy';
        xlabel('microns');
        title('employed kernel');
        
        subplot(1,2,2);
        stem(xx-dx,apertureKernel(m,:), 'filled', 'Basevalue', 0, 'Color', 'r', 'LineWidth', 1.5); hold on;
        set(gca, 'XTick', -5:1:5);
        xlabel('microns');
        colormap(brewermap(1024, 'greys'));
        drawnow
        pause
    end
    
end

function [apertureKernel, theoreticalAreaMetersSquared, theoreticalAreaMicronsSquared, actualAreaMicronsSquared] = pillboxKernel(apertureRows, apertureCols, blurApertureDiameterMicrons)
    
    % Generate aperture kernel
    apertureRadii = sqrt(apertureRows .^ 2 + apertureCols .^ 2);
    apertureKernel = apertureRadii*0;
    apertureKernel(apertureRadii <= 0.5*blurApertureDiameterMicrons) = 1;
    
    % Unit volume
    actualAreaMicronsSquared = sum(apertureKernel(:));
    apertureKernel = apertureKernel ./ actualAreaMicronsSquared;
    
    % Theoretical area based on continuous space
    theoreticalAreaMetersSquared = ((pi * (0.5*blurApertureDiameterMicrons*1e-6).^2));
    theoreticalAreaMicronsSquared = theoreticalAreaMetersSquared * 1e12;
    
    % Actual area based on discritized space
    dx = apertureCols(2)-apertureCols(1);
    actualAreaMicronsSquared = actualAreaMicronsSquared * dx^2;
end

function [apertureKernel, theoreticalAreaMetersSquared, theoreticalAreaMicronsSquared, actualAreaMicronsSquared] = ...
    gaussianKernel(apertureRows, apertureCols, gaussianSigma)

    % Generate aperture kernel
    apertureRadii = sqrt(apertureRows .^ 2 + apertureCols .^ 2);
    apertureKernel = exp(-0.5*(apertureRadii/gaussianSigma).^2);
    
    % Unit volume
    actualAreaMicronsSquared = sum(apertureKernel(:));
    apertureKernel = apertureKernel ./ actualAreaMicronsSquared;
    
    % Theoretical area based on continous space
    characteristicRadiusMicrons = gaussianSigma * sqrt(2.0);
    theoreticalAreaMetersSquared = ((pi * (characteristicRadiusMicrons*1e-6).^2));
    theoreticalAreaMicronsSquared = theoreticalAreaMetersSquared * 1e12;
    
    %Actual area based on discritized space
    dx = apertureCols(2)-apertureCols(1);
    actualAreaMicronsSquared = actualAreaMicronsSquared * dx^2;
end
