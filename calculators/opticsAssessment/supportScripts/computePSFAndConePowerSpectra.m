function [psfImageSpectrum, psfImageSpectrumRoatated, coneImageSpectrum, ...
    spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, coneSpectrumSlice, ...
    spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, ...
    spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, ...
    effectivePSFSpectrumSliceX, effectivePSFSpectrumSliceY, ...
    coneCutoffSF, psfXCutoffSF, psfYCutoffSF] = ...
        computePSFAndConePowerSpectra(psfImage, coneImage, psfSupport)
    

    fftSize = 1024;
    psfImageSpectrum = fftshift(abs(fft2(psfImage, fftSize,fftSize)));
    psfImageSpectrum = psfImageSpectrum / max(psfImageSpectrum(:));
    
    coneImageSpectrum = fftshift(abs(fft2(coneImage, fftSize,fftSize)));
    coneImageSpectrum = coneImageSpectrum / max(coneImageSpectrum(:));
    
    binaryImage = psfImageSpectrum;
    m1 = min(binaryImage(:));
    m2 = max(binaryImage(:));
    binaryImage = imbinarize((binaryImage - m1)/(m2-m1));
    s = regionprops(binaryImage,'orientation');
    psfRotation = s.Orientation;
    
    psfImageSpectrumRoatated = imrotate(psfImageSpectrum, -psfRotation, 'bilinear', 'crop');
    psfImageSpectrumRoatated = psfImageSpectrumRoatated / max(psfImageSpectrumRoatated(:));
    
    dx = psfSupport(2)-psfSupport(1);
    maxSF = 1/(2*dx);
    deltaSF = maxSF / (fftSize/2);
    spectralSupportCyclesPerDegree = ((-maxSF+deltaSF):deltaSF:maxSF) * 60;
    
    N = size(psfImageSpectrumRoatated,1);
    midPoint = [N/2 N/2];
    spectralSupportCyclesPerDegreePositive = spectralSupportCyclesPerDegree(midPoint(1):end);
    
    
    [maxPSF,maxPoint] = max(psfImageSpectrumRoatated(:),[], 1);
    [maxPoint(1), maxPoint(2)] = ind2sub(size(psfImageSpectrumRoatated), maxPoint);

%     if (midPoint(1) ~= maxPoint(1))
%         fprintf('midPoint (%d) differs from maxPointY (%d)\n', midPoint(1),maxPoint(1));
%     end
%     if (midPoint(2) ~= maxPoint(2))
%         fprintf('midPoint (%d) differs from maxPointX (%d)\n', midPoint(1),maxPoint(2));
%     end
    %midPoint = maxPoint;
    

    spectralSupportCyclesPerDegreePositiveY = spectralSupportCyclesPerDegree(midPoint(1):end-1);
    psfSpectrumSliceY = psfImageSpectrumRoatated(midPoint(1):end-1,midPoint(2));
    
    spectralSupportCyclesPerDegreePositiveX = spectralSupportCyclesPerDegree(midPoint(2):end-1);
    psfSpectrumSliceX = psfImageSpectrumRoatated(midPoint(1), midPoint(2):end-1);
    
    coneSpectrumSlice = coneImageSpectrum(midPoint(1), midPoint(2):end-1);
    

    effectivePSFSpectrumSliceX = squeeze(coneSpectrumSlice) .* squeeze(psfSpectrumSliceX);
    effectivePSFSpectrumSliceY = squeeze(coneSpectrumSlice) .* squeeze(psfSpectrumSliceY');
    
    % -15 dB in amplitude. At -15dB, power is 3.16%, amplitude is 17.78%
    cornerFrequencyAttenuationDB = -15;
    cornerFrequencyAttenuation = 10^(cornerFrequencyAttenuationDB/20);
    
    
    k = numel(effectivePSFSpectrumSliceX);
    notFound = true;
    psfXCutoffSF = [];
    while ((k > 1)&&(notFound))
        if (effectivePSFSpectrumSliceX(k) < cornerFrequencyAttenuation) && (effectivePSFSpectrumSliceX(k-1) >= cornerFrequencyAttenuation)
            notFound = false;
            idx = [k-1 k];
            d = [abs(effectivePSFSpectrumSliceX(idx(1))-cornerFrequencyAttenuation) abs(effectivePSFSpectrumSliceX(idx(2))-cornerFrequencyAttenuation)];
            w(1) = d(2)/(d(1)+d(2));
            w(2) = d(1)/(d(1)+d(2));
            psfXCutoffSF = spectralSupportCyclesPerDegreePositiveX(idx(1))*w(1) + spectralSupportCyclesPerDegreePositiveX(idx(2))*w(2);
        end
        k = k-1;
    end
    if (isempty(psfXCutoffSF))
        [~, idx] = min(abs(effectivePSFSpectrumSliceX-cornerFrequencyAttenuation));
        psfXCutoffSF = spectralSupportCyclesPerDegreePositiveX(idx);
    end
    
    k = numel(effectivePSFSpectrumSliceY);
    notFound = true;
    psfYCutoffSF = [];
    while (( k> 1)&&(notFound))
        if (effectivePSFSpectrumSliceY(k) < cornerFrequencyAttenuation) && (effectivePSFSpectrumSliceY(k-1) >= cornerFrequencyAttenuation)
            notFound = false;
            idx = [k-1 k];
            d = [abs(effectivePSFSpectrumSliceY(idx(1))-cornerFrequencyAttenuation) abs(effectivePSFSpectrumSliceY(idx(2))-cornerFrequencyAttenuation)];
            w(1) = d(2)/(d(1)+d(2));
            w(2) = d(1)/(d(1)+d(2));
            psfYCutoffSF = spectralSupportCyclesPerDegreePositiveY(idx(1))*w(1) + spectralSupportCyclesPerDegreePositiveY(idx(2))*w(2);
        end
        k = k-1;
    end
    if (isempty(psfYCutoffSF))
        [~, idx] = min(abs(effectivePSFSpectrumSliceY-cornerFrequencyAttenuation));
        psfYCutoffSF = spectralSupportCyclesPerDegreePositiveY(idx);
    end
    
    
    k = 1;
    notFound = true;
    coneCutoffSF = [];
    while(k < numel(coneSpectrumSlice))&&notFound
        if (coneSpectrumSlice(k) > cornerFrequencyAttenuation) && (coneSpectrumSlice(k+1) <= cornerFrequencyAttenuation)
            notFound = false;
            idx = [k k+1];
            d = [abs(coneSpectrumSlice(idx(1))-cornerFrequencyAttenuation) abs(coneSpectrumSlice(idx(2))-cornerFrequencyAttenuation)];
            w(1) = d(2)/(d(1)+d(2));
            w(2) = d(1)/(d(1)+d(2));
            coneCutoffSF = spectralSupportCyclesPerDegreePositiveX(idx(1))*w(1) + spectralSupportCyclesPerDegreePositiveX(idx(2))*w(2);
        end
        k = k+1;
    end
end