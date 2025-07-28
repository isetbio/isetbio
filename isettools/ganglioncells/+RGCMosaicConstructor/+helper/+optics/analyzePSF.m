function [peakAmplitude, wavelengthIndexOfPeakAmplitude, wavelengthOfPeakAmplitude] = ...
	analyzePSF(thePSF, wavelengthIndexOfPeakAmplitude)
	
	if (isempty(wavelengthIndexOfPeakAmplitude))
		wavesNum = size(thePSF.data,3);
		wavePeakAmplitude = zeros(1, wavesNum);
		for iWave = 1:wavesNum
			the2DPSF = squeeze(thePSF.data(:,:,iWave));
			wavePeakAmplitude(iWave) = max(the2DPSF(:));
		end
		[peakAmplitude, wavelengthIndexOfPeakAmplitude] = max(wavePeakAmplitude);
	else
		the2DPSF = squeeze(thePSF.data(:,:,wavelengthIndexOfPeakAmplitude));
		peakAmplitude = max(the2DPSF(:));
	end
	wavelengthOfPeakAmplitude = thePSF.supportWavelength(wavelengthIndexOfPeakAmplitude);

	visualizeThePSFslices = false;
	if (visualizeThePSFslices)
		for iWave = 1:wavesNum
			the2DPSF = squeeze(thePSF.data(:,:,iWave));
			figure(iWave);
			imagesc(the2DPSF, [0 peakAmplitude])
			axis 'image'
			title(sprintf('%2.0fnm (max: %f/%f)', thePSF.supportWavelength(iWave), max(the2DPSF(:)), peakAmplitude));
			drawnow
		end
	end
end