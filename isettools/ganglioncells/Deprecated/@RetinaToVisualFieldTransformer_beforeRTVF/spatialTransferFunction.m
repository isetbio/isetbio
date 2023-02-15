function [oneSidedSpatialFrequencySupport, oneSidedSTF] = spatialTransferFunction(spatialSupportDegs, theRFprofile)

    fftSize = 2^(2+(nextpow2(numel(theRFprofile))));
    dx = spatialSupportDegs(2)-spatialSupportDegs(1);
    maxSF = 1/(2*dx);
    deltaSF = maxSF / (fftSize/2);
    spatialFrequencySupport = ((-maxSF+deltaSF):deltaSF:maxSF);

    theAmplitudeSpectrum = fftshift(abs(fft(theRFprofile, fftSize)));
    [~,zeroFrequencyBin] = min(abs(spatialFrequencySupport));
    oneSidedSpatialFrequencySupport = spatialFrequencySupport(zeroFrequencyBin:end-1) + deltaSF/10;
    oneSidedSTF = theAmplitudeSpectrum(zeroFrequencyBin:end-1);

    % Sample at 32 frequencies
    frequencySamplesNum = 32;
    sfSampling = logspace(log10(oneSidedSpatialFrequencySupport(1)), log10(oneSidedSpatialFrequencySupport(end)), frequencySamplesNum);
    oneSidedSTFSampled = interp1(oneSidedSpatialFrequencySupport, oneSidedSTF, sfSampling, 'linear');

    oneSidedSpatialFrequencySupport = sfSampling(2:(end-1));
    oneSidedSTF = oneSidedSTFSampled(2:(end-1));
end