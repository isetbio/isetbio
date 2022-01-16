function [subSampledWavelengthSampling, subSampledSPDs] = SubSampleSPDs(originalS, originalSPDs, targetS, lowPassSigma, showFig)
% Subsample spectral power distributions.
% 
% Synopsis:
%    [subSampledWavelengthSampling, subSampledSPDs] = ptb.SubSampleSPDs(originalS, originalSPDs, targetS, lowPassSigma, showFig)
%
% Description:
%   Subsample the SPDs by a given sampling interval (given in nanometers)
%   after first low-passing them with a Gaussian kernel with sigma =
%   lowPassSigma (given in nanometers).
%
%   This routine also converts power units from PTB convention of power per
%   wavelength band to ISETBio convention of power per nm. That happens in
%   the call to SplineSpd when the input is splined to 1 nm interval.
%
%   If the showFig flag is set to true a figure showing all 
%   the intermediate steps of this operation is displayed.
%
% Inputs:
%    originalS           - Wavelength sampling of input as S vector
%                          [startWl deltaWl nWlSamples].
%    originalSpds        - Input spds in columns of a matrx.
%    targetS             - New wavelength sampling as S vector. New
%                          wavelengths must be within range of original
%                          wavelengths, but can specify a shorter
%                          wavelength range. The deltaWl is rounded to an
%                          integer, which could cause all kinds of pain if
%                          you don't pass it in as an integer.
%    lowPassSigma        - SD of Gaussian kernal in nm.
%    showFig             - True for a plot, false otherwise.
%
% Outputs:
%    subSampledWavelengthSampling - List of wavelengths (not S vector
%                         format) for the resampled data.
%    subSampledSPDs      - SPDs on targetS, smoothed, and converted to
%                          units of power per nm.
%
% Optional key/value pairs:
%   None.
%
% See also: ptb.GenerateIsetbioDisplayObjectFromPTBCalStruct,
%           ptb.ExtraCalData.

% History:
%   2/26/2015     npc   Wrote it.
%   3/10/2015     xd    Changed newSamplingInterval input to targetS
%                       This allows targetS to have a different range
%                       that is still contained within orginialS.
%   1/16/22       dhb   Comments, particularly about power conversion.
%                       Also removed rounding of lowPassSigma to an integer,
%                       becuase I couldn't see why this was being done.

    % New delta wavelength
    newSamplingInterval = targetS(2);

    % Ensure that newSamplingInterval, lowPassSigma are integers
    newSamplingInterval = round(newSamplingInterval);

    % This used to be rounded to an integer, but I can't see why that is a
    % good idea. Commented out for now.
    %   lowPassSigma = round(lowPassSigma);
    
    % Interpolate to 1 nm resolution.  The call to SplineSpd converts power
    % from assumed PTB convention of power per wavelength band to power per nm.
    originalWavelengthSampling = SToWls(originalS);
    maxResWavelengthSampling = (originalWavelengthSampling(1):1:originalWavelengthSampling(end))';
    newS = WlsToS(maxResWavelengthSampling);
    maxResSPDs = SplineSpd(originalS, originalSPDs, newS);
     
    % Generate subsampling vector containing the indices of the samples to
    % keep. Do this by finding start and end intervals for targetS and
    % where these wls fall in the 1 nm resolution version.
    targetWavelengthSampling = SToWls(targetS);
    startIndex = targetS(1) - originalS(1) + 1;
    endIndex = targetWavelengthSampling(end) - originalS(1) + 1;
    subSamplingVector = (startIndex:newSamplingInterval:endIndex);
    subSampledWavelengthSampling = maxResWavelengthSampling(subSamplingVector);

    % Check that resampling worked as expected.
    if (any(subSampledWavelengthSampling ~= targetWavelengthSampling))
        error('Something funny about way wavelength sampling was computed');
    end
    
    % Rreallocate memory for the subsampled SPDs
    channelsNum = size(originalSPDs,2);
    subSampledSPDs = zeros(numel(subSampledWavelengthSampling), channelsNum);
    lowpassedSPDs = zeros(numel(maxResWavelengthSampling), channelsNum);
    
    % Generate the lowpass kernel
    lowPassKernel = generateGaussianLowPassKernel(newSamplingInterval, lowPassSigma, maxResWavelengthSampling);
        
    % Zero pad lowpass kernel
    FFTsize = 1024;
    paddedLowPassKernel = zeroPad(lowPassKernel, FFTsize);
    
    if (showFig)
        hFig = figure();
        set(hFig, 'Position', [200 200 1731 1064]);
        clf;
    end

    % Do the subsampling/smoothing
    for channelIndex = 1:channelsNum
        % Zero pad SPD
        paddedSPD = zeroPad(squeeze(maxResSPDs(:, channelIndex)), FFTsize);
        
        % Filter SPD with kernel
        FFTkernel = fft(paddedLowPassKernel);
        FFTspd    = fft(paddedSPD);
        tmp       = FFTspd .* FFTkernel;
            
        % Back in original domain
        tmp = ifftshift(ifft(tmp));
        lowpassedSPDs(:, channelIndex) = extractSignalFromZeroPaddedSignal(tmp, numel(maxResWavelengthSampling));
        
        % Subsample the lowpassed SPD
        subSampledSPDs(:,channelIndex) = lowpassedSPDs(subSamplingVector,channelIndex);
    end

    % Figure, optionally  
    if (showFig)
        % Find max for plotting. Sometimes the passed vector is all 0, so
        % we make something up for that case.
        maxY = max([max(subSampledSPDs(:)) max(originalSPDs(:)) max(lowpassedSPDs(:))]);
        if (maxY == 0)
            maxY = 2e-3;
        end

        % Compare power, taking into account the change in power units for the resampled spectrum.
        originalSPDpower   = sum(originalSPDs,1);
        subSampledSPDpower = sum(subSampledSPDs,1)*targetS(2);

        % Plot results for each resampled spectrum
        for channelIndex = 1:channelsNum
            % Set the subplot
            subplot(3, 7, [1 2 3 4 5 6]+(channelIndex-1)*7);
            hold on;

            % Plot the lowpass kernel as a stem plot
            hStem = stem(maxResWavelengthSampling, maxY/2 + lowPassKernel*maxY/3, 'Color', [0.5 0.5 0.90], 'LineWidth', 1, 'MarkerFaceColor', [0.7 0.7 0.9]);
            hStem.BaseValue = maxY/2;

            % Plot the subSampledSPD in red
            plot(subSampledWavelengthSampling, subSampledSPDs(:,channelIndex), 'ro-', 'MarkerFaceColor', [1.0 0.7 0.7], 'MarkerSize', 14);

            % Plot the lowpass version of the original SPD in gray
            plot(maxResWavelengthSampling, lowpassedSPDs(:, channelIndex), 'k.-', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 8);

            % Plot the the original SPD in black, and put it in units of
            % power per nm.
            plot(originalWavelengthSampling, originalSPDs(:, channelIndex)/originalS(2), 'ks-', 'MarkerFaceColor', [0.1 0.1 0.1], 'MarkerSize', 6);

            % Axis and labels.
            hold off;
            set(gca, 'YLim', [0 maxY], 'XLim', [min(maxResWavelengthSampling) max(maxResWavelengthSampling)]);
            h_legend = legend('lowpass kernel','subsampled SPD', 'lowpassedSPD', 'originalSPD');
            box on;
            xlabel('wavelength (nm)'); ylabel('power (per nm units');
            title(sprintf('power: %2.4f (original SPD) vs. %2.4f (subsampled SPD)', originalSPDpower(channelIndex), subSampledSPDpower(channelIndex)));
        end
        drawnow;
    end
end

% Method to zero pad a vector to desired size
function paddedF = zeroPad(F, padSize)
    ix = floor(numel(F)/2);
    paddedF = zeros(1,padSize);
    paddedF(padSize/2+(-ix:ix)) = F;
end

% Method to extract a signal from a zero-padded version of it
function F = extractSignalFromZeroPaddedSignal(paddedF, desiredSize)
    ix = floor(desiredSize/2);
    xo = numel(paddedF)/2-1;
    F = paddedF(xo+(-ix:-ix+desiredSize-1));
end

% Method to generate a Gaussian LowPass kernel
function gaussF = generateGaussianLowPassKernel(newSamplingInterval, sigma, samplingAxis) 

    samplingAxis = (0:(numel(samplingAxis)-1))-(numel(samplingAxis)/2)+0.5;
    
    if (newSamplingInterval <= 1)
        gaussF = zeros(size(samplingAxis));
        gaussF(samplingAxis == 0) = 1;
    else
        gaussF = exp(-0.5*(samplingAxis/sigma).^2);
        gaussF = gaussF / sum(gaussF);
    end
end