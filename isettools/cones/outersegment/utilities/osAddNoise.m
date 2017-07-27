function [adaptedCur, freq] = osAddNoise(curNF, varargin)
% Add noise to noise-free membrane current (curNF)
%
%   [adaptedCur, freq] = osAddNoise(curNF, varargun)
%
% Cone noise is additive and independent of cone signal. The noise is
% Gaussian and the spectral power distribution can be characterized by the
% sum of two Lorentz functions.  (Discuss this.  The SPD of a Gaussian is a
% Gaussian, BW thought).
% 
% See Angueyra and Rieke (2013, Nature Neuroscience) for details.
%
% Inputs:
%    curNF  - noise free cone adapted membrane photo current.  The units
%             are pA (picoAmps) (CHECK!).
%   sampTime - sample time interval (secs)
%   **NYI**  seed     - noise seed for reproducibility
%
%  Outputs:
%    adaptedCur - membrane current with noise added
%    freq       - temporal frequency for checking the noise
%
%  Example:
%    nSamp = 10000; deltaT = 1/5000;
%    [noise, freq] = osAddNoise(zeros(10,10,nSamp), 'sampTime',deltaT);
%    fprintf('Should be close to zero:  Mean noise %f (pA)\n',mean(noise(:)));
%
%  This formula for the noise frequency response amplitude (from Fred R)
%  The 2 is there for some negative frequency, and the correction for the
%  number of time samples is the usual FFT correction.
%    noiseF = squeeze(mean(mean(abs(fft(noise, [], 3)*(2/nSamp)).^2)));
%    vcNewGraphWin; loglog(freq, noiseF(1:length(freq)));
%    line('Xdata',[0.1 freq(end)],'Ydata',[.205 .205]);
%    xlabel('Frequency(Hz)'); ylabel('Power Spectrum (pA^2/Hz)'); 
%    
%
%  See also:
%    coneAdapt, osAdaptSteadyState, osAdaptTemporal
%
%  (HJ) ISETBIO, 2014

%% Init
p = inputParser;
p.addRequired('curNF',@isnumeric);
p.addParameter('sampTime',0.001,@isscalar);
p.addParameter('seed', [] ,@isscalar);
p.parse(curNF,varargin{:});

sampTime = p.Results.sampTime;
seed = p.Results.seed;
if (~isempty(seed))
    rng(seed);
end

%% Build model and generate noise

%  Make sure curNF is [row, col, time] 3D matrix
if isvector(curNF), curNF = reshape(curNF, [1 1 length(curNF)]); end

if (ndims(curNF) == 3),     temporalDimIndex = 3;
elseif (ismatrix(curNF)),   temporalDimIndex = 2;
else error('unexpected curNF format\n');
end

% Generate the noise according to the noise spectral distribution
k = ceil((size(curNF, temporalDimIndex)-1)/2);
freq = (0:k)/ sampTime / size(curNF, temporalDimIndex);

LorentzCoeffs = [0.16  55  4;
                 0.045 190 2.5];
noiseSPD = lorentzSum(LorentzCoeffs, freq);

% Make-up the negative frequency part
noiseSPD = [noiseSPD noiseSPD(end:-1:1)];

noiseSPD = noiseSPD(1:size(curNF, temporalDimIndex));
if (temporalDimIndex == 3)
    noiseSPD = reshape(noiseSPD, [1 1 length(noiseSPD)]);
elseif (temporalDimIndex == 2)
    noiseSPD = reshape(noiseSPD, [1 length(noiseSPD)]);
end
% Have a look at the noise in the frequency domain
% vcNewGraphWin;loglog(squeeze(sqrt(noiseSPD)));

% Generate white Gaussian noise
noise = randn(size(curNF));
noiseFFT = fft(noise, [], temporalDimIndex); % / sqrt(size(noise, 3));

% adjust the spectral power distribution of the noise
noiseFFT = bsxfun(@times, noiseFFT, sqrt(noiseSPD));

noise = real(ifft(noiseFFT, [], temporalDimIndex)) / sqrt(2*sampTime); % take real part

% vcNewGraphWin; plot(squeeze(noise));

% add to noise-free signal
adaptedCur = curNF + noise;

end
