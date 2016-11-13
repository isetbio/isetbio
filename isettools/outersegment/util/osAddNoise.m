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
%   sampTime - sample time interval (secs), see sensorGet(s, 'time interval');
%   NYI:  seed     - noise seed for reproducibility
%
%  Outputs:
%    adaptedCur - membrane current with noise added
%    freq       - temporal frequency for checking the noise
%
%  Example:
%    
%    [noise, freq] = osAddNoise(zeros(10,10,10000), 'sampTime',1/5000);
%    noiseF = squeeze(mean(mean(abs(fft(noise, [], 3)).^2)));
%    vcNewGraphWin; loglog(freq, noiseF(1:length(freq)));
%    xlabel('Frequency(Hz)'); ylabel('Power Spectrum (pA^2/Hz)'); 
%
%   These y units are way off, it seems to BW.
%      xlim([1 1e3]); ylim([1e-4 1]);
%
%  See also:
%    coneAdapt, osAdaptSteadyState, osAdaptTemporal
%
%  (HJ) ISETBIO, 2014

%% Init
p = inputParser;
p.addRequired('curNF',@isnumeric);
p.addParameter('sampTime',0.001,@isscalar);
p.addParameter('seed',rng,@isstruct);
p.parse(curNF,varargin{:});

sampTime = p.Results.sampTime;

% We do not seem to be handling this correctly!
% seed = p.Results.seed;

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
