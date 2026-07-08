function [noisyCurrent, freq] = osAddNoise(noiseFreeCurrent, varargin)
% Add noise to noise-free membrane current
%
% Syntax:
%    [noisyCurrent, freq] = osAddNoise(noiseFreeCurrent, [varargin])
%
% Description:
%	 Cone photocurrent noise is additive and independent of the current.
%	 The noise is Gaussian and the noise spectral power distribution is
%	 characterized by the sum of two Lorentz functions.
%
%	 This function is called by the os object.
% 
%	 See Angueyra and Rieke (2013, Nature Neuroscience) for details.
%
%    There are examples contained in the code. To access, type 'edit
%    osAddNoise.m' into the Command Window.
%
% Inputs:
%    noiseFreeCurrent - Noise free cone adapted membrane photo current. The
%                       units are pA (picoAmps). [Note: DHB - There was a
%                       comment here that said to double check that these
%                       are in fact the units.]
%
%  Outputs:
%    noisyCurrent     - Membrane current with noise added
%
%    freq             - Temporal frequency domain on which noise is
%                       generated. Could be useful for checking the noise
%
% Optional key/value pairs:
%    'sampTime'       - Number. Sample time interval in seconds
%                       (default 0.001).
%
%    'seed'           - Set random number generator seed to this value
%                       (default empty, don't set).
%
% Notes:
%    * [Note: DHB, BW - The SPD of a Gaussian is a Gaussian, so something
%      about the comment above is therefore incorrect. Maybe the noise is
%      additive but not Gaussian? Or maybe we just don't understand
%      something fundamental.]
%
%  See Also:
%    osAdaptSteadyState, t_osNoise
%

% History:
%    xx/xx/14  HJ   ISETBIO, 2014
%    08/07/17  dhb  Make header comments a closer match to reality.
%    02/09/18  jnm  Formatting

% Examples:
%{
    nSamp = 10000;
    deltaT = 1 / 5000;
	[noise, freq] = osAddNoise(zeros(10, 10, nSamp), 'sampTime', deltaT);
	fprintf('Should be close to zero:  Mean noise %f (pA)\n', ...
        mean(noise(:)));
%}

%% Init
p = inputParser;
p.addRequired('noiseFreeCurrent', @isnumeric);
p.addParameter('sampTime', 0.001, @isscalar);
p.addParameter('seed', [] , @isscalar);
p.parse(noiseFreeCurrent, varargin{:});
sampTime = p.Results.sampTime;

%% Freeze rng seed if desired.
seed = p.Results.seed;
if (~isempty(seed)), rng(seed); end

%% Build model and generate noise
%  Make sure noiseFreeCurrent is [row, col, time] 3D matrix
if isvector(noiseFreeCurrent)
    noiseFreeCurrent = reshape(noiseFreeCurrent, ...
        [1 1 length(noiseFreeCurrent)]);
end

if (ndims(noiseFreeCurrent) == 3)
    temporalDimIndex = 3;
elseif (ismatrix(noiseFreeCurrent))
    temporalDimIndex = 2;
else
    error('Unexpected noiseFreeCurrent format\n');
end

% Generate the noise according to the noise spectral distribution
k = ceil((size(noiseFreeCurrent, temporalDimIndex) - 1) / 2);
freq = (0:k) / sampTime / size(noiseFreeCurrent, temporalDimIndex);
LorentzCoeffs = [0.16  55  4;
                 0.045 190 2.5];
noiseSPD = lorentzSum(LorentzCoeffs, freq);

% Make-up the negative frequency part
noiseSPD = [noiseSPD noiseSPD(end:-1:1)];
noiseSPD = noiseSPD(1:size(noiseFreeCurrent, temporalDimIndex));
if (temporalDimIndex == 3)
    noiseSPD = reshape(noiseSPD, [1 1 length(noiseSPD)]);
elseif (temporalDimIndex == 2)
    noiseSPD = reshape(noiseSPD, [1 length(noiseSPD)]);
end

% Generate white Gaussian noise
noise = randn(size(noiseFreeCurrent));
noiseFFT = fft(noise, [], temporalDimIndex); % / sqrt(size(noise, 3));

% Adjust the spectral power distribution of the noise
noiseFFT = bsxfun(@times, noiseFFT, sqrt(noiseSPD));
% take the real part of noise
noise = real(ifft(noiseFFT, [], temporalDimIndex)) / sqrt(2 * sampTime);

% Add to noise-free signal
noisyCurrent = noiseFreeCurrent + noise;

end
