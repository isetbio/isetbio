%%t_osNoise  Illustrate the noise spectrum from the outer segment biophysical model
%
% Description:
%    
%    Illustrate the noise spectrum from the outer segment biophysical model we
%    use for converting photon absorptions to photo current.
%
%    Notes from FR on obtaining the noise frequency response amplitude.
%    The formula used is for the noise frequency response amplitude.
%    The 2 is there for some negative frequency, and the correction for the
%    number of time samples is the usual FFT correction.
%    [DHB NOTE: I don't understand this, really.  If possible, understand and
%    expand the comment.]
%
% See also: osAddNoise.

% BW ISETBIO Team, 2016
%
% 08/07/17  dhb  Commenting pass.

%%  Show that the mean is near zero
%
% We are putting in absorption rate of zero for a 10 x 10 array of cones,
% over the nSamp sample times. 
% The timing is 0.2 ms (1/5000)
nSamp = 10000; deltaT = 1/5000;
[noise, freq] = osAddNoise(zeros(10,10,nSamp), 'sampTime',deltaT);
fprintf('Mean noise should be close to zero: %f (pA)\n',mean(noise(:)));

%% Calculate the spectrum
noiseF = squeeze(mean(mean(abs(fft(noise, [], 3)*(2/nSamp)).^2)));

% Plot
vcNewGraphWin; loglog(freq, noiseF(1:length(freq)));
line('Xdata',[0.1 freq(end)],'Ydata',[.205 .205]);
xlabel('Frequency(Hz)'); ylabel('Power Spectrum (pA^2/Hz)');
grid on


