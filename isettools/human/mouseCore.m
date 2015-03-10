function [otf, williamsFactor] = mouseCore(wave,sampleSF,p,D0)
% Compute the mouse optical transfer function. This is closely copying the
% human OTF in humanCore, with mouse parameters.
% The Williams factor is 1 : we don't use it.
%
%%%%%%%%%%%%%%%%%%%%%%%%%% comments for humanCore
%   [otf, williamsFactor] = humanCore(wave,nWave,sampleSF,p,D0)
%
% The computation of the human OTF is performed here.  It is placed in this
% routine to permit compilation.  The OTF is computed at the wavelengths in
% wave and at sampleSF spatial frequencies in units of cycles/deg
%
% The returned williamsFactor is a general loss of contrast, independent of
% wavelength, that we apply to bring the data into alignment with human
% optics.  That factor comes from empirical measurements out of the
% Williams' lab.
%
%  p:         Pupil radius in meters
%  D0:        Dioptric power (accomodation), usually around 60
%  sampleSF:  Spatial frequencies in cycles/deg
%  wave:      Wavelength in nanometers
%
% Example:
%   wave = 400:10:700; sampleSF = [0:1:30]; 
%   p = 0.0015; D0 = 60;
%   otf = humanCore(wave,sampleSF,p,D0);
%   mesh(sampleSF,wave,otf)
%
% Copyright ImagEval Consultants, LLC, 2005.
%%%%%%%%%%%%%%%%%%%%%%%%%%%% end comments from humanCore

%% Defocus
% The defocus values come from : A schematic eye for the mouse, and
% comparisons with the rat, Remtulla and Hallett, 1984.
% The article gives mesures of defocus (ametropia, which I think is the
% same thing) for 4 wavelength. We extend the curve to get the full
% spectrum. 
D488 = -9.4;
D544 = 0.4;
D598 = 6.6;
D655 = 10.8;

% Interpolate with bits of straight lines
baseWave = 488:655;
Dmouse = zeros(size(baseWave));

w1 = 488; w2 = 544; D1 = D488; D2 = D544; 
a = (D1-D2)/(w1-w2);
b = (w1*D2 - D1*w2)/(w1-w2);
Dmouse(1:57) = a*baseWave(1:57) + b;
aSmall = a; bSmall = b;

w1 = 544; w2 = 598; D1 = D544; D2 = D598;
a = (D1-D2)/(w1-w2);
b = (w1*D2 - D1*w2)/(w1-w2);
Dmouse(57:111) = a*baseWave(57:111) + b;

w1 = 598; w2 = 655; D1 = D598; D2 = D655;
a = (D1-D2)/(w1-w2);
b = (w1*D2 - D1*w2)/(w1-w2);
Dmouse(111:end) = a*baseWave(111:end) + b;
aBig = a; bBig = b;

% Extend with same slope to fill full scope of wavelengths
smallWave = wave(1):487;
bigWave = 656:wave(end);
Dmouse = [aSmall*smallWave + bSmall, Dmouse, aBig*bigWave + bBig];
baseWave = [smallWave, baseWave, bigWave];

% plot
%figure; hold on; 
%plot(baseWave, Dmouse); plot([488, 544, 598, 655],[D488, D544, D598, D655], 'o');
%xlabel('wavelength in nm'); ylabel('Defocus in diopters');
%q1 = 1.7312; q2 = 0.63346; q3 = 0.21410;
%Dhuman = q1-q2./(baseWave*1e-3-q3);     
%plot(baseWave,Dhuman, 'r');
%legend('mouse','','human')
%title('Defocus for the mouse and human eye')

% Take values at wave points
idx = floor(wave-wave(1) + 1);
D = Dmouse(idx);

%% Compute the OTF
% We use Hopkins' model for a diffraction-limited lens with defocus and
% chromatic aberration. 
% See "Matching color images: the effects of axial chromatic aberration", 
% Marimont and Wandell, 1994.

% w20 = defocus with respect to optical path length error.
% D0 is dioptric power of the unnacommodated eye, at in-focus wavelength (approx 544nm for the mouse.)
w20 = p^2/2*(D0.*D)./(D0+D);

% We skip the Williams factor. It is typical of the human eye, and so not 
% useful for the mouse eye. (TODO : what is it?)
% % There is a typical human OTF scaling we use from the work at Dave
% % Williams' lab.  Here is a smooth fit to their data.  This was provided by
% % Dave Brainard and could be updated or drawn from the literature in some
% % other way.  Perhaps from Ijspeert?
% a =  0.1212;		%Parameters of the fit
% w1 = 0.3481;		%Exponential term weights
% w2 = 0.6519;
% williamsFactor =  w1*ones(size(sampleSF)) + w2*exp( - a*sampleSF );
williamsFactor = 1;

% We use this factor to convert from the input spatial frequency units
% (cycles/deg) to cycles/meter needed for the Hopkins eye
fLength = 0.001756;  % mouse focal length : 1.756 mm. See discussion in opticsCreate.
metersPerDegree = fLength*tan(1/180*pi);
degreesPerMeter = 1/metersPerDegree;
c = degreesPerMeter;            % degrees per meter for mouse eye

% Variables for Hopkins' formula
% s is the reduced spacial frequency.
s = zeros(length(wave), length(sampleSF));
% alpha is the variable called "a" in the formula (no physical significance?)
alpha = zeros(length(wave), length(sampleSF));
% otf, computed in opticsDefocusedMTF
otf = zeros(length(wave), length(sampleSF));
for ii = 1:length(wave)
    
    s(ii,:) = (c * wave(ii)*1e-9 /(D0*p)) * sampleSF;
    
    alpha(ii,:) = 4*pi./(wave(ii)*1e-9 ).*w20(ii).*s(ii,:);
    
    % We put the sample SF into this array.
    % Then we interpolate to the full 2D array outside of this loop.
    otf(ii,:) = opticsDefocusedMTF(s(ii,:),abs(alpha(ii,:)));
    
    otf(ii,:) = otf(ii,:).*williamsFactor; % williamsFactor is 1 for the mouse.
    
end

end