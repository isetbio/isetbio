function [OTF2D, frequencySupport, wave] = mouseOTF(pupilRadius, dioptricPower, frequencySupport, wave)
% Calculate the mouse OTF, including chromatic aberration
% This very similar to the human OTF, but with different defocus
% parameters, and no Williams factor.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Comments for humanOTF
%  [OTF2D, frequencySupport, wave] = ...
%        humanOTF([pupilRadius = 0.0015], [dioptricPower = 59.9404], [frequencySupport = :],[wave = :])
%
% The spatial frequency range is determined by the spatial extent and
% sampling density of the original scene.
% 
% Reference: Marimont & Wandell (1994 --  J. Opt. Soc. Amer. A,  v. 11, p.
% 3113-3122 -- see also Foundations of Vision by Wandell, 1995.
%
% See Also:  humanLSF, sceneGet(scene,'frequencyresolution')
%
% Example:
%   [OTF2D, frequencySupport, wave] = humanOTF(0.0015,60);
%   vcNewGraphWin;
%   subplot(1,2,1), 
%   mesh(frequencySupport(:,:,1),frequencySupport(:,:,2),abs(OTF2D(:,:,10))); 
%   title('500 nm'); xlabel('Frequency (cyc/deg)'), zlabel('Relative amp')
%   subplot(1,2,2), mesh(frequencySupport(:,:,1),frequencySupport(:,:,2), abs(OTF2D(:,:,3))); 
%   set(gca,'zlim',[-.2,1]);
%   xlabel('Frequency (cyc/deg)'), zlabel('Relative amp'); title('400 nm')
%
% Reference and discussion
%
% We build the otf by first using Hopkins' formula of an eye with only
% defocus and chromatic aberration.  Then, we multiply in an estimate of
% the other aberrations.  At present, we are using some data from Dave
% Williams and colleagues measured using double-pass and threshold data.  
%
% Williams et al. (19XX) predict the measured MTF at the infocus wavelength by 
% multiplying the diffraction limited OTF by a weighted exponential.  
% We perform the analogous calculation at every wavelength.  That is, we
% multiply the aberration-free MTF at each wavelength by the weighted
% exponential in the Williams measurements.  Speaking with Dave last month,
% he said his current experimental observations confirmed that this was
% an appropriate correction.  (BW 05.24.96).
%
% As a further simplification, the human measurements are all 1D.  We build
% a 1D function and then we assume that the true function is circularly
% symmetric.  That is how we fill in the full 2D MTF.  We call it an OTF
% and assume there is no phase shift ...  All an approximation, but
% probably OK for these types of calculations.  Further details could be
% sought out in the recent papers from the Spanish group (e.g. Artal) and
% from Williams and the other Hartmann Shack people.
%
% Copyright ImagEval Consultants, LLC, 2003.
%%%%%%%%%%%%%%%%%%%%%%%%% end comments for humanOTF
% For the mouse, we'll use different defocus values, and no williams
% factor.
% The values come from : A schematic eye for the mouse, and comparisons
% with the rat, Remtulla and Hallett, 1984.
%
% Note : The article talks about ametropia, whereas the code talks about
% defocus. Both are a problem of focus out of the focal plane, and both are
% in diopters, so I'm guessing they're the same thing.


% Default mouse pupil radius is 0.00059 m = 0.59 mm.  
if notDefined('pupilRadius'), p = 0.00059; else p = pupilRadius; end 

% dioptric power of unaccomodated mouse eye (0.001756m = 1.756 mm focal length, at 544nm)
if notDefined('dioptricPower')
    D0 = 1/0.001756; 
else
    D0 = dioptricPower; 
end

% Wavelength in nanometers
if notDefined('wave'), wave = (400:700)'; end
nWave = length(wave);

% We use a frequency support that covers 2 cyc/deg.
% (the human code uses 60 cycles/deg, but that's way too generous for a mouse. 
% The peak spatial frequency sensitivity is around 0.2 cyc/deg, and the 
% cutoff frequency at less than 1 cyc/deg.)
maxFrequency = 3; % cyc/deg
%maxFrequency = 20; % cyc/deg
if notDefined('frequencysupport'), 
  %  fList = unitFrequencyList(maxFrequency);
    fList = -1:.05:1;
    fList = fList*maxFrequency;
    [X,Y] = meshgrid(fList,fList);
    frequencySupport(:,:,1) = X; frequencySupport(:,:,2) = Y;
end

% We treat the OTF as a circularly symmetric function.  We treat the
% effective frequency as the distance from the origin. 
dist = sqrt((frequencySupport(:,:,1).^2 + frequencySupport(:,:,2).^2));
t = max(frequencySupport(:,:,1)); maxF1 = max(t(:));
t = max(frequencySupport(:,:,2)); maxF2 = max(t(:));

%  We don't want to allow any output frequencies beyond the circle defined
%  by the minimum of the two largest frequencies.  We will zero those terms
%  later.
maxDist = min(maxF1,maxF2);
% mesh(dist)

% The human OTF is smooth (we're guessing the mouse's is too).
% To speed up calculations, we use 40 samples and interpolate the other
% values. The sample spatial frequencies are in cycles per degree.
sampleSF = ((0:39)/39)*max(dist(:));
otf      = mouseCore(wave,sampleSF,p,D0);
% WARNING : using [0;60] cyc/deg makes
% complex values apear at the high frequency, and the abs(otf) has huge
% peaks at these values. Maybe they are due to rounding errors. In any
% case, don't go this high in frequency. The mouse's frequency cutoff is
% lower than 1 cyc/deg anyway.
% NOTE : There are a bunch of negative values in the undulations on the
% sides. Does this imply contrast inversion?

% plot the otf :
% figure; [x,y]=meshgrid(wave,sampleSF); 
% surf(x',y',otf);
% xlabel('wavelength'); ylabel('cycles/degree')
% title('Mouse OTF')
% plot a small part of the spacial frequencies only
%figure; [x,y] = meshgrid(wave, sampleSF(1:6));
%surf(x',y',abs(otf(:,1:6)));
%xlabel('wavelength'); ylabel('cycles/degree')


% Interpolate the full 2D OTF from the individual values.
[r,c] = size(frequencySupport(:,:,1)); 
OTF2D = zeros(r,c,nWave);
l = (dist > maxDist);
%wBar = waitbar(0,'humanOTF');
for ii=1:nWave
%    waitbar(ii/nWave,wBar);
    % We have small imaginary values sometimes.  Probably rounding error
    % in some calculation above.  We remove them here.
    tmp = abs(interp1(sampleSF,otf(ii,:),dist,'spline'));
    
    % We don't want any frequencies beyond the sampling grid.  Here we
    % zero them out.
    tmp(l) = 0;
    
    % This is the proper storage format for the OI-ShiftInvariant case.
    OTF2D(:,:,ii) = fftshift(tmp);
end
%close(wBar);

% plot the OTF2D
end
