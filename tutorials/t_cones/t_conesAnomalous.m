%% t_conesAnomalous
%
% How to shift the photopigment absorbance
%
% Built on PsychToolbox and Asano/Lamb pigment shifts
%
% Notes:
%   The photopigment properties are derived from the internal
%   'absorbance_' and 'wave_' variables.  When the public 'wave'
%   variable is changed, the available 'absorbance' variable is
%   re-interpolated.
%
%   I am not sure this is a great approach.  It could be that we
%   should have 'absorbance_' be specified by a file on disk that is
%   re-interpolated when 'wave' is changed.
% 
% 

%%
cm = cMosaic;

lAbsorbance = cm.pigment.absorbance_(:,1);
wave = cm.pigment.wave_;

cm.plot('spectral qe');

% Shift the red to the right
absorbance = ShiftPhotopigmentAbsorbance(wave(:),lAbsorbance',-20,'log');

% ieNewGraphWin; plot(wave,lAbsorbance(:),'--',wave,absorbance,'go');
cm.pigment.absorbance_(:,1) = absorbance;

cm.plot('spectral qe');

%% END