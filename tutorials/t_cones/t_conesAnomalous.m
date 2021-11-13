%% t_conesAnomalous
%
% How to shift the photopigment absorbance
%
% Built on PsychToolbox and Asano/Lamb pigment shifts

% The photopigment properties are derived from the absorbance_, and density
% terms   

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