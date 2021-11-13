%% t_conesAnomalous
%
% How to shift the photopigment absorbance
%
% Built on PsychToolbox and Asano/Lamb pigment shifts


%%
cm = cMosaic;

lAbsorbance = cm.pigment.absorbance(:,1);
wave = cm.wave;
cm.plot('spectral qe');

% Shift the red to the right
absorbance = ShiftPhotopigmentAbsorbance(wave(:),lAbsorbance',10,'log');

% ieNewGraphWin; plot(wave,lAbsorbance(:),'--',wave,absorbance,'go');
cm.pigment.absorbance(:,1) = absorbance;
cm.plot('spectral qe');


%% END