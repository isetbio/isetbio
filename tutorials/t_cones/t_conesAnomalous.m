%% t_conesAnomalous
%
% How to shift the photopigment absorbance
%
% Built on PsychToolbox and Asano/Lamb pigment shifts


%%
cm = cMosaic;

lAbsorbance = cm.pigment.absorbance(:,1);
wave = cm.wave;

ieNewGraphWin;
plot(wave,lAbsorbance);

% We will shift the Stockman fundamentals.
% [c,w] = ieReadSpectra('stockmanQuanta');

% PTB format
S = WlsToS(wave(:));

absorbance = ShiftPhotopigmentAbsorbance(S,lAbsorbance',10,'log');
ieNewGraphWin; plot(wave,lAbsorbance(:),'--',wave,absorbance,'go');

