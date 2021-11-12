%% t_conesAnomalous
%
% How to shift the photopigment absorbance
%
% Built on PsychToolbox and Asano/Lamb pigment shifts


%%
cm = cMosaic;

[c,w] = ieReadSpectra('stockmanQuanta');
S = WlsToS(w(:));
absorbance = ShiftPhotopigmentAbsorbance(S,c(:,1)',10,'log');
ieNewGraphWin; plot(w,c(:,1),'--',w,absorbance,'go');

