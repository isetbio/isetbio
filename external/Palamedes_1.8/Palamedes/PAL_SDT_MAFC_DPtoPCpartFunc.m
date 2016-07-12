%
% PAL_SDT_MAFC_DPtoPCpartFunc
%   
%   Internal function
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.6.0 (NP)

function Y=PAL_SDT_MAFC_DPtoPCpartFunc(X,dP,M)

N=PAL_pdfNormal(X-dP,0,1);
P=PAL_ZtoP(X);
Y=N.*P.^(M-1);