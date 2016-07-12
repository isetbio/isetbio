%
% PAL_SDT_PS_DPtoPCpartFunc(X,dP,M,Q,n)
% 
% Internal function
%
% Introduced: Palamedes version 1.7.0 (FK&NP)

function Y=PAL_SDT_PS_DPtoPCpartFunc(X,dP,M,Q,n)

NPD=PAL_pdfNormal(X-dP,0,1);
P1=PAL_ZtoP(X);
P1=P1.^(Q*M-n);
P2=PAL_ZtoP(X-dP);
P2=P2.^(n-1);
Y1=n.*NPD.*P1.*P2;

NPD=PAL_pdfNormal(X,0,1);
P1=PAL_ZtoP(X);
P1=P1.^(Q*M-n-1);
P2=PAL_ZtoP(X-dP);
P2=P2.^n;
Y2=(Q-n).*NPD.*P1.*P2;

Y=(Y1+Y2);