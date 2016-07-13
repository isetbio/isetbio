%
% PAL_SDT_PS_uneqDPtoPCpartFunc2(X,dP,i,M,Q,n)
% 
% Internal function
%
% Introduced: Palamedes version 1.7.0 (FK&NP)

function Y=PAL_SDT_PS_uneqDPtoPCpartFunc2(X,dP,M,Q,n)

NPD=PAL_pdfNormal(X,0,1);
P1=PAL_ZtoP(X);
P1=P1.^(Q*M-n-1);

P2=1.0;
for j=1:n
    P2=P2.*PAL_ZtoP(X-dP(j));
end

Y=(Q-n).*NPD.*P1.*P2;