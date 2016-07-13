%
% PAL_SDT_1AFCsameDiff_DiffMod_DPtoPHF converts d' (d-prime) with 
% criterion k into proportion hits and proportion false alarms for a 1AFC 
% (one-alternative-forced-choice) same-different task under the 
% Differencing Observer model
%
% Syntax: [pHF]=PAL_SDT_1AFCsameDiff_DiffMod_DPtoPHF(dP,k);
%
% returns a Nx2 matrix of N proportion hits and proportion false alarms 
% ('pHF') for a scalar or N-length vector of d' ('dP') and criterion k 
% ('k'), defined in the ranges 0<d'<inf and 0<k<inf
%
% Example:
%
% [pHF]=PAL_SDT_1AFCsameDiff_DiffMod_DPtoPHF([.5 1.5 3],[0 1 2])
%
% returns:
%
% pHF =
%
%    1.0000    1.0000
%    0.6767    0.4795
%    0.7605    0.1573
%
% The example input arguments are two N=3 vectors of d' and C.  The
% first column in the 3 x 2 matrix output gives the resulting proportion 
% of hits and the second column the corresponding proportion of false 
% alarms
%
% Introduced: Palamedes version 1.0.0 (FK)

function pHF=PAL_SDT_1AFCsameDiff_DiffMod_DPtoPHF(dP,k)

pH=PAL_ZtoP((-k+dP)./sqrt(2)) + PAL_ZtoP((-k-dP)./sqrt(2));
pF=2.*PAL_ZtoP(-k./sqrt(2));
pHF(:,1)=pH;
pHF(:,2)=pF;