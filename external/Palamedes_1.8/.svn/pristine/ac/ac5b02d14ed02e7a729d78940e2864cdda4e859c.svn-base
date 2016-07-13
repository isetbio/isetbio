%
% PAL_SDT_1AFCsameDiff_DiffMod_DPtoPH
%   
%   Internal function: note this is a different function from 
%   the user function PAL_SDT_1AFCsameDiff_DiffMod_DPtoPHF which
%   outputs proportion hits AND proportion false alarms
%
% Introduced: Palamedes version 1.0.0 (FK)

function pH=PAL_SDT_1AFCsameDiff_DiffMod_DPtoPH(dP,k)

pH=PAL_ZtoP((-k+dP)./sqrt(2)) + PAL_ZtoP((-k-dP)./sqrt(2));