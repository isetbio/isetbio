%
% PAL_SDT_1AFCsameDiff_DiffMod_PHFtoDP converts proportion hits and 
% proportion false alarms into d'(d-prime) and criterion k for a 1AFC 
% (one-alternative-forced-choice) same-different task under the
% Differencing Observer model
%
% Syntax: [dP k]=PAL_SDT_1AFCsameDiff_DiffMod_PHFtoDP(pHF);
% 
% returns a scalar or N-length vector of d' ('dP') and criterion k ('k'), 
% for an Nx2 input matrix of N proportion hits and proportion false 
% alarms ('pHF') defined in the range 0<p<1 (p=proportion)
%
% Example:
%
% [dP k]=PAL_SDT_1AFCsameDiff_DiffMod_PHFtoDP([.6 .1; .7 .1; .8 .5])
%
% returns:
% 
% dP =
% 
%     2.6837
%     3.0675
%     2.0629
% 
% 
% k =
% 
%     2.3262
%     2.3262
%     0.9539
%
% The example input argument is a 3 x 2 matrix in which each row 
% (demarcated by a semi-colon) consists of a proportion of hits and a 
% corresponding proportion of false alarms.  The columns in the output 
% are the resulting N=3 vectors of dP and k
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.4.0, 1.6.3 (see History.m)


function [dP, k]=PAL_SDT_1AFCsameDiff_DiffMod_PHFtoDP(pHF)

[rows, cols]=size(pHF);

k=-sqrt(2).*PAL_PtoZ(pHF(:,2)./2);
pH=pHF(:,1);

func=@PAL_SDT_1AFCsameDiff_DiffMod_DPtoPH;

dP = zeros(rows,1);

for r=1:rows
    dP(r) = PAL_minimize(@PAL_sqDistanceYfuncX,1,[],pH(r),func,k(r));
end