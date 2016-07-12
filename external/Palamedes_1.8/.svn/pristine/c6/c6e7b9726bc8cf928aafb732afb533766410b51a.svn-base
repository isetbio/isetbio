%
% PAL_SDT_2AFC_DPtoPHF converts d'(d-prime) for a criterion C into 
% proportion hits and proportion false alarms for a 2AFC 
% (two-alternative-forced-choice) task
%
% Syntax: [pHF]=PAL_SDT_2AFC_DPtoPHF(dP,C);
% 
% returns a Nx2 matrix of N proportion hits and proportion false alarms 
% ('pHF') for a scalar or N-length vector of d' ('dP') and criterion C 
% ('C'), defined in the ranges 0<d'<inf and -inf<C<inf
%
% Example:
%
% [pHF]=PAL_SDT_2AFC_DPtoPHF([.0 1 10],[-2.5 0 4.0])
%
% returns:
% 
% pHF =
% 
%     0.9615    0.9615
%     0.7602    0.2398
%     1.0000         0
%
% The example input consists of two N=3 vectors (each demarcated by 
% brackets) of d' and C. The first column of the 3 x 2 output matrix 
% gives the proportion of hits and the second column the corresponding 
% proportion of false alarms
%
% Introduced: Palamedes version 1.0.0 (FK)

function pHF = PAL_SDT_2AFC_DPtoPHF(dP,C)

zH=(dP-C)./sqrt(2);
zF=-(dP+C)./sqrt(2);

pHF(:,1)=PAL_ZtoP(zH);
pHF(:,2)=PAL_ZtoP(zF);