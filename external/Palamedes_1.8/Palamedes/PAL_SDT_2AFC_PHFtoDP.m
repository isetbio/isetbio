%
% PAL_SDT_2AFC_PHFtoDP converts proportion hits and proportion false 
% alarms into d'(d-prime) and criterion C for a 2AFC 
% (two-alternative-forced-choice) task
%
% Syntax: [dP C lnB pC]=PAL_SDT_2AFC_PHFtoDP(pHF);
% 
% returns a scalar or N-length vector of d' ('dP'), criterion C ('C'), 
% criterion lnBeta ('lnB') and proportion correct ('pC'), for an Nx2 input 
% matrix of N proportion hits and proportion false alarms ('pHF') defined 
% in the range 0<p<1 (p=proportion) 
%
% Example:
%
% [dP C lnB pC]=PAL_SDT_2AFC_PHFtoDP([.5 .3; .7 .2; .9 .1])
%
% returns:
% 
% dP =
% 
%     0.3708
%     0.9659
%     1.8124
% 
% 
% C =
% 
%     0.3708
%     0.2243
%          0
% 
% 
% lnB =
% 
%     0.1375
%     0.2167
%          0
% 
% 
% pC =
% 
%     0.6000
%     0.7500
%     0.9000
%
% In the example input, each row of the 3 x 2 matrix (demarcated by 
% a semi-colon) consists of a proportion of hits and a corresponding 
% proportion of false alarms.   The columns in the output give the  
% resulting N=3 vectors of dP, C, lnB and pC
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.6.3 (see History.m)

function [dP, C, lnB, pC]=PAL_SDT_2AFC_PHFtoDP(pHF)

zH=PAL_PtoZ(pHF(:,1));
zF=PAL_PtoZ(pHF(:,2));

dP=(zH-zF)./sqrt(2);
C=-(zH+zF)./sqrt(2);
lnB=(zF.^2-zH.^2)./2;
pC=(pHF(:,1)+(1.0-pHF(:,2)))./2;