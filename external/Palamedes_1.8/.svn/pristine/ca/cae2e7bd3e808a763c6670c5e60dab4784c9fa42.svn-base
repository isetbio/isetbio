%
% PAL_SDT_1AFC_PHFtoDP_Demo is a script that demonstrates how to output
% in more user-friendly format the results of the routine 
% PAL_SDT_1AFC_PHFtoDP, which converts proportion of hits and proportion 
% of false alarms into d' (d-prime), criterion C, criterion lnBeta and
% proportion correct for a 1AFC (one-alternative-forced-
% choice task, such as Yes/No or a symmetric single-interval task
%
% Syntax: PAL_SDT_1AFC_PHFtoDP_Demo
%
% asks for a Nx2 matrix or N proportion hits and proportion false
% alarms and returns vectors of proportion hits (pHit), proportion false 
% alarms (pFA) d' (d-prime), proportion correct (p Corr), criterion C
% (crit C) and criterion ln Beta (crit lnB)
%
% Example:
%
% PAL_SDT_1AFC_PHFtoDP_Demo
% Enter a Nx2 matrix of N proportion Hits and N False Alarms [0.6 0.1; 0.8 0.3; 0.9 0.4]
%
% returns:
% 
%      pHit      pFA      d-prime   p Corr    crit C    crit lnB
%     0.6000    0.1000    1.5349    0.7500    0.5141    0.7891
%     0.8000    0.3000    1.3660    0.7500   -0.1586   -0.2167
%     0.9000    0.4000    1.5349    0.7500   -0.5141   -0.7891
%
% The example input argument consists of a 3 x 2 matrix, with each
% row (demarcated by a semi-colon) consisting of a proportion of hits and 
% a corresponding proportion of false alarms 
%
%FK (September 2009)

clear all;

pHF=input('Enter a Nx2 matrix of N proportion Hits and N False Alarms ');

[dP C lnB pC]=PAL_SDT_1AFC_PHFtoDP(pHF);

pH=pHF(:,1);
pF=pHF(:,2);

SDM=[pH';pF';dP';pC';C';lnB'];
SDM=SDM';

fprintf('\n');
disp('     pHit      pFA      d-prime   p Corr    crit C    crit lnB');
disp(SDM);