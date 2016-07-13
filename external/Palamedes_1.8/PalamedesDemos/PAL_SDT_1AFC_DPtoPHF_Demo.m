%
% PAL_SDT_1AFC_DPtoPHF_Demo is a script that demonstrates how to output
% in more user-friendly format the results of the routine 
% PAL_SDT_1AFC_DPtoPHF, which converts d' and criterion C into proportion
% hits and proportion false alarms for a 1AFC (one-alternative-frced-
% choice task, such as Yes/No or a symmetric single-interval task
%
% Syntax: PAL_SDT_1AFC_DPtoPHF_Demo
%
% asks for a scalar or vector of d' and criterion C, and outputs 
% scalars or vectors of d' (dprime), C (crit C), proportion of hits (pHit), 
% proportion of false alarms (pFA) and proportion correct (p Corr)
%
% Example:
%
% PAL_SDT_1AFC_DPtoPHF_Demo
% Enter a vector of Dprime values [0 1.5 3]
% Enter a vector of Criterion C values [1 1 1]
%
% returns:
% 
%     d-prime   crit C     pHit      pFA      p Corr
%          0    1.0000    0.1587    0.1587    0.5000
%     1.5000    1.0000    0.4013    0.0401    0.6806
%     3.0000    1.0000    0.6915    0.0062    0.8426
%
% The example inputs are N=3 vectors of d' and C.
%
%FK (September 2009)

clear all;

dP=input('Enter a vector of Dprime values ');
C=input('Enter a vector of Criterion C values ');

pHF = PAL_SDT_1AFC_DPtoPHF(dP,C);
pH=pHF(:,1);
pF=pHF(:,2);

pC=(pH+1.0-pF)./2;

pHFC=[dP;C;pH';pF';pC'];
pHFC=pHFC';

fprintf('\n');
disp('    d-prime   crit C     pHit      pFA      p Corr');
disp(pHFC);