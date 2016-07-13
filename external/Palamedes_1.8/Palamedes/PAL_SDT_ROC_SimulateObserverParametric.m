%
%PAL_SDT_ROC_SimulateObserverParametric  Simulate observer characterized by
%   an ROC curve.
%
%syntax: cumNumHF = PAL_SDT_ROC_SimulateObserverParametric(paramsValues, ...
%   OutOfNum, invSDTF)
%
%Input:
%
%   'paramsValues': vector consisting of a d-prime value, a ratio-of-SDs
%       value and N values of criterion C.  These may be estimates
%       derived from PAL_SDT_ROCML_Fit
%
%   'OutOfNum' is a Nx2 matrix of number of trials corresponding to each
%       value in cumNumHF.  PAL_SDT_cumulateHF returns OutOfNum as well as
%       cumNumHF
%
%   'invSDTF' is the SDT function that converts a proportion of 
%       Hits and False Alarms into a d-prime and criterion C value, 
%       for a given SD ratio of signal-to-noise.  For example:
%
%       invSDTF = @PAL_SDT_1AFC_DPtoPHF;
%
%Output:
%
%   'cumNumHF': Nx2 matrix of cumulative number of Hits and False Alarms
%
%Example:
%  
%   paramsValues = [1 1 -1:.5:1];
%   OutOfNum = ones(5,2)*70;
%   invSDTF = @PAL_SDT_1AFC_DPtoPHF;
%
%   %Simulate observer:
%
%   cumNumHF = PAL_SDT_ROC_SimulateObserverParametric(paramsValues,...
%       OutOfNum, invSDTF)
%
%returns something like:
%
%cumNumHF =
%
%    64    52
%    62    37
%    43    17
%    31    17
%    25     1
%
%Introduced: Palamedes version 1.6.0 (FK & NP)

function cumNumHF = PAL_SDT_ROC_SimulateObserverParametric(paramsValues, OutOfNum, invSDTF)

dP = paramsValues(1);
R = paramsValues(2);
C = paramsValues(3:length(paramsValues));

DPvec = repmat(dP,1,length(C)); % converts a scalar d-prime into a vector with identical values
Rvec = repmat(R,1,length(C));

[pHF] = invSDTF(DPvec,C,'ratioSD',Rvec);

pHF3 = repmat(pHF, [1 1 OutOfNum(1,1)]);

Pos = rand(size(pHF,1),size(pHF,2), OutOfNum(1,1));
Pos(Pos < pHF3) = 1;
Pos(Pos ~= 1) = 0;
cumNumHF = sum(Pos,3);