%
%PAL_SDT_ROCML_negLL     (negative) Log Likelihood associated with fit 
%   of ROC curve
%
%
%Internal Function
%
%Introduced: Palamedes version 1.6.0 (FK & NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function negLL = PAL_SDT_ROCML_negLL(paramsFreeVals, paramsFixedVals, paramsFree, cumNumHF, OutOfNum, invSDTF)

paramsValues(paramsFree == 1) = paramsFreeVals;
paramsValues(paramsFree == 0) = paramsFixedVals;

dP = paramsValues(1);
R = paramsValues(2);
C = paramsValues(3:length(paramsValues));

DPvec = repmat(dP,1,length(C)); % converts a scalar dP into a vector of identical values
Rvec = repmat(R,1,length(C));

pHF = invSDTF(DPvec,C,'ratioSD',Rvec); %takes three vectors of equal length as input and returns a matrix

negLL = -sum(PAL_nansum(cumNumHF.*log(pHF))+PAL_nansum((OutOfNum-cumNumHF).*log(1 - pHF)));