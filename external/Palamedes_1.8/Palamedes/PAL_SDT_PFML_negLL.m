%
%PAL_SDT_PFML_negLL     (negative) Log Likelihood associated with fit 
%   of SDT (signal-detection-theory) psychometric function 
%
%Internal Function
%
%Introduced: Palamedes version 1.8.0 (FK & NP)

function negLL = PAL_SDT_PFML_negLL(params,StimLevels,NumPos,OutOfNum,SDTfunc,M)

DP=(params(1).*StimLevels).^params(2);

if isempty(M) 
    PC = SDTfunc(DP); 
else
    PC = SDTfunc(DP,M);
end

negLL = -sum(PAL_nansum(NumPos.*log(PC))+PAL_nansum((OutOfNum-NumPos).*log(1 - PC)));