%
% PAL_SDT_3AFCoddity_IndMod_DPtoPC converts d'(d-prime) to proportion 
% correct for a 3-AFC (3-alternative-forced-choice) task, 
% assuming an independent observation model and an unbiased observer. 
% 
% Uses the method described in  Versfeld, Dai & Green (1996) Perception & 
% Psychophysics, 58, 10-51.
%
% Syntax: [pC]=PAL_SDT_3AFCoddity_IndMod_DPtoPC(dP);
%
% returns a scalar, vector or matrix of proportion correct ('pC') for a 
% scalar, vector or matrix of d' ('dP'), defined in the ranges 0<d'<inf.  
% 
% Note that the Matlab function quadgk, which performs a numerical 
% integration between the bounds -Inf and Inf, is used as the default in 
% the routine. For versions of Matlab that do not have quadgk, quadl is 
% used instead with the bounds set to -12 and 12.  The results using quadl 
% are accurate except in extreme cases where d' is very high (>8)
%
% Note also that unlike PAL_SDT_MFCoddity_IndMod_DPtoPC, which is the 
% M-AFC version of the same task but which uses Monte Carlo simulation, 
% this routine is deterministic and hence much faster
%
%
% Example:
%
% [pC] = PAL_SDT_3AFCoddity_IndMod_DPtoPC([0 1 2 3 4 5])
%
% returns:
% 
% pC =
%
%    0.3333    0.4468    0.6774    0.8611    0.9534    0.9875
%
% The example input consists of a N=6 vector of d', and the output is the 
% corresponding N=6 proportion correct
%
% Introduced: Palamedes version 1.6.0 (FK & NP)
% Modified: Palamedes version 1.6.3(see History.m)

function pC = PAL_SDT_3AFCoddity_IndMod_DPtoPC(dP)

[rows, cols] = size(dP);

pC = zeros(rows, cols);

if exist('quadgk.m','file') == 2
    
    for r = 1:rows
        for c = 1:cols
            funcA = PAL_ZtoP(dP(r,c)./2).^3;
            integralB = quadgk(@(X)PAL_SDT_3AFCoddity_IndMod_DPtoPCpartFuncA(X,dP(r,c)),-Inf,-dP(r,c)./2);
            funcC = (1-PAL_ZtoP(dP(r,c)./2)).^3;
            integralD = quadgk(@(X)PAL_SDT_3AFCoddity_IndMod_DPtoPCpartFuncB(X,dP(r,c)),-dP(r,c)./2,Inf);
            pC(r,c) = funcA+integralB+funcC+integralD;
        end
    end
    

else
    
    if (max(dP)>8)
        warning('d-prime values >8 may not yield accurate results');
    end
    
    for r = 1:rows
        for c = 1:cols
            funcA = PAL_ZtoP(dP(r,c)./2).^3;
            integralB = quadl(@(X)PAL_SDT_3AFCoddity_IndMod_DPtoPCpartFuncA(X,dP(r,c)),-12,-dP(r,c)./2);
            funcC = (1-PAL_ZtoP(dP(r,c)./2)).^3;
            integralD = quadl(@(X)PAL_SDT_3AFCoddity_IndMod_DPtoPCpartFuncB(X,dP(r,c)),-dP(r,c)./2,12);
            pC(r,c) = funcA+integralB+funcC+integralD;
        end
    end
    
end