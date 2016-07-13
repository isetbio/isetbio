%
% PAL_SDT_MAFC_DPtoPC converts d'(d-prime) to proportion correct for a 
% standard M-AFC (M-alternative-forced-choice) task, assuming an unbiased 
% observer
%
% Syntax: [pC]=PAL_SDT_MAFC_DPtoPC(dP,M);
% 
% returns a scalar, vector or matrix of proportion correct ('pC') for a 
% scalar, vector or matrix of d' ('dP') and a scalar value of M (number of 
% alternatives), defined in the ranges 0<d'<inf and 2<M<inf.  Note that 
% the Matlab function quadgk, which performs a numerical integration 
% between the bounds -Inf and Inf, is used as the default in the routine.  
% For versions of Matlab that do not have quadgk, quadl is used instead 
% with the bounds set to -12 and 12.  The results using quadl are
% accurate except in extreme cases where d' is very high (>8)
%
%
% Example:
%
% [pC]=PAL_SDT_MAFC_DPtoPC([0 1 2 3 4 5],25)
%
% returns:
% 
% pC =
% 
%     0.0400    0.1997    0.5218    0.8263    0.9648    0.9961
%
% The example input consists of a N=6 vector of d' with M set to 25, and 
% the output is the corresponding N=6 proportion correct
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.4.0, 1.6.3 (see History.m)

function pC=PAL_SDT_MAFC_DPtoPC(dP,M)

[rows, cols]=size(dP);

pC = zeros(rows,cols);

if exist('quadgk.m','file')==2
    
    limit = Inf;
    if exist('OCTAVE_VERSION','builtin')
        limit = 100;    %Octave does not like Inf as limit
    end
    
    for r = 1:rows
        for c = 1:cols
            pC(r,c)=quadgk(@(X)PAL_SDT_MAFC_DPtoPCpartFunc(X,dP(r,c),M),-limit,limit);
        end
    end
    

else
    
    if (max(dP)>8)
        warning('d-prime values >8 will not yield accurate results');
    end
    
    for r = 1:rows
        for c = 1:cols
            pC(r,c)=quadl(@(X)PAL_SDT_MAFC_DPtoPCpartFunc(X,dP(r,c),M),-12,12);
        end
    end
    
end