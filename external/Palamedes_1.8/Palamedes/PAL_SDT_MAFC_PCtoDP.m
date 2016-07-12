%
% PAL_SDT_MAFC_PCtoDP converts proportion correct into d'(d-prime) for a 
% standard M-AFC (M-alternative-forced-choice) task, assuming an unbiased 
% observer
%
% Syntax: [dP]=PAL_SDT_MAFC_PCtoDP(pC,M);
% 
% returns a scalar, vector or matrix of d' ('dP') for a scalar, vector or
% matrix of proportion correct ('pC') and a scalar value of M (number of 
% alternatives), defined in the ranges 0<p<1 (p=proportion) and 2<M<inf.
% Note that this routine performs an iterative search on the inverse 
% routine PAL_SDT_MAFC_DPtoPC, which uses the Matlab function quadgk as 
% the default to perform a numerical integration between the bounds -Inf 
% and Inf.  For versions of Matlab that do not have quadgk, quadl is used 
% instead with the bounds set to -12 and 12.  The results using quadl are
% accurate except in extreme cases where d' is very high (>8)
%
% Example:
%
% dP=PAL_SDT_MAFC_PCtoDP([.1 .3 .5 .7 .9],10)
% 
% dP =
% 
%    -0.0000    0.8680    1.4748    2.0872    2.9829
%
% The example input is an N=5 vector of proportion correct with M set to 
% 10, and the output is the corresponding N=5 d's
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.4.0, 1.6.3 (see History.m)

function dP=PAL_SDT_MAFC_PCtoDP(pC,M)

[rows, cols]=size(pC);

func=@PAL_SDT_MAFC_DPtoPC;

dP = zeros(rows,cols);

for r=1:rows
    for c=1:cols
        dP(r,c)=PAL_minimize(@PAL_sqDistanceYfuncX,1,[],pC(r,c),func,M);
        if pC(r,c)==1.0
            dP(r,c)=1/0;
        end
    end
end