%
% PAL_SDT_MAFCmatchSample_DiffMod_PCtoDP converts proportion correct into 
% d'(d-prime) for a M-AFC (M-alternative-forced-choice) match-to-sample
% task, assuming a Differencing Observer model and an unbiased observer
%
% Syntax: [dP]=PAL_SDT_MAFCmatchSample_DiffMod_PCtoDP(pC,M);
% 
% returns a scalar, vector or matrix of d' ('dP') for an input scalar, 
% vector or matrix of proportion correct p ('pC') and an input value of M 
% (number of alternatives), defined in the ranges 0<p<1 and 2<M<inf. Note 
% that the routine may take several seconds or even minutes to
% execute depending on the size of the input and speed of computer, as
% the routine performs an iterative search on the inverse version of the 
% routine which employs Monte Carlo simulation
%
% Example:
%
% [dP]=PAL_SDT_MAFCmatchSample_DiffMod_PCtoDP([.5 .6 .7 .8 .9],6)
%
% returns:
% 
% dP =
% 
%     2.0188    2.4047    2.8258    3.3777    4.2375
%
% The example input arguments are an N=5 vector of proportion correct and
% a scalar M with a value of 6, and the output is an N=5 vector of d's
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.4.0, 1.6.3 (see History.m)


function dP=PAL_SDT_MAFCmatchSample_DiffMod_PCtoDP(pC,M)

[rows, cols]=size(pC);

func=@PAL_SDT_MAFCmatchSample_DiffMod_DPtoPC;

dP = zeros(rows,cols);

for r=1:rows
    for c=1:cols
        dP(r,c)=PAL_minimize(@PAL_sqDistanceYfuncX,1,[],pC(r,c),func,M);
        if pC(r,c)==1.0
            dP(r,c)=1/0;
        end
    end
end