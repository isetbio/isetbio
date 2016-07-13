%
% PAL_SDT_MAFCoddity_DiffMod_PCtoDP converts proportion correct into 
% d'(d-prime) for a M-AFC (M-alternative-forced-choice) oddity
% task, assuming a Differencing Observer model and an unbiased observer
%
% Syntax: [dP]=PAL_SDT_MAFCoddity_DiffMod_PCtoDP(pC,M);
% 
% returns a scalar, vector or matrix of d' ('dP') for an input scalar, 
% vector or matrix of proportion correct p ('pC') and an input scalar of M 
% (number of alternatives), defined in the ranges 0<p<1 and 3<M<inf.  Note 
% that the routine may take several seconds or even minutes to
% execute depending on the size of the input and speed of computer, as
% the routine performs an iterative search on the inverse version of the 
% routine which employs Monte Carlo simulation
%
% Example:
%
% [dP]=PAL_SDT_MAFCoddity_DiffMod_PCtoDP([.5 .6 .7 .8 .9],5)
%
% returns something like:
% 
% dP =
% 
%     1.6922    2.0469    2.4281    2.8797    3.5258
%
% The example input arguments are an N=5 vector of proportion correct and
% a scalar M with a value of 5, and the output is an N=5 vector of d's.  
% 
% Introduced: Palamedes version 1.8.0 (FK)


function dP=PAL_SDT_MAFCoddity_DiffMod_PCtoDP(pC,M)

[rows, cols]=size(pC);

func=@PAL_SDT_MAFCoddity_DiffMod_DPtoPC;

dP = zeros(rows,cols);

for r=1:rows
    for c=1:cols
        dP(r,c)=PAL_minimize(@PAL_sqDistanceYfuncX,1,[],pC(r,c),func,M);
        if pC(r,c)==1.0
            dP(r,c)=1/0;
        end
    end
end