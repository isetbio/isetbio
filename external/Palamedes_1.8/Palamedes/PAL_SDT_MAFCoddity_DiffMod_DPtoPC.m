%
% PAL_SDT_MAFCoddity_DiffMod_DPtoPC converts d'(d-prime) to 
% proportion correct for a M-AFC (M-alternative-forced-choice) 
% oddity task, assuming a Differencing model and an unbiased 
% observer
%
% Syntax: [pC]=PAL_SDT_MAFCoddity_DiffMod_DPtoPC(dP,M);
% 
% returns a scalar, vector or matrix of proportion correct ('pC') for a 
% scalar, vector or matrix of d' ('dP') and a scalar value of M (number of 
% alternatives), defined in the ranges 0<d'<inf and 3<M<inf. Note however
% that there will be an upper limit on M depending on the amount of 
% computer memory. Note also that because the routine uses Monte Carlo 
% simulation, no two outputs will be identical.  The speed/accuracy/memory 
% tradeoff can be altered by changing the parameter numReps, currently 
% set at 100000.  
%
% Example:
%
% [pC]=PAL_SDT_MAFCoddity_DiffMod_DPtoPC([0 1 2 3 4 5],5)
%
% returns something like:
% 
% pC =
% 
%     0.2004    0.3208    0.5860    0.8224    0.9446    0.9876
%
% The example input consists of a N=6 vector of d' with M set to 5, and 
% the output is the corresponding N=6 proportion correct.  
%
% Introduced: Palamedes version 1.8.0 (FK & NP)

function pC = PAL_SDT_MAFCoddity_DiffMod_DPtoPC(dP,M)

[rows, cols]=size(dP);
numReps=100000;

pC = zeros(rows,cols);

for r = 1:rows
    for c = 1:cols

        R = randn(M,numReps);
        R(1,:) = R(1,:)+dP(r,c);

        sqDiff = (R-repmat(mean(R),M,1)).^2;
        [maxVal, maxIndex]=max(sqDiff);

        pC(r,c) = length(maxIndex(maxIndex == 1))./numReps;

    end
end