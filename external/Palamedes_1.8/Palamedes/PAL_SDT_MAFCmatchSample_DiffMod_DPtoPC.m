%
% PAL_SDT_MAFCmatchSample_DiffMod_DPtoPC converts d'(d-prime) to 
% proportion correct for a M-AFC (M-alternative-forced-choice) 
% match-to-sample task, assuming a Differencing model and an unbiased 
% observer
%
% Syntax: [pC]=PAL_SDT_MAFCmatchSample_DiffMod_DPtoPC(dP,M);
% 
% returns a scalar, vector or matrix of proportion correct ('pC') for a 
% scalar, vector or matrix of d' ('dP') and a scalar value of M (number of 
% alternatives), defined in the ranges 0<d'<inf and 2<M<inf. Note however
% that there will be an upper limit on M depending on the amount of 
% computer memory. Note also that because the routine uses Monte Carlo 
% simulation, no two outputs will be identical.  The speed/accuracy/memory 
% tradeoff can be altered by changing the parameter numReps, currently 
% set at 100000.  
%
% Example:
%
% [pC]=PAL_SDT_MAFCmatchSample_DiffMod_DPtoPC([0 1 2 3 4 5],7)
%
% returns: 
% 
% pC =
% 
%     0.1431    0.2293    0.4669    0.7152    0.8693    0.9443
%
% The example input consists of an N=6 vector of d' with M set to 7, and 
% the output is the corresponding N=6 proportion correct. 
%
% Introduced: Palamedes version 1.0.0 (FK & NP)
% Modified: Palamedes version 1.0.1, 1.6.3 (see History.m)

function pC = PAL_SDT_MAFCmatchSample_DiffMod_DPtoPC(dP,M)

[rows, cols]=size(dP);
numReps=100000;

pC = zeros(rows,cols);

for r = 1:rows
    for c = 1:cols

        RT = randn(1,numReps)+dP(r,c);
        RM = randn(M,numReps);
        RM(1,:) = RM(1,:)+dP(r,c);

        sqDiff = (repmat(RT,M,1)-RM).^2;
        [minVal, minIndex]=min(sqDiff,[],1);
       
        pC(r,c) = length(minIndex(minIndex == 1))./numReps;

    end
end