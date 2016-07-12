%
% PAL_SDT_MAFCoddity_IndMod_DPtoPC converts d'(d-prime) to 
% proportion correct for a M-AFC (M-alternative-forced-choice) 
% oddity task, assuming an Independent Observer model and an 
% unbiased observer
%
% Implements the method described in Versfeld, Dai & Green (1996) 
% Perception & Psychophysics, 58, 10-51.  See Table 1, column e2 (epsilon
% squared) to check output.
%
% Syntax: [pC] = PAL_SDT_MAFCoddity_IndMod_DPtoPC(dP,M);
% 
% returns a scalar, vector or matrix of proportion correct ('pC') for a 
% scalar, vector or matrix of d' ('dP') and a scalar value of M (number of 
% alternatives), defined in the ranges 0<d'<inf and 2<M<inf. Note however
% that there will be an upper limit on M depending on the amount of 
% computer memory. Note also that because the routine uses Monte Carlo 
% simulation, no two outputs will be identical.  The speed/accuracy/memory 
% tradeoff can be altered by changing the parameter numReps using the
% optional argument 'numReps', followed by a positive integer indicating
% the number of Monte Carlo simulations to be used. Default value of
% numReps is 100000.
%
% Note that if M=3, you can use PAL_SDT_3AFCoddity_IndMod_DPtoPC, which
% is deterministic and is hence more accurate and faster
%
% Example:
%
% [pC] = PAL_SDT_MAFCoddity_IndMod_DPtoPC([0.5 1 1.5 2 2.5 3],5,...
%   'numReps',50000)
%
% might return: 
% 
% pC =
% 
%     0.2523    0.3903    0.5691    0.7324    0.8596    0.9355
%
% The example input consists of an N=6 vector of d' with M set to 5, and 
% the output is the corresponding N=6 proportion correct. 
%
% Introduced: Palamedes version 1.6.0 (FK & NP)
% Modified: Palamedes version 1.6.3 (see History.m)

function pC = PAL_SDT_MAFCoddity_IndMod_DPtoPC(dP,M,varargin)

[rows, cols] = size(dP);
numReps = 100000;

if ~isempty(varargin)
    valid = 0;
    if strncmpi(varargin{1}, 'numReps',4)
        numReps = varargin{2};
        valid = 1;
    end
    if valid == 0
        warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{1});
    end        
end            

Diag = logical(repmat(eye(M),[1,1,numReps]));
Diag = permute(Diag,[3,1,2]);

pC = zeros(rows,cols);

for r = 1:rows
    for c = 1:cols
        
        if dP(r,c) < 1e-4;
            pC(r,c) = 1/M;
        else
        
            xM = randn(numReps,M);
            xM = repmat(xM,[1,1,M]);
            xM(:,1,:) = xM(:,1,:) + dP(r,c);           

            L = zeros(size(Diag));
            L(Diag) = PAL_pdfNormal(xM(Diag),0,1);
            L(~Diag) = PAL_pdfNormal(xM(~Diag),dP(r,c),1);
            
            prodL = prod(L,2);

            L(~Diag) = PAL_pdfNormal(xM(~Diag),0,1);
            L(Diag) = PAL_pdfNormal(xM(Diag),dP(r,c),1);

            prodL = prodL + prod(L,2);

            [maxVal, maxIndex] = max(prodL,[],3);
            pC(r,c) = length(maxIndex(maxIndex == 1))./numReps;
            
        end 
    end
end