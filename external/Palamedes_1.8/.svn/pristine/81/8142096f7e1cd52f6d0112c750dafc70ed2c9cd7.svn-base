%
%PAL_SDT_ROCML_RatioSDcomparison  tests whether the estimated ratio of
%   signal to noise SDs from an ROC rating scale experiment differs 
%   significantly from a user specified value
%
%Syntax:    
%   [TLR pTLR TLRSim converged] = PAL_SDT_ROCML_RatioSDcomparison(cumNumHF,... 
%       OutOfNum, Rtest, SDTF, invSDTF, B, {optional arguments});
%
%Input:
%
%   'cumNumHF': Nx2 matrix of cumulative numbers of 'hits' and 'false 
%       alarms', where each pair of hits and false alarms is the cumulative 
%       number up to that rating, for N+1 ratings.  The last N+1 
%       cumulative value, which equals the total number of trials for each 
%       stimulus alternative, should be omitted.  PAL_SDT_cumulateHF can be 
%       used to convert a matrix of number of hits and false alarms to 
%       their cumulative values
%       
%   'OutOfNum' is a Nx2 matrix of number of trials corresponding to each
%       value in cumNumHF.  PAL_SDT_cumulateHF returns OutOfNum as well as
%       cumNumHF
%
%   'Rtest': value of the ratio of signal-to-noise SDs to be tested
%
%   'SDTF' is the name of the appropriate Palamedes SDT function 
%       for the task for which the rating scale data has been obtained.  
%       Here it is an input argument to PAL_SDT_ROCML_Fit. At present the 
%       only Palamedes SDT function that allows a non-unity SD ratio as 
%       input is PAL_SDT_1AFC_PHFtoDP, which estimates d-prime for a Yes/No 
%       or symmetric 1AFC task.  SDTF needs to be defined before calling 
%       the present routine, for example by executing the following 
%       command:
%
%       SDTF = @PAL_SDT_1AFC_PHFtoDP;
%
%    'invSDTF' is the inverse of SDTF, that is the SDT function that 
%       converts a d-prime, criterion C and SD ratio into a proportion of 
%       hits and false alarms.  invSDTF is used as an input to the
%       parameter estimation function PAL_SDT_ROCML_Fit.  invSDTF needs to 
%       be defined before calling the present routine, by executing the 
%       following command:
%
%       invSDTF = @PAL_SDT_1AFC_DPtoPHF;
%
%   'B' is the number of simulations 
%
%Output:
%
%   'TLR': Transformed likelihood ratio 
%       [-2 x ln(Likelihood(lesser model)/Likelihood(fuller model)]
%
%   'pTLR': proportion of the B Transformed Likelihood Ratio (TLR) values 
%       from simulations that were greater than TLR value of data. The 
%       larger the value of pTLR, the better the fit to the model based
%       on the value of Rtest.  Thus if pTLR is low, say <0.05, one can
%       reject the hypothesis at the 0.05 level that the ratio-of-SDs
%       of the ROC data was no different than the value of Rtest
%
%   'TLRSim': vector of length B containing the TLR values of all
%       simulations (i.e., empirical TLR sampling distribution).
%
%   'converged': For each simulation contains a 1 in case the fit was
%       succesfull (i.e., converged) or a 0 in case it did not.
%
%   PAL_SDT_ROCML_RatioSDcomparison will generate a warning if not all 
%       simulations were fit succesfully.
%
%   PAL_SDT_ROCML_RatioSDcomparison uses Nelder-Mead Simplex method as 
%       implemented in PAL_minimize. The default search options may be
%       changed by using the optional argument 'searchOptions' followed by 
%       an options structure created using options = 
%       PAL_minimize('options'). For example of usage see 
%       PAL_SDT_ROCML_Demo. For more information type: 
%       PAL_minimize('options','help');
%
%Example:
%
%   cumNumHF = [159 2; 200 5; 219 26; 256 106];  
%   OutOfNum = [288 288; 288 288; 288 288; 288 288];
%   Rtest = 0.5;
%   SDTF = @PAL_SDT_1AFC_PHFtoDP;
%   invSDTF = @PAL_SDT_1AFC_DPtoPHF;
%   B = 400;
%
%   %Determine whether SD ratio is different from 0.5
%
%   [TLR pTLR] = PAL_SDT_ROCML_RatioSDcomparison(cumNumHF, OutOfNum, ...
%       Rtest, SDTF, invSDTF, B);
%
%returns something like:
%
%   TLR = 
%
%    0.2167
%
%   pTLR =
%
%    0.6450
%
%Introduced: Palamedes version 1.6.0 (FK & NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function [TLR, pTLR, TLRSim, converged] = PAL_SDT_ROCML_RatioSDcomparison(cumNumHF, OutOfNum, Rtest, SDTF, invSDTF, B, varargin)

options = [];

converged = false(B,1);
TLRSim = zeros(B,1);
pTLR = [];

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'searchOptions',7)
            options = varargin{n+1};
            valid = 1;
        end
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});
        end
    end            
end

[dP, R, C, negLLlesser] = PAL_SDT_ROCML_Fit(cumNumHF, OutOfNum, SDTF, invSDTF, 'ratioSDvalue', Rtest, 'ratioSDfree', 0, 'searchOptions', options);
paramsValuesLesser = [dP R C];

[trash, trash, trash, negLLfuller] = PAL_SDT_ROCML_Fit(cumNumHF, OutOfNum, SDTF, invSDTF, 'ratioSDvalue', 1, 'ratioSDfree', 1, 'searchOptions', options);

TLR = 2*(negLLlesser-negLLfuller);

for b = 1:B
    
    cumNumHF = PAL_SDT_ROC_SimulateObserverParametric(paramsValuesLesser, OutOfNum, invSDTF);
    
    [trash, trash, trash, negLLlesser, converged(b)] = PAL_SDT_ROCML_Fit(cumNumHF, OutOfNum, SDTF, invSDTF, 'ratioSDvalue', Rtest, 'ratioSDfree', 0, 'searchOptions', options);
    [trash, trash, trash, negLLfuller, converged(b)] = PAL_SDT_ROCML_Fit(cumNumHF, OutOfNum, SDTF, invSDTF, 'ratioSDvalue', Rtest, 'ratioSDfree', 1, 'searchOptions', options);
    
    if ~converged(b)
        warning('PALAMEDES:convergeFail','Fit to simulation %s of %s did not converge.',int2str(b), int2str(B));
    end
    
    TLRSim(b)=2*(negLLlesser-negLLfuller);    
    
end

exitflag = sum(converged) == B;
if exitflag ~= 1
    warning('PALAMEDES:convergeFail','Only %s of %s simulations converged.',int2str(sum(converged)), int2str(B));
end

if B > 0
    pTLR = length(TLRSim(TLRSim>TLR))/B;
end