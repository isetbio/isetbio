%
%PAL_SDT_ROCML_GoodnessOfFit     Determines Goodness-of-Fit of ROC
% (receiver operating chracteristic)
%
%Syntax: [Dev pDev DevSim converged] = ...
%           PAL_SDT_ROCML_GoodnessOfFit(paramsValues, cumNumHF, OutOfNum,.. 
%           SDTF, invSDTF, B, varargin);
%
%Input:
%
%   'paramsValues': vector consisting of a d-prime value, a ratio-of-SDs
%       value and N values of criterion C.  These are the estimates
%       derived from PAL_SDT_ROCML_Fit

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
%   'SDTF' is the name of the appropriate Palamedes SDT function 
%       for the task for which the rating scale data has been obtained.  It 
%       is only used to calculate the initial guesses of d'prime and the
%       criteria C values that the parameter estimation function requires 
%       as input. The function specified by SDTF may take as input a 
%       non-unity ratio of signal and noise SDs. At present the 
%       only Palamedes SDT function that allows a non-unity SD ratio as 
%       input is PAL_SDT_1AFC_PHFtoDP, which estimates d-prime for a Yes/No 
%       or symmetric 1AFC task.  SDTF needs to be defined before calling 
%       the present routine, for example by executing the following command:
%
%       SDTF = @PAL_SDT_1AFC_PHFtoDP;
%
%    'invSDTF' is the inverse of SDTF, that is the SDT function that 
%       converts a d-prime, criterion C and SD ratio into a proportion of 
%       hits and false alarms.  invSDTF is used as an input to the
%       parameter estimation function.  invSDTF needs to be defined before 
%       calling the present routine, by executing the following 
%       command:
%
%       invSDTF = @PAL_SDT_1AFC_DPtoPHF;
%
%   'B' is the number of simulations 
%
%Output: 
%
%   'Dev': Deviance (transformed likelihood ratio comparing fit of
%       psychometric function to fit of saturated model)
%
%   'pDev': proportion of the B Deviance values from simulations that were
%       greater than Deviance value of data. The greater the value of pDev,
%       the better the fit.
%
%   'DevSim': vector containing all B simulated Deviance values.
%
%   'converged': For each simulation contains a 1 in case the fit was
%       succesfull (i.e., converged) or a 0 in case it did not.
%
%   PAL_SDT_ROCML_GoodnessOfFit will generate a warning if not all 
%       simulations were fit succesfully.
%
%   PAL_SDT_ROCML_GoodnessOfFit uses Nelder-Mead Simplex method as 
%       implemented in PAL_minimize. The default search options may be 
%       changed by using the optional argument 'searchOptions' followed by 
%       an options structure created using options = 
%       PAL_minimize('options'). For example of usage see 
%       PAL_SDT_ROCML_Demo. For more information type: 
%       PAL_minimize('options','help');
%
%Example:
% 
%   paramsValues = [1.7645 0.4663 0.7828 0.4867 0.1984 -0.2808];
%   cumNumHF = [159 2; 200 5; 219 26; 256 106];  
%   OutOfNum = [288 288; 288 288; 288 288; 288 288];
%   SDTF = @PAL_SDT_1AFC_PHFtoDP;
%   invSDTF = @PAL_SDT_1AFC_DPtoPHF;
%   B = 400;
%
%   %Determine Goodness-Of-Fit
%
%   [Dev pDev DevSim converged] = ...
%       PAL_SDT_ROCML_GoodnessOfFit(paramsValues, cumNumHF, OutOfNum,... 
%       SDTF, invSDTF, B);
%
%Introduced: Palamedes version 1.6.0 (FK & NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function [Dev, pDev, DevSim, converged] = PAL_SDT_ROCML_GoodnessOfFit(paramsValues, cumNumHF, OutOfNum, SDTF, invSDTF, B, varargin)

options = []; %default
Rfree = 1; %default
Rval = 1;  %default
converged = false(B,1);
DevSim = zeros(B,1);
pDev = [];

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'ratioSDvalue',9)
            Rval = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'ratioSDfree',9)
            Rfree = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'searchOptions',7)
            options = varargin{n+1};
            valid = 1;
        end
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});
        end
    end            
end

paramsFree = [1 Rfree ones(1,length(paramsValues)-2)];
paramsFreeVals = paramsValues(paramsFree == 1);
paramsFixedVals = paramsValues(paramsFree == 0);

negLLCon = PAL_SDT_ROCML_negLL(paramsFreeVals, paramsFixedVals, paramsFree, cumNumHF, OutOfNum, invSDTF);
negLLAug = PAL_SDT_ROCML_negLLNonParametric(cumNumHF, OutOfNum);

Dev = 2*(negLLCon-negLLAug);

for b = 1:B
    
    cumNumHF = PAL_SDT_ROC_SimulateObserverParametric(paramsValues, OutOfNum, invSDTF);
    
    [trash, trash, trash, negLLConSim, converged(b)] = PAL_SDT_ROCML_Fit(cumNumHF, OutOfNum, SDTF, invSDTF, 'ratioSDvalue', Rval, 'ratioSDfree', Rfree,'searchOptions',options);
    negLLAugSim = PAL_SDT_ROCML_negLLNonParametric(cumNumHF, OutOfNum); 
    
    DevSim(b) = 2*(negLLConSim-negLLAugSim);
    
    if ~converged(b)
        warning('PALAMEDES:convergeFail','Fit to simulation %s of %s did not converge.',int2str(b), int2str(B));
    end

end

exitflag = sum(converged) == B;
if exitflag ~= 1
    warning('PALAMEDES:convergeFail','Only %s of %s simulations converged.',int2str(sum(converged)), int2str(B));
end

if B > 0
    pDev = length(DevSim(DevSim>Dev))/B;
end