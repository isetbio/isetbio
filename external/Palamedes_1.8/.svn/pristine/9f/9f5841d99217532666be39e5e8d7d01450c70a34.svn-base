%
%PAL_ROC_BootstrapParametric    Perform parametric bootstrap to
%   determine standard errors on parameters of fitted ROC curve
%   
%syntax: [DPsd Rsd paramsSim LLsim converged] = ...
%   PAL_SDT_ROCML_BootstrapParametric(paramsValues, OutOfNum, B, SDTF,... 
%   invSDTF, varargin)
%
%Input:
%   'paramsValues': vector containing the fitted parameters provided by
%       PAL_SDT_ROCML_Fit.
%
%   'OutOfNum': matrix containing the total number of trials for each
%       rating and for both 'hits' and 'false alarms'.  Can be provided by
%       PAL_SDT_cumulateHF.
%
%   'B': number of bootstrap simulations to perform.
%
%   'SDTF': name of the function for calculating d-prime for the task.  
%       The function should take as input arguments proportions of 'hits' 
%       and 'false alarms' and have the option of specifying the ratio of 
%       noise-to-signal (or signalA-to-signalB) SDs.  At present the only 
%       Palamedes SDT function that permits non-unity ratio-of-SDs as input 
%       is PAL_SDT_1AFC_PHFtoDP, which estimates d-prime for a Yes/No 
%       or symmetric 1AFC task.  SDTF needs to be defined before calling 
%       the present routine, at present by executing the following 
%       command:  
% 
%       SDTF = @PAL_SDT_1AFC_PHFtoDP;
%
%   'invSDTF': the name of inverse of SDTF, that is the SDT function that 
%       converts a d-prime, criterion C and ratio-of-SDs into a proportion 
%       of hits and false alarms.  invSDTF needs to be defined before 
%       calling the present routine, at present by executing the 
%       following command:
%
%       invSDTF = @PAL_SDT_1AFC_DPtoPHF;
%
%  'B': number of bootstrap simulations to perform.
%
%Output:
%   
%   'dPsd': scalar containing the standard deviation of the d-prime 
%       parameters across the B fits to simulated data. These are 
%       estimates of the standard errors of the parameter estimates.
%
%   'Rsd':  scalar containing the standard deviation of the signal-to-noise
%       SD ratio estimated across the B fits to simulated data.  
%
%   'paramsSim': Listing of fitted parameter values across all B
%       simulations. May be useful e.g., to determine symmetry of
%       distribution, to find outliers, etc. B x 2 matrix, first column
%       lists simulated dP values, second column lists simulated R values.
%
%   'LLSim': vector containing Log Likelihoods associated with all B fits
%       to simulated data.
%
%   'converged': For each simulation contains a 1 in case the fit was
%       succesfull (i.e., converged) or a 0 in case it did not.
%
%   PAL_ROC_BootstrapParametric will generate a warning if not all 
%       simulations were fit succesfully.
%
%Optional input arguments:
%
%   User may provide an initial guess of the ratio of noise-to-signal SDs, 
%       in the range 0<R<Inf.  If not specified, 'R' defaults to 1. To 
%       include this argument use the following syntax, where R gives the 
%       initial guess:
%
%  [DPsd Rsd paramsSim LLsim converged] = ...
%   PAL_SDT_ROCML_BootstrapParametric(paramsValues, OutOfNum, B, SDTF,... 
%     invSDTF,'ratioSDvalue',R)
%
%   User may constrain the ratio of noise-to-signal SDs to equal a fixed 
%       value (e.g. 1). This may be useful to test whether the ratio 
%       differs significantly from 1 (or some other set value). To fix the 
%       ratio-of-SD use the following syntax:
%
%   [DPsd Rsd paramsSim LLsim converged] = ...
%     PAL_SDT_ROCML_BootstrapParametric(paramsValues, OutOfNum, B, SDTF,... 
%     invSDTF,'ratioSDfree',0).
%
%   In order to fix the SD ratio to some value other than 1, use the
%       'ratioSDfree'= 0 option in conjunction with the 'ratioSDvalue' 
%       option.
%
%   PAL_ROCML_BootstrapParametric uses the Nelder-Mead Simplex method as 
%       implemented in PAL_minimize. The default search options may be 
%       changed by using the optional argument 'searchOptions' followed by 
%       an options structure created using options = 
%       PAL_minimize('options'). For example of usage see 
%       PAL_SDT_ROCML_Demo. For more information type: 
%       PAL_minimize('options','help');
%
%  
%Example:
%
%   SDTF = @PAL_SDT_1AFC_PHFtoDP;
%   invSDTF = @PAL_SDT_1AFC_DPtoPHF;
%
%   NumHF = [159 2; 41 3; 19 21; 37 80; 32 182]; %input data
%
%   [cumNumHF OutOfNum pHF] = PAL_SDT_cumulateHF(NumHF); %cumulate data
%
%   [dP R C negLL exitflag] = PAL_SDT_ROCML_Fit(cumNumHF,OutOfNum,SDTF,...
%     invSDTF,'ratioSDvalue',0.75);  %fit ROC parameters
%
%   paramsValues = [dP R C]; %set up params for bootstrap routine
%   B=400;
%
%   % perform bootstrap
%
%   [dPsd Rsd] = ...
%     PAL_SDT_ROCML_BootstrapParametric(paramsValues, OutOfNum, SDTF,...
%     invSDTF, B)
%
%returns something like this:
%
%   dPsd =
%
%    0.1032
%
%   Rsd =
%
%    0.0760
%
%Introduced: Palamedes version 1.6.0 (FK & NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function [dPsd, Rsd, paramsSim, LLsim, converged] = PAL_SDT_ROCML_BootstrapParametric(paramsValues, OutOfNum, SDTF, invSDTF, B, varargin)

options = [];

Rfree=1; %default
Rval=1; %default

LLsim = zeros(B,1);
converged = false(B,1);
dP = zeros(B,1);
R = zeros(B,1);
C = zeros(B,length(paramsValues)-2);

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

for b = 1:B

    %Simulate experiment
    cumNumHF = PAL_SDT_ROC_SimulateObserverParametric(paramsValues, OutOfNum, invSDTF);
    
    [dP(b), R(b), C(b,:), LLsim(b), converged(b)] = PAL_SDT_ROCML_Fit(cumNumHF, OutOfNum, SDTF, invSDTF,'ratioSDvalue', Rval, 'ratioSDfree', Rfree, 'searchOptions',options);
    
    if ~converged(b)
        warning('PALAMEDES:convergeFail','Fit to simulation %s of %s did not converge.',int2str(b), int2str(B));
    end
    
end

exitflag = sum(converged) == B;
if exitflag ~= 1
    warning('PALAMEDES:convergeFail','Only %s of %s simulations converged.',int2str(sum(converged)), int2str(B));
end

paramsSim(:,1) = dP;
paramsSim(:,2) = R;

[dPmean, dPsd] = PAL_MeanSDSSandSE(dP);
[Rmean, Rsd] = PAL_MeanSDSSandSE(R);