%
%PAL_SDT_ROCML_Fit   Fits the parameters of an ROC curve from 
% rating scale data using a Maximum Likelihood criterion.
%
% Syntax: 
% [dP R C negLL exitflag output] = PAL_SDT_ROCML_Fit(cumNumHF,OutOfNum,...
%   SDTF,invSDTF,{optional arguments})
%
%Input:
%
%   'cumNumHF': Nx2 matrix of cumulative numbers of 'hits' and 'false 
%       alarms'. Each pair of hits and false alarms is the cumulative 
%       number for each rating, for N+1 ratings, where the rating starts 
%       at high certainty that the trial was a signal (or signal 1) and 
%       ends at high certainty that the trial was noise (or signal2).  
%       The last N+1 cumulative value, which equals the total number of 
%       trials for each stimulus alternative, must be omitted.  
%       PAL_SDT_cumulateHF can be used to convert a matrix of number hits 
%       and false alarms to their cumulative values
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
%Output:
%
%   'dP'  Parameter estimate of d-prime (d') for the ROC curve.
%   
%   'R': Parameter estimate of the slope of the ROC curve, also the ratio
%       of noise (or signal 1) to signal (or signal 2) SDs (standard 
%       deviations). If 'ratioSDfree' is set to 'fixed' (see below) 
%       the returned 'R' will be the fixed parameter value.
%
%   'C' is a N-length vector of criterion values corresponding to each
%       pair of hits and false alarms
%
%   'negLL': Negative Log likelihood associated with the fit.
%
%   'exitflag': 1 indicates a succesful fit, 0 indicates fit did not
%       converge (trying again using new initial guesses might help).
%
%   'output': message containing some information concerning fitting
%       process.
%
%Optional input arguments:
%
%   User may provide an initial guess of the ratio of noise (or signal 1) 
%       to signal (or signal 2) standard deviations, in the range 0<R<inf.  
%       If not specified, 'R' defaults to 1. To include this argument use 
%       the following syntax, where R gives the initial guess:
%
%   [DP R C negLL exitflag output] = PAL_SDT_ROCML_Fit(cumNumHF,OutOfNum,..
%     SDTF,invSDTF,'ratioSDvalue',R)
%
%   User may constrain the ratio of noise-to-signal (or signalA-to-signalB) 
%       SDs to equal a fixed value (e.g. 1). This may be useful to test 
%       whether the ratio differs significantly from 1 or some other set 
%       value. To fix the SD ratio use the following syntax:
%
%   [dP R C LL exitflag output] = SDT_1AFC_ROCML_Fit(cumNumHF,OutOfNum...
%     'ratioSDfree',0)
%
% Note: In order to fix the SD ratio to some value other than 1, use the
% 'ratioSDfree'= 0 option in conjunction with the 'ratioSDvalue' 
% option.
%
% PAL_SDT_ROCML_Fit uses Nelder-Mead Simplex method as implemented in 
%   PAL_minimize. The default search options may be changed 
%   by using the following syntax:
%
% [dP R C LL exitflag output] = PAL_SDT_ROCML_Fit(cumNumHF,OutOfNum,...
%   SDTF,invSDTF,'searchOptions', options) 
%
%   where 'options' is a structure created using: options = 
%   PAL_minimize('options'). For example of usage see PAL_SDT_ROCML_Demo.
%   For more information type PAL_minimize('options','help').
%
% Example: The cumulative hit and false alarm data below are taken from 
%   Table 5.4 in McNicol, D. (2005) A Primer of Signal Detection Theory: 
%   Lawrence Erlbaum Associates. The input ratioSDvalue is the initial 
%   guess of the ratio of signal to noise SDs
%    
%   cumNumHF=[159 2; 200 5; 219 26; 256 106];  
%   OutOfNum=[288 288; 288 288; 288 288; 288 288];
%   SDTF=@PAL_SDT_1AFC_PHFtoDP;
%   invSDTF=@PAL_SDT_1AFC_DPtoPHF;
%   ratioSDvalue=0.6;
%
%   Fit data:
%
%   [dP R C negLL exitflag] = PAL_SDT_ROCML_Fit(cumNumHF,OutOfNum,...
%       SDTF,invSDTF,'ratioSDvalue',0.6)
%
% returns:
% 
% dP =
% 
%     1.7645
% 
% 
% R =
% 
%     0.4663
% 
% 
% C =
%
%    0.7828    0.4867    0.1984   -0.2808
% 
% 
% 
% negLL =
% 
%  -949.1988
% 
% 
% exitflag =
% 
%      1
%   
% In the output the two fitted parameters are d-prime ('dP') and the ratio
% of SDs ('R').  'C' gives the fitted values of criterion C for the 4 data
% points along the ROC curve. negLL gives the Maximum negative Log 
% Likelihood of the fit, and exitflag is 1 if the fit was successful, 
% 0 otherwise
%
%Introduced: Palamedes version 1.6.0 (FK & NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function [dP, R, C, negLL, exitflag, output] = PAL_SDT_ROCML_Fit(cumNumHF,OutOfNum,SDTF,invSDTF,varargin)

options = []; %default
Rval = 1; %default
Rfree = 1; %default


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
        if strncmpi(varargin{n}, 'SearchOptions',7)
            options = varargin{n+1};
            valid = 1;
        end
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});
        end
    end            
end


% Make initial guess of the d-prime value by geometrically averaging 
% the d-primes for each pHF.  Make initial guesses of criterion Cs for each 
% pHF. Use default guess of R=1 unless otherwise specified in the input
Rvec=repmat(Rval,1,length(cumNumHF)); % convert scalar R to a vector for input to routine below
pHF=cumNumHF./OutOfNum;
pHF(pHF == 0) = 0.0001;
pHF(pHF == 1) = 0.9999;
[dP, C]=SDTF(pHF,'ratioSDvalue',Rvec);
dP(dP<=0)=0.0001;
guessDP=10.^mean(log10(dP'));
guessC=C';
guessR=Rval;

% Set up free parameters with initial guesses
paramsValues = [guessDP guessR guessC];
paramsFree = [1 Rfree ones(size(guessC))];

paramsFreeVals = paramsValues(paramsFree == 1);
paramsFixedVals = paramsValues(paramsFree == 0);

[paramsFreeVals, negLL, exitflag, output] = PAL_minimize(@PAL_SDT_ROCML_negLL, paramsFreeVals, options, paramsFixedVals, paramsFree, cumNumHF, OutOfNum, invSDTF); 

paramsValues(paramsFree == 1) = paramsFreeVals;
paramsValues(paramsFree == 0) = paramsFixedVals;

dP = paramsValues(1);
R = paramsValues(2);
C = paramsValues(3:length(paramsValues));



