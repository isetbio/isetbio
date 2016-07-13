%
%PAL_PFLR_ModelComparison   Statistically compare two models fitted to 
%   multi-condition experiment using Monte-Carlo simulations.
%
%Syntax: [TLR pTLR paramsL paramsF TLRSim converged] = ...
%   PAL_PFLR_ModelComparison(StimLevels, NumPos, OutOfNum, ...
%   paramsValues, B, PF, {optional arguments})
%
%PAL_PFLR_ModelComparison compares two models which specify constraints on
%   the parameter values of Psychometric Functions (PFs) between multiple
%   conditions in experiment, the 'lesser' model and the 'fuller' model.
%   The lesser model should be nested under the fuller model. Briefly,
%   PAL_PFLR_ModelComparison calculates a transformed likelihood ratio for
%   the model fits to the data and compares this against an empirical
%   sampling distribution of transformed likelihood ratios obtained through
%   Monte Carlo simulations. Simulated datasets are generated according to 
%   the assumptions of the lesser model.
%
%Input:
%   'StimLevels': matrix containing stimulus levels used. Each row of the
%       matrix corresponds to a condition.
%
%   'NumPos': matrix containing for each of the entries of 'StimLevels' the 
%       number of trials a positive response (e.g., 'yes' or 'correct') was
%       given.
%
%   'OutOfNum': matrix containing for each of the entries of 'StimLevels' 
%       the total number of trials. (May contain zeros in case not all 
%       conditions utilized identical numbers of stimulus levels).
%
%   'paramsValues': Matrix containing parametervalues. Each row should 
%       correspond to a condition. Matrix should have four columns 
%       [threshold slope guess-rate lapse-rate]. Entries should indicate 
%       initial guesses for free parameters, fixed values for fixed 
%       parameters. Alternatively, you may supply a 1 x 4 vector in which
%       case the entries will be used for all conditions as initial search
%       values.
%
%   'B': Number of Monte Carlo simulations to perform.
%
%   'PF': psychometric function to be fitted. Passed as an inline function.
%       Options include:    
%           @PAL_Logistic
%           @PAL_Weibull
%           @PAL_Gumbel (i.e., log-Weibull)
%           @PAL_Quick
%           @PAL_logQuick
%           @PAL_CumulativeNormal
%           @PAL_Gumbel
%           @PAL_HyperbolicSecant
%
%Output:
%   'TLR': Transformed likelihood ratio 
%       [-2 x ln(Likelihood(lesser model)/Likelihood(fuller model)]
%
%   'pTLR': proportion of the B simulated TLRs that were greater than the
%       TLR of the dataset. In other words, it is an (estimate of) the
%       probability that an observer characterized by the lesser model 
%       would generate data resulting in as high a TLR as the human 
%       observer's TLR. By convention, if pTLR < 0.05 we conclude that the
%       lesser model does not describe human observer well.
%
%   'paramsL': parameter estimates under the lesser model.
%
%   'paramsF': parameter estimates under the fuller model.
%
%   'TLRSim': vector of length B containing the TLR values of all
%       simulations (i.e., empirical TLR sampling distribution).
%
%   'converged': For each simulation contains a 1 in case the fit was
%       succesfull (i.e., converged) or a 0 in case it did not.
%
%   'funcParamsL': parametrization structure containing parameter
%       estimates of custom-defined parameters in lesser model (if any). 
%       Type 'help PAL_PFML_CustomDefine' for more information.
%
%   'funcParamsF': parametrization structure containing parameter
%       estimates of custom-defined parameters in fuller model (if any). 
%       Type 'help PAL_PFML_CustomDefine' for more information.
%
%PAL_PFLR_ModelComparison will generate a warning if not all 
%   simulations were fit succesfully.
%
%   In case not all fits converge succesfully, use optional argument 
%   'maxTries' to set the maximum number of times each fit will be 
%   attempted. The first try uses initial search values equal to 
%   paramsValues provided in function call, subsequent tries use these 
%   search values plus some random 'jitter'. The range of the random 
%   jitter can be set for each of the four parameters separately using 
%   the optional argument 'rangeTries'. 'rangeTries' should be set to
%   a 1 x 4 vector containing the range of jitter to be applied to 
%   initial guesses of parameters [alpha beta gamma lambda]. Default 
%   value for 'maxTries' is 1, default value for 'rangeTries' is the 
%   entirely arbitrary and in most cases inappropriate [1 1 1 1].
%   Jitter will be selected randomly from rectangular distribution
%   centered on guesses in 'paramsValues' with range in 'rangeTries'. In
%   case some or all parameters are custom-reparametrized, please type
%   'help PAL_PFML_CustomDefine' for information on proper use of
%   'rangeTries' argument.
%
%   (Note that some simulated data sets may never be fit succesfully 
%       regardless of value of 'maxTries' and 'rangeTries') 
%
%User may specify the lesser and fuller model using optional arguments 
%   in much the same way as the model to be fitted is to be specified in 
%   PAL_PFML_FitMultiple (type 'help PAL_PFML_FitMultiple' or refer to 
%   example below). However, since here we need to define both a lesser and 
%   a fuller model the arguments need to specify which model they apply to 
%   and they are (default values indicated in curly brackets):
%       'lesserThresholds'      {'constrained'}
%       'lesserSlopes'          {'constrained'}
%       'lesserGuessrates'      {'fixed'}
%       'lesserLapserates'      {'fixed'}
%       'fullerThresholds'      {'unconstrained'}
%       'fullerSlopes'          {'unconstrained'}
%       'fullerGuessrates'      {'fixed'}
%       'fullerLapserates'      {'fixed'}
%
%If highest entries in 'StimLevels' are so high that it can be assumed that
%   errors observed there can be due only to lapses, use 'lapseFit'
%   argument to specify alternative fitting scheme. Options: 'nAPLE' 
%   (default), 'iAPLE', and 'jAPLE'. Type help PAl_PFML_FitMultiple for 
%   more information.
%
%The guess rate and lapse rate parameter can be constrained to be equal, as
%   would be appropriate, for example, in a bistable percept task. To
%   accomplish this, use optional argument 'gammaEQlambda', followed by a 
%   1. Both the guess rate and lapse rate parameters will be fit according
%   to options set for the lapse rate parameter. Value for guess rate in
%   'paramsValues' will be ignored.
%
%User may constrain the lapse rate to fall within limited range using 
%   the optional argument 'lapseLimits', followed by two-element vector 
%   containing lower and upper limit respectively. See full example below.
%
%User may constrain the guess rate to fall within limited range using 
%   the optional argument 'guessLimits', followed by two-element vector 
%   containing lower and upper limit respectively.
%
%PAL_PFLR_ModelComparison uses Nelder-Mead Simplex method to find 
%   the maximum in the likelihood function. The default search options may 
%   be changed by using the optional argument 'searchOptions' followed by 
%   an options structure created using options = PAL_minimize('options'). 
%   See example of usage below. For more information type:
%   PAL_minimize('options','help');
%
%Example:
%
%   options = PAL_minimize('options');  %decrease tolerance (i.e., increase
%   options.TolX = 1e-09;              %precision). This is a good idea,
%   options.TolFun = 1e-09;            %especially in high-dimension
%   options.MaxIter = 10000;           %parameter space.
%   options.MaxFunEvals = 10000;
% 
%   StimLevels = [-2:1:2; -2:1:2; -2:1:2; -2:1:2];
%   NumPos = [158 177 222 268 275;
%           196 232 282 351 380;
%           159 161 218 253 283;
%           213 245 296 356 375];
%   OutOfNum = [300 300 300 300 300;
%            400 400 400 400 400;
%            300 300 300 300 300;
%            400 400 400 400 400];
%
%   PF = @PAL_Logistic;
% 
%   paramsValues = [0 1.5 .5 .02]; %or: repmat([0 1.5 .5 .02],[4 1])
% 
%   maxTries = 4;                 
%   rangeTries = [1 1 0 0];
%   B = 400;
%   
%   %default comparison (thresholds AND slopes equal, while guess rates and
%   %lapse rates fixed
%   
%   [TLR pTLR paramsL paramsF TLRSim converged] = ...
%       PAL_PFLR_ModelComparison(StimLevels, NumPos, OutOfNum, ...
%       paramsValues, B, PF,'maxTries',maxTries,'rangeTries',rangeTries,...
%       'searchOptions',options);   
%   
%   %Same but now fitting a common lapse rate based solely on errors at 
%   highest stimulus intensity, assuming these errors can only be due to 
%   lapses. Type help PAL_PFML_Fit for more information on this.
%   
%   [TLR pTLR paramsL paramsF TLRSim converged] = ...
%       PAL_PFLR_ModelComparison(StimLevels, NumPos, OutOfNum, ...
%       paramsValues, B, PF, 'lesserlapse','constrained','fullerlapse',...
%       'constrained','maxTries',maxTries,'rangeTries',rangeTries,...
%       'searchOptions',options,'lapseFit','iAPLE');
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.0.2, 1.1.0, 1.2.0, 1.3.0, 1.4.0, 1.6.0, 
%   1.6.3 (see History.m)

function [TLR, pTLR, paramsL, paramsF, TLRSim, converged, funcParamsL, funcParamsF] = PAL_PFLR_ModelComparison(StimLevels, NumPos, OutOfNum, paramsValues, B, PF, varargin)

[StimLevels, NumPos, OutOfNum] = PAL_PFML_GroupTrialsbyX(StimLevels, NumPos, OutOfNum);

if size(paramsValues,1) == 1
    paramsValues = repmat(paramsValues,[size(StimLevels,1) 1]);
end

options = [];

alphas = paramsValues(:,1)';
betas = paramsValues(:,2)';
gammas = paramsValues(:,3)';
lambdas = paramsValues(:,4)';

gammaEQlambda = logical(false);
lapseFit = 'nAPLE';
lapseLimits = [];
guessLimits = [];

TLRSim = zeros(B,1);
converged = false(B,1); 

maxTries = 1;
rangeTries = [1 1 .1 .1];

MC = PAL_PFLR_setupMC(paramsValues);

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'SearchOptions',6)
            options = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'maxTries',4)            
            maxTries = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'rangetries',6)
            rangeTries = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'lapseLimits',6)
            lapseLimits = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'guessLimits',6)
            guessLimits = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'LesserThre',10)
            MC.argsAlesser = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'LesserSlop',10)
            MC.argsBlesser = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'LesserGues',10)
            MC.argsGlesser = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'LesserLaps',10)
            MC.argsLlesser = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'fullerThre',10)
            MC.argsAfuller = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'fullerSlop',10)
            MC.argsBfuller = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'fullerGues',10)
            MC.argsGfuller = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'fullerLaps',10)
            MC.argsLfuller = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'lapseFit',6)
            lapseFit = varargin{n+1};
            valid = 1;
        end                    
        if strncmpi(varargin{n}, 'gammaEQlambda',6)
            gammaEQlambda = logical(varargin{n+1});            
            valid = 1;
        end                                        
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});    
        end        
    end            
end

if (strncmpi(lapseFit,'iap',3) || strncmpi(lapseFit,'jap',3)) && ((PAL_whatIs(MC.argsLlesser) == 2 && strncmpi(MC.argsLlesser,'fix',3)) || (PAL_whatIs(MC.argsLfuller) == 2 && strncmpi(MC.argsLfuller,'fix',3)) || isempty(MC.argsLlesser) || isempty(MC.argsLfuller))
    warning('PALAMEDES:invalidOption','Lapse rates are not free in lesser or fuller model: ''lapseFit'' argument ignored');    
    lapseFit = 'nap';
end
if (strncmpi(lapseFit,'iap',3) || strncmpi(lapseFit,'jap',3)) && (PAL_whatIs(MC.argsLlesser) == 4 || PAL_whatIs(MC.argsLfuller) == 4)
    warning('PALAMEDES:invalidOption','Lapse rates custom-reparameterized in lesser or fuller model: ''lapseFit'' argument ignored');    
    lapseFit = 'nap';
end

if gammaEQlambda
    MC.argsGlesser = [];
    MC.argsGfuller = [];
    if ~isempty(guessLimits)
        warning('PALAMEDES:invalidOption','Guess rates constrained to equal lapse rates: ''guessLimits'' ignored');        
        guessLimits = [];
    end    
end

[TLR, exitflagDat, paramsL, paramsF, funcParamsL, funcParamsF] = PAL_PFLR_TLR(StimLevels, NumPos, OutOfNum, paramsValues, PF, MC,'searchoptions',options,'maxtries',maxTries,'rangetries',rangeTries,'lapseLimits',lapseLimits,'guessLimits',guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);

if ~exitflagDat
    TLRSim = [];
    converged = []; 
    pTLR = [];
    exitflag = 0;
    warning('PALAMEDES:convergeFailAbort','Fit to data did not converge. Exiting.');
else
    if isstruct(MC.argsAlesser)
        MC.argsAlesser = funcParamsL;
    end
    if isstruct(MC.argsAfuller)
        MC.argsAfuller = funcParamsF;
    end
    if isstruct(MC.argsBlesser)
        MC.argsBlesser = funcParamsL;
    end
    if isstruct(MC.argsBfuller)
        MC.argsBfuller = funcParamsF;
    end
    if isstruct(MC.argsGlesser)
        MC.argsGlesser = funcParamsL;
    end
    if isstruct(MC.argsGfuller)
        MC.argsGfuller = funcParamsF;
    end
    if isstruct(MC.argsLlesser)
        MC.argsLlesser = funcParamsL;
    end
    if isstruct(MC.argsLfuller)
        MC.argsLfuller = funcParamsF;
    end
    
    for b = 1:B
        
        for Cond = 1:size(StimLevels,1) 
            NumPos(Cond,:) = PAL_PF_SimulateObserverParametric(paramsL(Cond,:), StimLevels(Cond,:), OutOfNum(Cond,:), PF,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
        end

        [TLRSim(b), converged(b)] = PAL_PFLR_TLR(StimLevels, NumPos, OutOfNum, paramsValues, PF, MC,'searchoptions',options,'maxtries',maxTries,'rangetries',rangeTries,'lapseLimits',lapseLimits,'guessLimits',guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
        if ~converged(b)
            warning('PALAMEDES:convergeFail','Fit to simulation %s of %s did not converge.',int2str(b), int2str(B));
        end        
    end

    pTLR = length(TLRSim(TLRSim > TLR))./B;
    exitflag = (sum(converged) == B);
    if ~exitflag
        warning('PALAMEDES:convergeFail','Only %s of %s simulations converged.',int2str(sum(converged)), int2str(B));
    end
end