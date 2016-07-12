%
%PAL_PFML_GoodnessOfFitMultiple     Determine Goodness-of-Fit of set of 
%   psychometric functions (PFs) fitted to multiple conditions.
%
%Syntax:    [Dev pDev DevSim converged] = ...
%           PAL_PFML_GoodnessOfFitMultiple(StimLevels, NumPos, ...
%           OutOfNum, paramsValues, B, PF,{optional arguments});
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
%       [threshold slope guess-rate lapse-rate]. Use PAL_PFML_FitMultiple
%       to find best-fitting parameter values.
%
%   'B': number of bootstrap simulations to perform.
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
%   'Dev': Deviance (transformed likelihood ratio comparing fit of
%       psychometric functions to fit of saturated model)
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
%   PAL_PFML_GoodnessOfFitMultiple will generate a warning if not all 
%       simulations were fit succesfully.
%
%   PAL_PFML_GoodnessOfFitMultiple will accept a few optional arguments:
%
%       By default, PAL_PFML_GoodnessOfFitMultiple will test the goodness-
%       of-fit of a model which assumes that thresholds are equal between
%       datasets, slopes are equal between datasets and guess-rates and 
%       lapse-rates are fixed. User may specify alternative models in a 
%       manner similar to that in PAL_PFML_FitMultiple. Type 'help
%       PAL_PFML_FitMultiple' for more informations. Also see example
%       below.
%
%   If highest entries in 'StimLevels' are so high that it can be assumed 
%       that errors observed there can be due only to lapses, use 
%       'lapseFit' argument to specify alternative fitting scheme. Options: 
%       'nAPLE' (default), 'iAPLE', and 'jAPLE'. Type help 
%       PAl_PFML_FitMultiple for more information.
%
%   The guess rate and lapse rate parameter can be constrained to be equal, 
%       as would be appropriate, for example, in a bistable percept task. 
%       To accomplish this, use optional argument 'gammaEQlambda', followed 
%       by a 1. Both the guess rate and lapse rate parameters will be fit 
%       according to options set for the lapse rate parameter. Value for 
%       guess rate in 'paramsValues' will be ignored.
%
%       In case not all fits converge succesfully, use optional argument 
%       'maxTries' to set the maximum number of times each fit will be 
%       attempted. The first try uses initial search values equal to 
%       paramsValues provided in function call, subsequent tries use these 
%       search values plus some random 'jitter'. The range of the random 
%       jitter can be set for each of the four parameters separately using 
%       the optional argument 'rangeTries'. 'rangeTries' should be set to
%       a 1 x 4 vector containing the range of jitter to be applied to 
%       initial guesses of parameters [alpha beta gamma lambda]. Default 
%       value for 'maxTries' is 1, default value for 'rangeTries' is the 
%       entirely arbitrary and in most cases inappropriate [1 1 1 1].
%       Jitter will be selected randomly from rectangular distribution
%       centered on guesses in 'paramsValues' with range in 'rangeTries'.
%
%       (Note that some simulated data sets may never be fit succesfully 
%           regardless of value of 'maxTries' and 'rangeTries') 
%
%       User may constrain the lapse rate to fall within limited range 
%       using the optional argument 'lapseLimits', followed by a two-
%       element vector containing lower and upper limit respectively. See 
%       full example below.
%
%       User may constrain the guess rate to fall within limited range 
%       using the optional argument 'guessLimits', followed by a two-
%       element vector containing lower and upper limit respectively.
%
%   PAL_PFML_GoodnessOfFitMultiple uses Nelder-Mead Simplex method to find 
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
%   PF = @PAL_Logistic;
%   paramsValues = [0 1.5 .5 .04];
% 
% %Fitted using defaults (thresholds and slopes unconstrained, guess and 
%   lapse rates fixed):
% 
%   paramsFitted = PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, ...
%       paramsValues, PF,'searchOptions',options);
% 
%   %Determine Goodness-Of-Fit
% 
%   B = 400;
%   
%   [Dev pDev DevSim converged] = PAL_PFML_GoodnessOfFitMultiple(StimLevels, ...
%       NumPos, OutOfNum, paramsFitted, B, PF,'thresholds','uncons',...
%       'slopes','uncons','searchOptions',options);
%
%Fitted using different model:
% 
%   paramsFitted = PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, ...
%   paramsValues, PF,'lapserates',[1 1 1 1;1 1 -1 -1],'lapsefit','iAPLE',...
%   'searchOptions',options);
% 
% %Determine Goodness-Of-Fit (might take a while):
%
%   maxTries = 4;
%   rangeTries = [1 1 0 0];
% 
%   [Dev pDev DevSim converged] = PAL_PFML_GoodnessOfFitMultiple(StimLevels, ...
%       NumPos, OutOfNum, paramsFitted, B, PF,'lapserates',...
%       [1 1 1 1;1 1 -1 -1],'lapsefit','iAPLE','thresholds','uncons','slopes',...
%       'uncons','maxtries', maxTries,'rangeTries',rangeTries,...
%       'searchOptions',options);
% 
%Introduced: Palemedes version 1.0.0 (NP)
% Modified: Palamedes version 1.0.2, 1.1.0, 1.2.0, 1.3.0, 1.4.0, 1.6.0, 
%   1.6.3 (see History.m)

function [Dev, pDev, DevSim, converged] = PAL_PFML_GoodnessOfFitMultiple(StimLevels, NumPos, OutOfNum, paramsValues, B, PF, varargin)

[StimLevels, NumPos, OutOfNum] = PAL_PFML_GroupTrialsbyX(StimLevels, NumPos, OutOfNum);

options = [];
maxTries = 1;
rangeTries = [1 1 1 .1];
gammaEQlambda = logical(false);
lapseFit = 'nAPLE';

lapseLimits = [];
guessLimits = [];

converged = false(B,1);
DevSim = zeros(B,1);

MC = PAL_PFLR_setupMC(paramsValues);

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'SearchOptions',7)
            options = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'maxTries',4)
            maxTries = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'RangeTries',6)
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
        if strncmpi(varargin{n}, 'Thresh',6)
            MC.argsAlesser = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'Slopes',6)
            MC.argsBlesser = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'GuessR',6)
            MC.argsGlesser = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'LapseR',6)
            MC.argsLlesser = varargin{n+1};
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

if (strncmpi(lapseFit,'iap',3) || strncmpi(lapseFit,'jap',3)) && ((PAL_whatIs(MC.argsLlesser) == 2 && strncmpi(MC.argsLlesser,'fix',3)) || isempty(MC.argsLlesser))
    warning('PALAMEDES:invalidOption','Lapses rates are not free: ''LapseFit'' argument ignored');
    lapseFit = 'nap';
end
if (strncmpi(lapseFit,'iap',3) || strncmpi(lapseFit,'jap',3)) && PAL_whatIs(MC.argsLlesser) == 4
    warning('PALAMEDES:invalidOption','Lapse rates custom-reparameterized: ''lapseFit'' argument ignored');
    lapseFit = 'nap';
end

if gammaEQlambda
    MC.argsGlesser = [];
    if ~isempty(guessLimits)
        warning('PALAMEDES:invalidOption','Guess rates constrained to equal lapse rates: ''guessLimits'' ignored');
        guessLimits = [];
    end        
end

Dev = PAL_PFML_DevianceGoF(StimLevels,NumPos,OutOfNum,paramsValues,PF,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
 
for b = 1:B
    
    for cond = 1:size(StimLevels,1)
        NumPos(cond,:) = PAL_PF_SimulateObserverParametric(paramsValues(cond,:), StimLevels(cond,:), OutOfNum(cond,:), PF,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
    end
    [paramsValuesSim, LL, converged(b)] = PAL_PFML_FitMultiple(StimLevels,NumPos, OutOfNum, paramsValues, PF, 'Thresh', MC.argsAlesser, 'Slopes', MC.argsBlesser,'GuessR',MC.argsGlesser,'LapseR',MC.argsLlesser,'searchoptions',options, 'lapseLimits',lapseLimits,'guessLimits',guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);

    tries = 1;
    
    while converged(b) == 0 && tries < maxTries
        paramsTry = paramsValues;
        ArgsATry = MC.argsAlesser;
        ArgsBTry = MC.argsBlesser;
        ArgsGTry = MC.argsGlesser;
        ArgsLTry = MC.argsLlesser;
        if ~isstruct(rangeTries)
            multiplier = PAL_PFML_rangeTries(MC.argsAlesser,MC.argsBlesser,MC.argsGlesser,MC.argsLlesser, repmat(rangeTries,[size(paramsValues,1) 1]));
            paramsTry = paramsValues + multiplier.*repmat(rangeTries,size(paramsValues,1),1);
        else
            if isfield(rangeTries,'rangeTries')
                multiplier = PAL_PFML_rangeTries(MC.argsAlesser,MC.argsBlesser,MC.argsGlesser,MC.argsLlesser, repmat(rangeTries.rangeTries,[size(paramsValues,1) 1]));
                paramsTry = paramsValues + multiplier.*repmat(rangeTries.rangeTries,size(paramsValues,1),1);
            end
            if isstruct(ArgsATry) && isfield(rangeTries,'paramsValuesA');
                ArgsATry.paramsValuesA = ArgsATry.paramsValuesA + ArgsATry.paramsFreeA.*(rand(1,length(rangeTries.paramsValuesA))-.5).*rangeTries.paramsValuesA;
            end
            if isstruct(ArgsBTry) && isfield(rangeTries,'paramsValuesB');
                ArgsBTry.paramsValuesB = ArgsBTry.paramsValuesB + ArgsBTry.paramsFreeB.*(rand(1,length(rangeTries.paramsValuesB))-.5).*rangeTries.paramsValuesB;
            end
            if isstruct(ArgsGTry) && isfield(rangeTries,'paramsValuesG');
                ArgsGTry.paramsValuesG = ArgsGTry.paramsValuesG + ArgsGTry.paramsFreeG.*(rand(1,length(rangeTries.paramsValuesG))-.5).*rangeTries.paramsValuesG;
            end
            if isstruct(ArgsLTry) && isfield(rangeTries,'paramsValuesL');
                ArgsLTry.paramsValuesL = ArgsLTry.paramsValuesL + ArgsLTry.paramsFreeL.*(rand(1,length(rangeTries.paramsValuesL))-.5).*rangeTries.paramsValuesL;
            end
        end
        [paramsValuesSim, LL, converged(b)] = PAL_PFML_FitMultiple(StimLevels,NumPos, OutOfNum, paramsTry, PF, 'Thresh', ArgsATry, 'Slopes', ArgsBTry,'GuessR',ArgsGTry,'LapseR',ArgsLTry,'searchoptions',options, 'lapseLimits',lapseLimits,'guessLimits',guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);      
        tries = tries + 1;
    end
    if ~converged(b)
        warning('PALAMEDES:convergeFail','Fit to simulation %s of %s did not converge.',int2str(b), int2str(B));
    end
    
    DevSim(b) = PAL_PFML_DevianceGoF(StimLevels,NumPos,OutOfNum,paramsValuesSim,PF,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
end

exitflag = sum(converged) == B;
if exitflag ~= 1
    warning('PALAMEDES:convergeFail','Only %s of %s simulations converged.',int2str(sum(converged)), int2str(B));
end

pDev = length(DevSim(DevSim>Dev))/B;