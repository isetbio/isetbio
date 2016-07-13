%
%PAL_PFML_BootstrapParametricMultiple    Perform parametric bootstrap to
%   determine standard errors on parameters of psychometric functions
%   fitted simultaneously to several conditions.
%
%syntax: [SD paramsSim LLSim converged SDfunc funcParamsSim] = ...
%   PAL_PFML_BootstrapParametricMultiple(StimLevels, OutOfNum, ...
%       paramsValues, B, PF,{optional arguments})
%
%Input: 
%   'StimLevels': matrix containing stimulus levels used. Each row of the
%       matrix corresponds to a condition.
%
%   'OutOfNum': matrix containing for each of the entries of 'StimLevels' 
%       the total number of trials. (May contain zeros in case not all 
%       conditions utilized identical numbers of stimulus levels).
%
%   'paramsValues': Matrix containing parametervalues to be used during
%       simulations. Each row should correspond to a condition. Matrix 
%       should have four columns [threshold slope guess-rate lapse-rate]. 
%       PAL_PFML_FitMultiple might be used to obtain parameter values 
%       characterizing observer's psychometric functions.
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
%   'B': number of bootstrap simulations to perform.
%
%
%Output:
%   'SD': Matrix containing standard deviations of the PF's parameters
%       across the B fits to simulated data. Rows correspond to different
%       conditions. These are estimates of the standard errors of the 
%       parameter estimates.
%
%   'paramsSim': Bx[number of conditions]x4 matrix containing the fitted 
%       parameters for all B fits to simulated data.
%
%   'LLSim': vector containing Log Likelihoods associated with all B fits
%       to simulated data.
%
%   'converged': For each simulation contains a 1 in case the fit was
%       succesfull (i.e., converged) or a 0 in case it did not.
%
%   'SDfunc': structure containing standard errors of parameter estimates 
%       of custom-defined reparametrizations (if any). Type 'help
%       PAL_PFML_CustomDefine' for more information.
%
%   'funcParamsSim': structure containing best-fitting parameter estimates 
%       of parameters in the reparametrization structure for each of the 
%       simulations. Type 'help PAL_PFML_CustomDefine' for more 
%       information.
%
%Unless specified otherwise, PAL_PFML_BootstrapParametricMultiple will fit 
%   all conditions separately treating threshold and slope as free 
%   parameters and the guess-rate and lapse-rate as fixed parameters. 
%   However, other fitting schemes may be specified by user. Refer to 
%   PAL_PFML_FitMultiple or example below for information on how to specify 
%   other models.
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
%   PAL_PFML_BootstrapParametricMultiple will generate a warning if not all 
%       simulations were fit succesfully.
%
%   PAL_PFML_BootstrapParametricMultiple will accept a few other optional 
%       arguments:
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
%   PAL_PFML_BootstrapParametricMultiple uses Nelder-Mead Simplex method to 
%   find the maximum in the likelihood function. The default search options 
%   may be changed by using the optional argument 'searchOptions' followed 
%   by an options structure created using options = PAL_minimize('options'). 
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
%   paramsValues = [0 1.5 .5 .04]; %Initial guesses
%
% %Fitted using defaults:
%
%   paramsFitted = PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, ...
%       paramsValues, PF,'searchOptions',options);
%
%   %Use fitted model in bootstrap simulations:
%
%   B = 400;
%   SD = PAL_PFML_BootstrapParametricMultiple(StimLevels, OutOfNum, ...
%       paramsFitted, B, PF,'searchOptions',options)
%
%might return:
%
%   SD =
% 
%     0.1156    0.1921         0    0.0000
%     0.0860    0.2147         0    0.0000
%     0.1052    0.2066         0    0.0000
%     0.0985    0.1558         0    0.0000
%
% %Fitted using different model:
%
%   paramsFitted = PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, ...
%       paramsValues, PF,'lapserates',[1 1 1 1;1 1 -1 -1],'lapsefit',...
%       'isp','searchOptions',options);
%
% %Bootstrap below might take a while to complete:
%
%   SD = PAL_PFML_BootstrapParametricMultiple(StimLevels, OutOfNum, ...
%       paramsFitted, B, PF,'lapserates',[1 1 1 1; 1 1 -1 -1],'lapseFit',...
%       'isp' ,'maxTries',4, 'rangeTries',[2 1 0 0],'searchOptions',options)
%
%Might return:
%
%   SD =
% 
%     0.1145    0.9843         0    0.0095
%     0.1033    0.3110         0    0.0095
%     0.1150    0.2813         0    0.0088
%     0.1063    0.2367         0    0.0088
%    
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.0.2, 1.1.0, 1.2.0, 1.3.0, 1.4.0, 1.4.2, 
%   1.6.0, 1.6.3 (see History.m)

function [SD, paramsSim, LLSim, converged, SDfunc, funcParamsSim] = PAL_PFML_BootstrapParametricMultiple(StimLevels, OutOfNum, paramsValues, B, PF, varargin)

[StimLevels, trash, OutOfNum] = PAL_PFML_GroupTrialsbyX(StimLevels, OutOfNum, OutOfNum);

options = [];
maxTries = 1;
rangeTries = [1 1 1 .1];
lapseFit = 'nAPLE';
gammaEQlambda = logical(false);
lapseLimits = [];
guessLimits = [];

SDfunc = [];
paramsSim = zeros(B,size(paramsValues,1),size(paramsValues,2));
funcParamsSim = PAL_PFML_setupParametrizationStruct;
LLSim = zeros(B,1); 
converged = false(B,1);

ArgsA = PAL_Contrasts(size(paramsValues,1),'identity');
ArgsB = PAL_Contrasts(size(paramsValues,1),'identity');
ArgsG = [];
ArgsL = [];

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
            ArgsA = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'Slopes',6)
            ArgsB = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'GuessR',6)
            ArgsG = varargin{n+1};
            valid = 1;            
        end
        if strncmpi(varargin{n}, 'LapseR',6)
            ArgsL = varargin{n+1};
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

if (strncmpi(lapseFit,'iap',3) || strncmpi(lapseFit,'jap',3)) && ((PAL_whatIs(ArgsL) == 2 && strncmpi(ArgsL,'fix',3)) || isempty(ArgsL))
    warning('PALAMEDES:invalidOption','Lapses rates are not free: ''LapseLimits'' argument ignored');
    lapseFit = 'nap';
end
if (strncmpi(lapseFit,'iap',3) || strncmpi(lapseFit,'jap',3)) && PAL_whatIs(ArgsL) == 4
    warning('PALAMEDES:invalidOption','Lapse rates custom-reparameterized: ''lapseFit'' argument ignored');
    lapseFit = 'nap';
end

if gammaEQlambda
    ArgsG = [];
    if ~isempty(guessLimits)
        warning('PALAMEDES:invalidOption','Guess rates constrained to equal lapse rates: ''guessLimits'' ignored');
        guessLimits = [];
    end    
end

SumExitFlags = 0;

for b = 1:B

    %Simulate experiment
    for cond = 1:size(StimLevels,1)
        NumPosSim(cond,:) = PAL_PF_SimulateObserverParametric(paramsValues(cond,:), StimLevels(cond,:), OutOfNum(cond,:), PF,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
    end

    [paramsSim(b,:,:), LLSim(b), converged(b), trash, temp] = PAL_PFML_FitMultiple(StimLevels, NumPosSim, OutOfNum, paramsValues, PF, 'Thresh', ArgsA, 'Slopes', ArgsB,'GuessR',ArgsG,'LapseR',ArgsL,'searchoptions',options, 'lapseLimits', lapseLimits,'guessLimits', guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);

    if isstruct(ArgsA) || isstruct(ArgsB) || isstruct(ArgsG) || isstruct(ArgsL)
        funcParamsSim(b) = temp;
    end
    tries = 1;

    while converged(b) == 0 && tries < maxTries
        paramsTry = paramsValues;
        ArgsATry = ArgsA;
        ArgsBTry = ArgsB;
        ArgsGTry = ArgsG;
        ArgsLTry = ArgsL;
        if ~isstruct(rangeTries)
            multiplier = PAL_PFML_rangeTries(ArgsA, ArgsB, ArgsG, ArgsL, repmat(rangeTries,size(paramsValues,1),1));
            paramsTry = paramsValues + multiplier.*repmat(rangeTries,size(paramsValues,1),1);
        else
            if isfield(rangeTries,'rangeTries')
                multiplier = PAL_PFML_rangeTries(ArgsA, ArgsB, ArgsG, ArgsL, repmat(rangeTries.rangeTries,size(paramsValues,1),1));
                paramsTry = paramsValues + multiplier.*repmat(rangeTries.rangeTries,size(paramsValues,1),1);
            end
            if isstruct(ArgsA) && isfield(rangeTries,'paramsValuesA');
                ArgsATry.paramsValuesA = ArgsATry.paramsValuesA + ArgsATry.paramsFreeA.*(rand(1,length(rangeTries.paramsValuesA))-.5).*rangeTries.paramsValuesA;
            end
            if isstruct(ArgsB) && isfield(rangeTries,'paramsValuesB');
                ArgsBTry.paramsValuesB = ArgsBTry.paramsValuesB + ArgsBTry.paramsFreeB.*(rand(1,length(rangeTries.paramsValuesB))-.5).*rangeTries.paramsValuesB;
            end
            if isstruct(ArgsG) && isfield(rangeTries,'paramsValuesG');
                ArgsGTry.paramsValuesG = ArgsGTry.paramsValuesG + ArgsGTry.paramsFreeG.*(rand(1,length(rangeTries.paramsValuesG))-.5).*rangeTries.paramsValuesG;
            end
            if isstruct(ArgsL) && isfield(rangeTries,'paramsValuesL');
                ArgsLTry.paramsValuesL = ArgsLTry.paramsValuesL + ArgsLTry.paramsFreeL.*(rand(1,length(rangeTries.paramsValuesL))-.5).*rangeTries.paramsValuesL;
            end
        end
        [paramsSim(b,:,:), LLSim(b), converged(b), trash, temp] = PAL_PFML_FitMultiple(StimLevels, NumPosSim, OutOfNum, paramsTry, PF, 'Thresh', ArgsATry, 'Slopes', ArgsBTry,'GuessR',ArgsGTry,'LapseR',ArgsLTry,'searchoptions',options, 'lapseLimits', lapseLimits,'guessLimits', guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);

        if isstruct(ArgsA) || isstruct(ArgsB) || isstruct(ArgsG) || isstruct(ArgsL)
            funcParamsSim(b) = temp;
        end
        tries = tries + 1;
    end    
    if ~converged(b)
        warning('PALAMEDES:convergeFail','Fit to simulation %s of %s did not converge.',int2str(b), int2str(B));
    end

end

exitflag = sum(converged) == B;
if exitflag ~= 1
    warning('PALAMEDES:convergeFail','Only %s of %s simulations converged.',int2str(sum(converged)), int2str(B));
end

for cond = 1:size(StimLevels,1)
    [Mean(cond,:), SD(cond,:)] = PAL_MeanSDSSandSE(squeeze(paramsSim(:,cond,:)));
end
if ~isempty(funcParamsSim(1).paramsValuesA)
    [Meanfunc.A, SDfunc.A] = PAL_MeanSDSSandSE(cat(1,funcParamsSim.paramsValuesA));
end
if ~isempty(funcParamsSim(1).paramsValuesB)
    [Meanfunc.B, SDfunc.B] = PAL_MeanSDSSandSE(cat(1,funcParamsSim.paramsValuesB));
end
if ~isempty(funcParamsSim(1).paramsValuesG)
    [Meanfunc.G, SDfunc.G] = PAL_MeanSDSSandSE(cat(1,funcParamsSim.paramsValuesG));
end
if ~isempty(funcParamsSim(1).paramsValuesL)
    [Meanfunc.L, SDfunc.L] = PAL_MeanSDSSandSE(cat(1,funcParamsSim.paramsValuesL));
end