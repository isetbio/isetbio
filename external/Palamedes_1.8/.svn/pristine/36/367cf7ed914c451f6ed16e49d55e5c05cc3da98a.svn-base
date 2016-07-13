%
%PAL_PFML_Fit   Fit a psychometric function to data using a Maximum 
%    Likelihood criterion.
%
%Syntax: [paramsValues LL exitflag output] = PAL_PFML_Fit(StimLevels, ...
%   NumPos, OutOfNum, searchGrid, paramsFree, PF, {optional arguments})
%
%Input: 
%   'StimLevels': vector containing stimulus levels used.
%
%   'NumPos':
%       vector containing for each of the entries of 'StimLevels' the 
%       number of trials a positive response (e.g., 'yes' or 'correct') was
%       given.
%
%   'OutOfNum': vector containing for each of the entries of 'StimLevels' 
%       the total number of trials.
%
%   'searchGrid': Either a 1x4 vector [threshold slope guess-rate 
%       lapse-rate] containing initial guesses for free parametervalues and 
%       fixed values for fixed parameters or a structure with vector fields 
%       .alpha, .beta, .gamma, .lambda collectively defining a 4D parameter 
%       grid through which to perform a brute-force search for initial 
%       guesses (using PAL_PFML_BruteForceFit). These initial values will
%       serve as seeds for the Nelder-Mead iterative search. Fields for 
%       fixed parameters should be scalars equal to the fixed value. Note 
%       that fields for free parameters may also be scalars in which case 
%       the provided value will serve as the initial value for Nelder-Mead 
%       search. Note that choices made here have a large effect on 
%       processing time and memory usage.
%
%   'paramsFree': 1x4 vector coding which of the four parameters of the PF 
%       [threshold slope guess-rate lapse-rate] are free parameters and 
%       which are fixed parameters (1: free, 0: fixed).
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
%   'paramsValues': 1x4 vector containing values of fitted and fixed 
%       parameters of the psychometric function [threshold slope guess-rate 
%       lapse-rate].
%
%   'LL': Log likelihood associated with the fit.
%
%   'exitflag': 1 indicates a succesful fit, 0 indicates fit did not
%       converge (trying again using new initial guesses might help).
%
%   'output': message containing some information concerning fitting
%       process.
%
%User may constrain the lapse rate to fall within limited range using 
%   the optional argument 'lapseLimits', followed by two-element vector 
%   containing lower and upper limit respectively. See full example below.
%
%User may constrain the guess rate to fall within limited range using 
%   the optional argument 'guessLimits', followed by two-element vector 
%   containing lower and upper limit respectively. See full example below.
%
%Different fitting schemes may be specified using the optional argument
%   'lapseFits'. Default value is 'nAPLE' (non-Asymptotic Performance Lapse 
%   Estimation), in which all free parameters are estimated jointly in the 
%   manner outlined by e.g., Wichmann & Hill (2001). In case the highest 
%   stimulus value in 'StimLevels' is so high that it can be assumed that 
%   errors at this intensity can only be due to lapses, alternative schemes 
%   may be specified. 'jAPLE' (joint APLE) assumes that errors at
%   highest 'StimLevel' are due exclusively to lapses, and observations at 
%   this intensity contribute only to estimate of lapse rate. Observations 
%   at other values in 'StimLevels' contribute to all free parameters 
%   (including lapse rate). 'iAPLE' (isolated APLE) is identical to 
%   'jAPLE', except that the lapse rate estimate is based on observations 
%   made at the highest 'StimLevel' only. Other parameters are estimated 
%   from observations made at the other values of 'StimLevel' using lapse 
%   rate estimate obtained at highest 'StimLevel' as fixed value. For more
%   information see: http://www.journalofvision.org/content/12/6/25
%
%The guess rate and lapse rate parameter can be constrained to be equal, as 
%   would be appropriate, for example, in a bistable percept task. To 
%   accomplish this, use optional argument 'gammaEQlambda', followed by a 
%   1. Both the guess rate and lapse rate parameters will be fit according 
%   to options set for the lapse rate parameter. Entry for guess rate in 
%   'searchGrid' needs to be made but will be ignored.
%
%PAL_PFML_Fit uses Nelder-Mead Simplex method to find the maximum 
%   in the likelihood function. The default search options may be changed 
%   by using the optional argument 'SearchOptions' followed by an options 
%   structure created using options = PAL_minimize('options'). See example 
%   of usage below. For more information type 
%   PAL_minimize('options','help').
%
%   [paramsValues LL exitflag output] = PAL_PFML_Fit(StimLevels, NumPos, 
%       OutOfNum, searchGrid, paramsFree, PF, 'searchOptions', options)
%
%Full example:
%
%   options = PAL_minimize('options')
%   PF = @PAL_Logistic;
%   StimLevels = [-3:1:3];
%   NumPos = [0 13 28 56 73 91 93];    %observer data
%   OutOfNum = 100.*ones(size(StimLevels));
%   searchGrid.alpha = [-1:.01:1];    %structure defining grid to
%   searchGrid.beta = 10.^[-1:.01:2]; %search for initial values
%   searchGrid.gamma = [0:.01:.06];
%   searchGrid.lambda = [0:.01:.06];
%
%   %or (not advised):
%   % searchGrid = [0 1 0.01 0.01];       %Guesses
%
%   paramsFree = [1 1 1 1];
%
%   %Fit data:
%
%   paramsValues = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, ...
%       searchGrid, paramsFree, PF,'lapseLimits',[0 1],'guessLimits',...
%       [0 1], 'searchOptions',options)
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.0.2, 1.2.0, 1.3.0, 1.3.1, 1.4.0, 1.6.3
%   (see History.m)

function [paramsValues, LL, exitflag, output] = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, searchGrid, paramsFree, PF, varargin)

options = [];
lapseLimits = [];
guessLimits = [];
lapseFit = 'default';
gammaEQlambda = logical(false);

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'SearchOptions',7)
            options = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'lapseLimits',6)
            if paramsFree(4) == 0 && ~isempty(varargin{n+1})
                warning('PALAMEDES:invalidOption','Lapse rate is not a free parameter: ''LapseLimits'' argument ignored');
            else
                lapseLimits = varargin{n+1};
            end
            valid = 1;
        end
        if strncmpi(varargin{n}, 'guessLimits',6)
            if paramsFree(3) == 0 && ~isempty(varargin{n+1})
                warning('PALAMEDES:invalidOption','Guess rate is not a free parameter: ''GuessLimits'' argument ignored');
            else
                guessLimits = varargin{n+1};
            end
            valid = 1;
        end
        if strncmpi(varargin{n}, 'lapseFit',6)
            if paramsFree(4) == 0 && ~strncmpi(varargin{n+1},'def',3)
                warning('PALAMEDES:invalidOption','Lapse rate is not a free parameter: ''LapseFit'' argument ignored');
            else
                lapseFit = varargin{n+1};
            end
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

if ~isempty(guessLimits) && gammaEQlambda
    warning('PALAMEDES:invalidOption','Guess rate is constrained to equal lapse rate: ''guessLimits'' argument ignored');
    guessLimits = [];
end

[StimLevels, NumPos, OutOfNum] = PAL_PFML_GroupTrialsbyX(StimLevels, NumPos, OutOfNum);

if isstruct(searchGrid)
    if gammaEQlambda
        searchGrid.gamma = 0;
    end
    searchGrid = PAL_PFML_BruteForceFit(StimLevels, NumPos, OutOfNum, searchGrid, PF, 'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
end

if gammaEQlambda
    searchGrid(3) = searchGrid(4);
    paramsFree(3) = 0;
end

paramsFreeVals = searchGrid(paramsFree == 1);
paramsFixedVals = searchGrid(paramsFree == 0);

if isempty(paramsFreeVals)
    negLL = PAL_PFML_negLL(paramsFreeVals, paramsFixedVals, paramsFree, StimLevels, NumPos, OutOfNum, PF,'lapseFit',lapseFit);
    exitflag = 1;
    output = [];
else
    switch lower(lapseFit(1:3))
        case 'iap'
            len = length(NumPos);        
            if ~gammaEQlambda
                lambda = 1 - NumPos(len)./OutOfNum(len);
                if ~isempty(lapseLimits)
                    lambda = min(lambda,lapseLimits(2));
                    lambda = max(lambda,lapseLimits(1));
                end
                paramsFixedVals(length(paramsFixedVals)+1) = lambda; %lapse rate
                paramsFreeVals = paramsFreeVals(1:length(paramsFreeVals)-1);    %set lapse rate estimate as fixed value
                paramsFree(4) = 0;               
                [paramsFreeVals, negLL, exitflag, output] = PAL_minimize(@PAL_PFML_negLL, paramsFreeVals, options, paramsFixedVals, paramsFree, StimLevels(1:len-1), NumPos(1:len-1), OutOfNum(1:len-1), PF,'gammaEQlambda',gammaEQlambda,'guessLimits',guessLimits);        
                negLL = negLL - log((1 - lambda).^NumPos(len)) - log(lambda.^(OutOfNum(len)-NumPos(len)));
            else
                lambda = 1 - (NumPos(len)+(OutOfNum(1)-NumPos(1)))./(OutOfNum(len)+OutOfNum(1));                
                if ~isempty(lapseLimits)
                    lambda = min(lambda,lapseLimits(2));
                    lambda = max(lambda,lapseLimits(1));
                end
                paramsFixedVals(length(paramsFixedVals)+1) = lambda;
                paramsFixedVals(length(paramsFixedVals)) = lambda;
                paramsFreeVals = paramsFreeVals(1:length(paramsFreeVals)-1);
                paramsFree(4) = 0;
                [paramsFreeVals, negLL, exitflag, output] = PAL_minimize(@PAL_PFML_negLL, paramsFreeVals, options, paramsFixedVals, paramsFree, StimLevels(2:len-1), NumPos(2:len-1), OutOfNum(2:len-1), PF,'gammaEQlambda',gammaEQlambda);        
                negLL = negLL - log((1 - lambda).^(NumPos(len)+(OutOfNum(1)-NumPos(1)))) - log(lambda.^(OutOfNum(len)-NumPos(len)+NumPos(1)));
            end            
        case {'nap', 'jap', 'def'}
            [paramsFreeVals, negLL, exitflag, output] = PAL_minimize(@PAL_PFML_negLL, paramsFreeVals, options, paramsFixedVals, paramsFree, StimLevels, NumPos, OutOfNum, PF,'lapseLimits',lapseLimits,'guessLimits',guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
    end
end

paramsValues = zeros(1,4);
paramsValues(paramsFree == 1) = paramsFreeVals;
paramsValues(paramsFree == 0) = paramsFixedVals;

if gammaEQlambda
    paramsValues(3) = paramsValues(4);
end

LL = -1*negLL;