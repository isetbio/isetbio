%
%PAL_PFML_BruteForceFit   Fit PF using a brute-force search through
%   4D parameter space.
%
%   syntax: [paramsValues LL LLspace] = ...
%       PAL_PFML_BruteForceFit(StimLevels, NumPos, OutOfNum, ...
%       searchGrid, PF)
%
%Input:
%   'StimLevels': vector containing stimulus levels used.
%
%   'NumPos': vector containing for each of the entries of 'StimLevels' the 
%       number of trials a positive response (e.g., 'yes' or 'correct') was
%       given.
%
%   'OutOfNum': vector containing for each of the entries of 'StimLevels' 
%       the total number of trials.
%
%   'searchGrid': structure defining the 4D parameter grid. Structure has
%       four fields: .alpha, .beta, .gamma, .lambda. Each specifies a
%       vector of values for the respective parameter to be contained in 
%       the grid. Use of a scalar fixes parameter to value of scalar. Note 
%       that choices made here have a large effect on processing time and 
%       memory usage.
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
%       parameters of the psychometric function. Note that this is a crude
%       fit only (crudeness is determined by grain of 'searchGrid'). For 
%       more precise fit, use PAL_PFML_Fit.
%
%   'LL': Log likelihood value associated with the fit.
%
%   'LLspace': Log Likelihood values for entire search grid.
%
%   If highest entry in 'StimLevels' is so high that it can be 
%   assumed that errors observed there can be due only to lapses, use 
%   'lapseFit' argument to specify alternative fitting scheme. Options: 
%   'nAPLE' (default), 'iAPLE', and 'jAPLE'. Type help 
%   PAl_PFML_FitMultiple for more information.
%
%   The guess rate and lapse rate parameter can be constrained to be 
%   equal, as would be appropriate, for example, in a bistable percept 
%   task. To accomplish this, use optional argument 'gammaEQlambda', 
%   followed by a 1. Entry for guess rate in searchGrid needs to be made 
%   but will be ignored.
%
%Full example:
%
%   PF = @PAL_Logistic;
%   StimLevels = [-3:1:3];
%   NumPos = [55 55 66 75 91 94 97];    %observer data
%   OutOfNum = 100.*ones(size(StimLevels));
%
%   searchGrid.alpha = [-1:.1:1];
%   searchGrid.beta = 10.^[-1:.1:2];
%   searchGrid.gamma = .5;
%   searchGrid.lambda = [0:.005:.06];
%
%   %Fit data:
%
%   [paramsValues LL] = PAL_PFML_BruteForceFit(StimLevels, NumPos, ...
%       OutOfNum, searchGrid, PF)
%
% Introduced: Palamedes version 1.2.0 (NP)
% Modified: Palamedes version 1.3.0, 1.3.1, 1.6.0, 1.6.3 (see History.m)

function [ paramsValues, maxim, LLspace] = PAL_PFML_BruteForceFit(StimLevels,NumPos, OutOfNum, searchGrid, PF, varargin)

    lapseFit = 'default';
    gammaEQlambda = logical(false);

    if ~isempty(varargin)
        NumOpts = length(varargin);
        for n = 1:2:NumOpts
            valid = 0;
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
    
    if gammaEQlambda
        searchGrid.gamma = 0;
    end
    
    [StimLevels, NumPos, OutOfNum] = PAL_PFML_GroupTrialsbyX(StimLevels, NumPos, OutOfNum);    

    switch lower(lapseFit(1:3))
        case {'nap', 'def'}
            [paramsGrid.alpha, paramsGrid.beta, paramsGrid.gamma, paramsGrid.lambda] = ndgrid(searchGrid.alpha,searchGrid.beta,searchGrid.gamma,searchGrid.lambda);
            if gammaEQlambda
                paramsGrid.gamma = paramsGrid.lambda;
            end            
            LLspace = zeros(size(paramsGrid.alpha,1),size(paramsGrid.alpha,2),size(paramsGrid.alpha,3),size(paramsGrid.alpha,4));
            for level = 1:length(StimLevels)
               LLspace = LLspace + NumPos(level).*log(PF(paramsGrid,StimLevels(level)))+(OutOfNum(level)-NumPos(level)).*log(1-PF(paramsGrid,StimLevels(level)));
            end
        case 'jap' 
            [paramsGrid.alpha, paramsGrid.beta, paramsGrid.gamma, paramsGrid.lambda] = ndgrid(searchGrid.alpha,searchGrid.beta,searchGrid.gamma,searchGrid.lambda);
            len = length(NumPos);            
            LLspace = zeros(size(paramsGrid.alpha,1),size(paramsGrid.alpha,2),size(paramsGrid.alpha,3),size(paramsGrid.alpha,4));
            if gammaEQlambda
                paramsGrid.gamma = paramsGrid.lambda;
                LLspace = log(paramsGrid.gamma.^NumPos(1))+log((1-paramsGrid.gamma).^(OutOfNum(1)-NumPos(1)));  %Unlike 0.*log(0), Matlab evaluates log(0.^0) as 0
            end
            LLspace = LLspace + log((1-paramsGrid.lambda).^NumPos(len))+log(paramsGrid.lambda.^(OutOfNum(len)-NumPos(len)));
            for level = 1+gammaEQlambda:len-1
                LLspace = LLspace + NumPos(level).*log(PF(paramsGrid,StimLevels(level)))+(OutOfNum(level)-NumPos(level)).*log(1-PF(paramsGrid,StimLevels(level)));
            end    
        case 'iap' 
            len = length(NumPos);
            if gammaEQlambda
                searchGrid.lambda = (OutOfNum(len)-NumPos(len)  +   NumPos(1))/(OutOfNum(len)+OutOfNum(1));
                searchGrid.gamma = searchGrid.lambda;
            else
                searchGrid.lambda = 1 - NumPos(len)/OutOfNum(len);
            end
            [paramsGrid.alpha, paramsGrid.beta, paramsGrid.gamma, paramsGrid.lambda] = ndgrid(searchGrid.alpha,searchGrid.beta,searchGrid.gamma,searchGrid.lambda);                                    
            LLspace = zeros(size(paramsGrid.alpha,1),size(paramsGrid.alpha,2),size(paramsGrid.alpha,3),size(paramsGrid.alpha,4));
            LLspace = LLspace + log((1-searchGrid.lambda).^NumPos(len)) + log(searchGrid.lambda.^(OutOfNum(len)-NumPos(len)));
            if gammaEQlambda
                LLspace = LLspace + log((1-searchGrid.gamma).^(OutOfNum(1)-NumPos(1))) + log(searchGrid.gamma.^NumPos(1));
            end
            for level = 1+gammaEQlambda:len-1
               LLspace = LLspace + NumPos(level).*log(PF(paramsGrid,StimLevels(level)))+(OutOfNum(level)-NumPos(level)).*log(1-PF(paramsGrid,StimLevels(level)));
            end
    end
        
    [maxim, I] = PAL_findMax(LLspace);
    
    while length(I) < 4
        I = [I 1];
    end
    
    paramsValues = [searchGrid.alpha(I(1)) searchGrid.beta(I(2)) searchGrid.gamma(I(3)) searchGrid.lambda(I(4))];

    if gammaEQlambda
        paramsValues(3) = paramsValues(4);
    end
    
end