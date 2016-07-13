%
%PAL_PFML_negLL     (negative) Log Likelihood associated with fit of
%   Psychometric Function
%
%Syntax: negLL = PAL_PFML_negLL(paramsFreeVals, paramsFixedVals, ...
%   paramsFree, StimLevels, NumPos, OutOfNum, PF, {optional arguments})
%
%Internal Function
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.3.0, 1.3.1, 1.4.0, 1.4.1, 1.4.2, 1.6.3
%   (see History.m)
%
% 5/27/2015  dhb, ar  Added code to check for and prevent solutions that
%                     result in complex predicted pcorrect.  See comment by inserted code
%                     below.

function negLL = PAL_PFML_negLL(paramsFreeVals, paramsFixedVals, paramsFree, StimLevels, NumPos, OutOfNum, PF, varargin)

lapseLimits = [];
guessLimits = [];
lapseFit = 'default';
gammaEQlambda = logical(false);

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'lapseLimits',6)
            lapseLimits = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'guessLimits',6)
            guessLimits = varargin{n+1};
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

params(paramsFree == 1) = paramsFreeVals;
params(paramsFree == 0) = paramsFixedVals;    

if gammaEQlambda
    params(3) = params(4);
end

pcorrect = PF(params, StimLevels);

if (~isempty(lapseLimits) && (params(4) < lapseLimits(1) || params(4) > lapseLimits(2))) || (~isempty(guessLimits) && (params(3) < guessLimits(1) || params(3) > guessLimits(2)))
    negLL = Inf;
    
%% This line added by DHB to prevent return of crazy parameters.
% We hit this when fitting a Weibull and found that the routine
% was returning a negative value for alpha and complex values
% for the returned probability correct.  We check for complex
% probability correct and push the error function to infinity
% when they occur.  This seems to work.
elseif (any(imag(pcorrect) ~= 0))
    negLL = Inf;
    
else    
    switch lower(lapseFit(1:3))
        case {'nap', 'def'}            
            negLL = sum(PAL_nansum(NumPos.*log(pcorrect)+(OutOfNum-NumPos).*log(1 - pcorrect)));
        case {'jap', 'iap'}        %F assumed to equal unity at highest stimulus level.            
            len = length(StimLevels);                   
            negLL = PAL_nansum(NumPos(len).*log(1-params(4)) + (OutOfNum(len)-NumPos(len)).*log(params(4)));
            if gammaEQlambda
                negLL = negLL + PAL_nansum((OutOfNum(1)-NumPos(1))*log(1-params(4)) + NumPos(1).*log(params(4)));
            end
            negLL = negLL + sum(PAL_nansum(NumPos(1+gammaEQlambda:len-1).*log(pcorrect(1+gammaEQlambda:len-1))+(OutOfNum(1+gammaEQlambda:len-1)-NumPos(1+gammaEQlambda:len-1)).*log(1 - pcorrect(1+gammaEQlambda:len-1))));
    end    
    negLL = -negLL;
end