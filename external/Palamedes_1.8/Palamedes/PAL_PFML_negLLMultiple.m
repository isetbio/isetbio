%
%PAL_PFML_negLLMultiple     (negative) Log Likelihood associated with 
%   simultaneous fit of multiple Psychometric Functions
%
%Syntax: negLL = PAL_negLLMultiple(thetas, thetasID, params, ...
%   StimLevels, NumPos, OutOfNum, FM, PF, {optional arguments})
%
%Internal Function
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.3.0, 1.4.0, 1.6.0, 1.6.3 (see History.m)

function negLL = PAL_PFML_negLLMultiple(thetas, thetasID, params, StimLevels, NumPos, OutOfNum, FM, PF, varargin)

lapseLimits = [];
guessLimits = [];
lapseFit = 'nAPLE';
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

params = PAL_PFML_TtoP(params, thetas, thetasID, FM);

if (~isempty(lapseLimits) && (min(params(:,4)) < lapseLimits(1) || max(params(:,4)) > lapseLimits(2))) || (~isempty(guessLimits) && (min(params(:,3)) < guessLimits(1) || max(params(:,3)) > guessLimits(2)))
    negLL = Inf;    
else
    for Cond = 1:size(StimLevels,1)
        SL = StimLevels(Cond,OutOfNum(Cond,:)~=0);
        NP = NumPos(Cond,OutOfNum(Cond,:)~=0);
        OoN = OutOfNum(Cond,OutOfNum(Cond,:)~=0);
        negLLcond(Cond) = PAL_PFML_negLL([params(Cond,1) params(Cond,2) params(Cond,3) params(Cond,4)], [], [1 1 1 1], SL, NP, OoN,PF,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
    end
    negLL = sum(negLLcond);
end