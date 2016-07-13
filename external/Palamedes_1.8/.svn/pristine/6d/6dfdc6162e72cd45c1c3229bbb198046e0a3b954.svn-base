%
%PAL_PF_SimulateObserverParametric   Simulate observer characterized by
%   psychometric function.
%
%syntax: NumPos = PAL_PF_SimulateObserverParametric(paramsValues, ...
%   StimLevels, OutOfNum, PF)
%
%Input:
%   'paramsValues': vector containing parameter values [alpha beta gamma 
%       lambda] of the psychometric function to be simulated. Typically,
%       these would be those characterizing human observer and obtained by
%       fitting a psychometric function by using a function such as 
%       PAL_PFML_Fit.
%
%   'StimLevels': vector containing stimulus levels.
%
%   'OutOfNum': Number of trials to be simulated for each of the entries in
%       'StimLevels'.
%
%   'PF': Form of the psychometric function to be simulated. Should be
%       passed as inline function. Options include:
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
%   'NumPos': simulated number of positive responses for each entry in
%       'StimLevels'.
%
%   PAL_PF_SimulateObserverParametric will accept a few optional
%   arguments:
%
%   'lapseFit': By default, 'lapseFit' will be 'nAPLE' (non-Asymptotic 
%       Performance Lapse Estimation) in which probability of positive 
%       response will be as in the standard PF equation:
%   
%       psi(x; a b g l) = g + (1 - g - l)F(x; a b), with
%
%       x: stimulus intensity, a: alpha, b: beta, g: gamma, l: lambda, and 
%       F a sigmoid with range [0 1]. 
%   
%       'lapseFit' may be changed to 'iAPLE' (isolated Asymptotic 
%       Performance Lapse Estimation) or 'jAPLE' (joint APLE), in which 
%       probability of positive response at the highest stimulus intensity 
%       equals 1 - the lapse rate (i.e., F is assumed to equal unity 
%       there). 'iAPLE' and 'jAPLE' act equivalently in this function.
%
%   'gammaEQlambda': In case 'lapseFit' is set to 'iAPLE' or 'jAPLE' and a
%       task is used in which the guess rate parameter and the lapse rate
%       parameter can be assumed to be equal (e.g., bistable percept),
%       'gammEQlambda' should be set to 1 (or logical TRUE). Probability
%       correct at lowest stimulus intensity will equal the lapse rate and
%       the value of paramsValues(3) (i.e., gamma or guess rate) will be 
%       ignored.
%
%Examples: 
%   NumPos = PAL_PF_SimulateObserverParametric([0 1 .5 0.05], ...
%    [-2:1:2],[100 100 100 100 100], @PAL_Logistic) might return:
%
%   NumPos =
%       
%       59  60  77  87  93
%
%   NumPos = PAL_PF_SimulateObserverParametric([0 1 .5 0.05], ...
%    [-2:1:2],[100 100 100 100 100], @PAL_Logistic,'lapseFit','iAPLE',...
%   'gammaEQlambda',logical(true)) might return:
%
%   NumPos =
%       
%       3   29  50  76  95
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.3.0, 1.4.0, 1.6.3 (see History.m)

function NumPos = PAL_PF_SimulateObserverParametric(paramsValues, StimLevels, OutOfNum, PF, varargin)

lapseFit = 'nAPLE';
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
   paramsValues(3) = paramsValues(4);
end

pcorrect = PF(paramsValues, StimLevels);

if strncmpi(lapseFit,'jap',3) || strncmpi(lapseFit,'iap',3)

    pcorrect(find(StimLevels == max(StimLevels(OutOfNum>0)))) = 1-paramsValues(4);
    if gammaEQlambda
        pcorrect(find(StimLevels == min(StimLevels(OutOfNum>0)))) = paramsValues(4);        
    end
end
NumPos = zeros(1,length(StimLevels));
for Level = 1:length(StimLevels)
    Pos = rand(OutOfNum(Level),1);
    Pos(Pos < pcorrect(Level)) = 1;
    Pos(Pos ~= 1) = 0;
    NumPos(Level) = sum(Pos);
end