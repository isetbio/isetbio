%
%PAL_MLDS_Bootstrap     Estimate standard errors of parameters in MLDS
%   setting using Bootstrap procedure.
%
%syntax: [SD_PsiValues SD_SDnoise paramsSim LLSim converged] = 
%   PAL_MLDS_Bootstrap(Stim, OutOfNum, PsiValues, SDnoise, B, ...
%   {optional arguments})
%
%Input:
%   'Stim': Stimulus set (see e.g., PAL_MLDS_GenerateStimList)
%
%   'OutOfNum': Number of trials for each row in 'Stim'
%
%   'PsiValues': PsiValues characterizing observer (use PAL_MLDS_Fit to
%       find best-fitting PsiValues)
%
%   'SDnoise': internal noise characterizing observer (use PAL_MLDS_Fit to
%       find best-fitting SDnoise)
%
%   'B': Number of bootstrap simulations to perform.
%
%Output:
%   'SD_PsiValues': Standard deviations of fitted Psi value parameters 
%       across the B simulations. These are estimates of the standard 
%       errors of the Psi value parameters.
%
%   'SD_SDnoise': Standard deviation of fitted internal noise parameter
%       across the B simulations. This is an estimate of the standard error
%       of the SDnoise parameter.
%
%   'paramsSim': Listing of fitted parameter values across all B
%       simulations. May be useful e.g., to determine symmetry of
%       distribution, to find outliers, etc.
%
%   'LLSim': Listing of Log Likelihoods for fits to all B simulations. May
%       be useful e.g., to compare human observer fit to simulated fits.
%
%   'converged': For each simulation contains a 1 in case the fit was
%       succesfull (i.e., converged) or a 0 in case it did not.
%
%   PAL_MLDS_Bootstrap will generate a warning if not all simulations were
%       fit succesfully.
%
%   PAL_MLDS_Bootstrap will accept a few optional arguments:
%
%       In case not all fits converge succesfully, use optional argument 
%       'maxTries' to set the maximum number of times each fit will be 
%       attempted. The first try uses initial search values equal to 
%       'PsiValues' provided in function call, subsequent tries use these 
%       search values plus some random jitter. The range of the random 
%       jitter can be set using the optional argument 'rangeTries'. Default 
%       value for 'maxTries' is 1, default value for 'rangeTries' is .1 
%       (which applies to all parameter values). 
%       
%       Example: [SD_PsiValues SD_SDnoise paramsSim LLSim converged] = ...
%           PAL_MLDS_Bootstrap(Stim, OutOfNum, PsiValues, SDnoise, B, ...
%           'maxTries', 5, 'rangeTries',.2)
%
%           will try each fit up to five times using a random jitter range
%           equal to .2 (for each free parameter, initial search values
%           will be drawn from rectangular distribution with range .2
%           centered on user-provided PsiValues).
%
%       (Note that some simulated data sets may never be fit succesfully 
%           regardless of value of 'maxTries' and 'rangeTries') 
%
%       The default search parameters used by the Nelder-Mead simplex
%           search may be changed by using the optional argument
%           'SearchOptions' as in PAL_MLDS_Fit. Type help PAL_MLDS_Fit for
%           example of use.
%
%Example:
%   Stim = PAL_MLDS_GenerateStimList(2, 6, 2, 10);
%   OutOfNum = ones(1,size(Stim,1));
%   PsiValues = [0:1/5:1];
%   SDnoise = .5;
%
%   %Generate hypothetical data:
%   NumGreater = PAL_MLDS_SimulateObserver(Stim, OutOfNum, PsiValues, ...
%      SDnoise);
%  
%   %Fit hypothetical data:
%   [PsiValues SDnoise] = PAL_MLDS_Fit(Stim, NumGreater, OutOfNum, ...
%       PsiValues, SDnoise);
%
%   %Run bootstrap (using fitted parameter values):
%   [SD_PsiValues SD_SDnoise] = PAL_MLDS_Bootstrap(Stim, OutOfNum, ...
%       PsiValues, SDnoise, 100, 'maxTries',4,'rangeTries',.1);
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.1.0, 1.2.0, 1.4.0, 1.6.3 (see History.m)

function [SD_PsiValues, SD_SDnoise, paramsSim, LLSim, converged] = PAL_MLDS_Bootstrap(Stim, OutOfNum, PsiValues, SDnoise, B, varargin)

options = [];
maxTries = 1;
rangeTries = .1;
converged = false(B,1);
LLSim = zeros(B,1);
paramsSim = zeros(B,length(PsiValues)-1);
size(paramsSim);

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
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n})    
        end        
    end            
end

NumLevels = length(PsiValues);

for b = 1:B

    NumGreater = PAL_MLDS_SimulateObserver(Stim, OutOfNum, PsiValues, SDnoise);

    [paramsSim(b,:), LLSim(b), converged(b)] = PAL_minimize(@PAL_MLDS_negLL,[PsiValues(2:NumLevels-1) SDnoise], options, Stim, NumGreater, OutOfNum);

    tries = 1;
    while converged(b) == 0 && tries < maxTries        
        NewSearchInitials = [PsiValues(2:NumLevels-1) SDnoise]+(rand(1,length(PsiValues)-1)-.5).*rangeTries;
        [paramsSim(b,:), LLSim(b,:), converged(b)] = PAL_minimize(@PAL_MLDS_negLL,NewSearchInitials, options, Stim, NumGreater, OutOfNum);
        tries = tries + 1;
    end    
    if ~converged(b)
        warning('PALAMEDES:convergeFail','Fit to simulation %s of %s did not converge.',int2str(b), int2str(B));
    end    
end

paramsSim = [zeros(B,1) paramsSim(:,1:size(paramsSim,2)-1) ones(B,1) paramsSim(:,size(paramsSim,2))];

[Mean, SD] = PAL_MeanSDSSandSE(paramsSim);
SD_PsiValues = SD(1:length(SD)-1);
SD_SDnoise = SD(length(SD));
LLSim = -LLSim;
exitflag = sum(converged) == B;
if exitflag ~= 1
    warning('PALAMEDES:convergeFail','Only %s of %s simulations converged',int2str(sum(converged)),int2str(B));
end