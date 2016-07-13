%
%PAL_MLDS_Fit   Fit scaling data using the method developed by Maloney & 
%   Yang (2003) Journal of Vision, 3, 573-585.
%
%syntax: [PsiValues SDnoise LL exitflag output] = PAL_MLDS_Fit(Stim, ...
%   NumGreater, OutOfNum, PsiValues, SDnoise, {optional arguments})
%
%Input:
%   'Stim': Listing of stimulus pairs, triads or quadruples used in
%       experiment. May contain repeats (i.e., like trials need not be
%       grouped). See PAL_MLDS_GenerateStimList for format.
%
%   'NumGreater': For each row of 'Stim', 'NumGreater' lists the number of
%       trials on which the response was 'greater' (for pairs i-j: j was
%       judged greater than i, for triads i-j-k: (k - j) was deemed greater
%       than (j - i), for quadruples i-j-k-l: (k - l) was deemed greater
%       than (j - i)).
%
%   'OutOfNum': For each row of 'Stim', 'OutOfNum' lists the number of
%       trials.
%
%   'PsiValues' is a vector containing initial guesses for the values of 
%       Psi. Must have as many elements as there are stimulus levels in
%       'Stim', first entry must be 0, last entry must be 1.
%
%   'SDnoise' is a scalar containing initial guess for magnitude of
%       internal noise.
%
%Output:
%   'PsiValues': Best-fitting values for Psi.
%
%   'SDnoise': Best-fitting value for internal noise.
%
%   'LL': Log likelihood associated with the fit.
%
%   'exitflag': 1 indicates a succesful fit, 0 indicates fit did not
%       converge (trying again using new initial guesses might help).
%
%   'output': message containing some information concerning fitting
%       process.
%
%   PAL_MLDS_Fit uses Nelder-Mead Simplex method. The default search 
%       options may be changed by using the following syntax:
%
%   [PsiValues SDnoise LL exitflag output] = PAL_MLDS_Fit(Stim, ...
%       NumGreater, OutOfNum, PsiValues, SDnoise, 'SearchOptions', options)
%
%   where 'options' is a structure that can be created using:
%       options = PAL_minimize('options');
%   type PAL_minimize('options','help'); to get a brief explanation of
%       options available and their default values.
%
%Example:
%
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
%   options = PAL_minimize('options');
%   options.TolX = 1e-9;    %increase desired precision
%   [PsiValues SDnoise] = PAL_MLDS_Fit(Stim, NumGreater, OutOfNum, ...
%       PsiValues, SDnoise,'SearchOptions',options);
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.4.0, 1.6.3 (NP): (see History.m)

function [PsiValues, SDnoise, LL, exitflag, output] = PAL_MLDS_Fit(Stim, NumGreater, OutOfNum, PsiValues, SDnoise, varargin)

options = [];

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'SearchOptions',7)
            options = varargin{n+1};
            valid = 1;
        end
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});
        end        
    end            
end

NumLevels = length(PsiValues);

[FreeParams, negLL, exitflag, output] = PAL_minimize(@PAL_MLDS_negLL,[PsiValues(2:NumLevels-1) SDnoise], options, Stim, NumGreater, OutOfNum);

LL = -negLL;

PsiValues(2:NumLevels-1) = FreeParams(1:NumLevels-2);
SDnoise = FreeParams(NumLevels-1);