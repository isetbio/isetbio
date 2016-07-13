%PAL_PFML_paramsTry     Generate jitter on values of guesses to be
%   supplied to PAL_PFML_Fit or PAL_PFML_FitMultiple as initial values in
%   search.
%
%Syntax: paramsTry = PAL_PFML_paramsTry(params, rangeTries, ...
%   {optional arguments});
%
%Input:
%   'params': n x 4 matrix (n may be 1) containing initial guesses for the 
%       4 parameters of PF in n conditions.
%   'rangeTries': matrix of same size as 'params' containing the range of
%       jitter to be applied to values in params. Jitter will be sampled 
%       from uniform distribution of width specified in rangeTries and
%       centered on values in params. If 'params' is n x 4 and rangeTries
%       is 1 x 4, the jitter rage specified in rangeTries will be applied 
%       to all rows in 'params'.
%
%Optional arguments may be used to constrain jitter to be constrained. For
%   example, if thresholds are constrained in the fit (see help in 
%   PAL_PFML_FitMultiple), one should constrain the jitter also. For each
%   of the four parameters of a PF, one may specify the following options:
%   'constrained', 'unconstrained','fixed'.
%
%Example: paramsTry = PAL_PFML_paramsTry([2 1 0 .06; 2 1 0 .06], ...
%   [2 1 0 .03],'thresholds','constrained','slopes','unconstrained',...
%   'lapserates','fixed') 
%   
%   might return:
%
%   paramsTry =
%
%    2.8315    1.2922         0    0.0600
%    2.8315    1.4595         0    0.0600
%
%Introduced: Palamedes version 1.1.1 (NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function [ paramsTry ] = PAL_PFML_paramsTry(params, rangeTries, varargin)

thresholds = 'unconstr';
slopes = 'unconstr';
guessrates = 'fixed';
lapserates = 'fixed';

paramsTry = zeros(size(params));

if size(rangeTries,1) == 1 && size(params,1) > 1
    rangeTries = repmat(rangeTries, [size(params,1) 1]);
end

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'thresholds',4) && (strncmpi(varargin{n+1}, 'con',3) || strncmpi(varargin{n+1}, 'unc',3) || strncmpi(varargin{n+1}, 'fix',3))
            thresholds = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'slopes',4) && (strncmpi(varargin{n+1}, 'con',3) || strncmpi(varargin{n+1}, 'unc',3) || strncmpi(varargin{n+1}, 'fix',3))
            slopes = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'guessrates',4) && (strncmpi(varargin{n+1}, 'con',3) || strncmpi(varargin{n+1}, 'unc',3) || strncmpi(varargin{n+1}, 'fix',3))
            guessrates = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'lapserates',4) && (strncmpi(varargin{n+1}, 'con',3) || strncmpi(varargin{n+1}, 'unc',3) || strncmpi(varargin{n+1}, 'fix',3))
            lapserates = varargin{n+1};
            valid = 1;
        end
        if ~valid
            warning('PALAMEDES:invalidOption','''%s'', ''%s'' is not a valid combination of options. Using default instead.',varargin{n},varargin{n+1});
        end
    end
end

switch lower(thresholds(1:3))
    case 'con'
        paramsTry(:,1) = ones(size(paramsTry(:,1)))*rand(1)-0.5;
    case 'unc'
        paramsTry(:,1) = rand(size(paramsTry(:,1)))-0.5;
    case 'fix'
        paramsTry(:,1) = zeros(size(paramsTry(:,1)));
end

switch lower(slopes(1:3))
    case 'con'
        paramsTry(:,2) = ones(size(paramsTry(:,1)))*rand(1)-0.5;
    case 'unc'
        paramsTry(:,2) = rand(size(paramsTry(:,1)))-0.5;
    case 'fix'
        paramsTry(:,2) = zeros(size(paramsTry(:,1)));
end

switch lower(guessrates(1:3))
    case 'con'
        paramsTry(:,3) = ones(size(paramsTry(:,1)))*rand(1)-0.5;
    case 'unc'
        paramsTry(:,3) = rand(size(paramsTry(:,1)))-0.5;
    case 'fix'
        paramsTry(:,3) = zeros(size(paramsTry(:,1)));
end

switch lower(lapserates(1:3))
    case 'con'
        paramsTry(:,4) = ones(size(paramsTry(:,1)))*rand(1)-0.5;
    case 'unc'
        paramsTry(:,4) = rand(size(paramsTry(:,1)))-0.5;
    case 'fix'
        paramsTry(:,4) = zeros(size(paramsTry(:,1)));
end

paramsTry = params + paramsTry.*rangeTries;