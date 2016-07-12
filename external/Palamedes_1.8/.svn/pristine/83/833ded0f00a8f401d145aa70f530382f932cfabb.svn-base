%
% PAL_SDT_MAFCoddity_IndMod_PCtoDP converts proportion correct into 
% d'(d-prime) for a M-AFC (M-alternative-forced-choice) oddity
% task, assuming an Independent Observer model and an unbiased observer
%
% Syntax: 
%   [dP] = PAL_SDT_MAFCoddity_IndMod_PCtoDP(pC,M,{optional arguments});
% 
% returns a scalar, vector or matrix of d' ('dP') for an input scalar, 
% vector or matrix of proportion correct p ('pC') and an input value of M 
% (number of alternatives), defined in the ranges 0<p<1 and 2<M<inf. Note 
% that the routine may take several seconds or even minutes to
% execute depending on the size of the input and speed of computer, as
% the routine performs an iterative search on the inverse version of the 
% routine which employs Monte Carlo simulation. Note also that because the 
% routine uses Monte Carlo simulation, no two outputs will be identical.  
% The speed/accuracy/memory tradeoff can be altered by changing the 
% parameter numReps using the optional argument 'numReps', followed by a 
% positive integer indicating the number of Monte Carlo simulations to be 
% used. Default value of numReps is 100000. Another manner in which to
% alter the speed/accuracy/memory tradeoff is to change the parameters of
% the Nelder-Mead Simplex search using the optional argument
% 'searchOptions', followed by an options structure created using options =
% PAL_minimize('options'). For more information, type
% PAL_minimize('options','help'). Example of usage shown below.
%
% Note that if M=3, you can use PAL_SDT_3AFCoddity_IndMod_PCtoDP, which
% is based on an inverse routine that is deterministic, making the
% routine more accurate and faster. 
%
% Example:
%
% [dP] = PAL_SDT_MAFCoddity_IndMod_PCtoDP([.5 .6 .7 .8 .9],6)
%
% returns something like:
% 
% dP =
% 
%     1.3437    1.6189    1.9000    2.2500    2.7436
%
% The example input arguments are an N=5 vector of proportion correct and
% a scalar M with a value of 6, and the output is an N=5 vector of d's
%
% The above example could be run with optional arguments to change the 
% speed/accuracy tradeoff.  For example, the following could be
% included in the body of the script prior to calling the routine
%
% options = PAL_minimize('options');
% options.TolX = 1e-4;      %default is 1e-6
% options.TolFun = 1e-4;    %ditto
%
% Now execute the routine, for example as follows:
%
% [dP] = PAL_SDT_MAFCoddity_IndMod_PCtoDP([.5 .6 .7 .8 .9],6,...
%   'numReps',50000,'searchOptions',options)
%
% will return something similar to the above but several times faster 
% (but with lower precision)
%
% Introduced: Palamedes version 1.6.0 (FK & NP)
% Modified: Palamedes version 1.6.3 (see History.m)

function dP = PAL_SDT_MAFCoddity_IndMod_PCtoDP(pC,M,varargin)

options = [];
va = {};

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'numReps',7)
            va{1} = varargin{n};
            va{2} = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'SearchOptions',7)
            options = varargin{n+1};
            valid = 1;
        end
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});
        end
    end            
end

[rows, cols] = size(pC);

func = @PAL_SDT_MAFCoddity_IndMod_DPtoPC;

dP = zeros(rows,cols);

for r = 1:rows
    for c = 1:cols
        dP(r,c) = PAL_minimize(@PAL_sqDistanceYfuncX,1,options,pC(r,c),func,M,va{:});
        if pC(r,c) == 1.0
            dP(r,c) = 1/0;
        end
    end
end