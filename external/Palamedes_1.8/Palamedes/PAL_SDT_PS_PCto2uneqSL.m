%
% PAL_SDT_PS_PCto2uneqSL 
%
% uses an iterative search procedure of PAL_SDT_PS_2uneqSLtoPC to 
% calculate two stimulus levels x with a ratio r from proportion correct PC, 
% a PC that would result from the PROBABILITY SUMMATION of the two 
% different stimuli, in Q monitored channels, with x subject to a scaling 
% factor g and transducer exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task under the assumptions of 
% signal-detection-theory, assuming an unbiased observer.
%
% Syntax: [x]=PAL_SDT_PS_PCto2uneqSL(PC,r,g,p,M,Q);
% 
% returns a N=2 vector of stimulus levels 'x' for a scalar value of 
% proportion correct 'PC' (1/M)<pC<1.0), a scalar 'r' (0<r<inf), a N=2 
% vector of 'g' (0<g<inf), a N=2 vector of 'p' (0<p<inf), a scalar integer 
% 'M' (2<M<inf) and a scalar integer 'Q' (2=<Q<inf).
%
% Example:
%
% [x]=PAL_SDT_PS_PCto2uneqSL(0.75,1.5,[1 1],[2 2],2,2)
%
% returns:
% 
% x =
%
%   0.6784    1.0176
%
% The example input consists of a scalar PC=0.75, scalar r=1.5, N=2 vector 
% of g, N=2 vector of p, scalar integer M=2 and scalar integer Q=2.  The 
% output is a N=2 vector of x
%
% Introduced: Palamedes version 1.8.0 (FK)

function [x]=PAL_SDT_PS_PCto2uneqSL(PC,r,g,p,M,Q)

options = PAL_minimize('options');   %PAL_minimize search options        


%Checking the input arguments
if nargin ~= 6
    message = 'Wrong number of input arguments: input requires a scalar prop. correct PC, a scalar ratio of stimulus levels r, a N=2 vector of scaling factors g, a N=2 vector of transducer exponents p, a scalar number of task alternatives M and a scalar number of monitored channels Q';
    error('PALAMEDES:nargin',message);
end

func=@PAL_SDT_PS_2uneqSLtoPC;

x(1)=PAL_minimize(@PAL_sqDistanceYfuncX,0,options,PC,func,r,g,p,M,Q); % note that guess x param is set to 0

if PC==1.0
    x=1/0;
end

x(2)=x(1).*r;

end