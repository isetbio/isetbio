%
% PAL_SDT_AS_2uneqSLtoPC(x,r,g,p,M,Q)
%
% uses numerical integration to calculate proportion correct based 
% on the ADDITIVE SUMMATION of two unequal stimulus intensities, one of 
% which is x the other r*x, in Q monitored channels, with each x subject 
% to a scaling factor g and transducer exponent p such that d'= (gx)^p, 
% for a M-AFC (M-alternative-forced-choice) task, under the assumptions of 
% SDT and assuming an unbiased observer
%
% Syntax: [PC]=PAL_SDT_AS_2uneqSLtoPC(x,r,g,p,M,Q);
% 
% returns a scalar proportion correct 'PC' for a scalar 'x' (0<x<inf) a 
% scalar 'r' (0<r<inf), a N=2 vector of 'g' (0<g<inf), a N=2 vector of 'p' 
% (0<p<inf), a scalar integer 'M' (2<M<inf) and a scalar integer 
% 'Q' (n=<Q<inf)
%
% Because there is only one value of x as an input argument, the other
% given by r*x, the routine is invertible, i.e. one can use 
% PAL_SDT_AS_PCto2uneqSL(PC,r,g,p,M,Q) to convert PC to x and r*x
%
% Note that the Matlab function quadgk, which performs a numerical 
% integration between the bounds -Inf and Inf, is used as the default in 
% the routine.  For versions of Matlab that do not have quadgk, quadl is 
% used instead with the bounds set to -12 and 12.  The results using quadl 
% are accurate except in extreme cases where d' is very high (>8)
%
%
% Example:
%
% [PC]=PAL_SDT_AS_2uneqSLtoPC(2.5,1.5,[0.3 0.3],[1.5 1.5],2,2)
%
% returns:
% 
% PC = 
%
%   0.8216
%
% The input consists of an scalar x=2.5, a scalar r=1.5, a N=2 vector of
% g, N=2 vector of p, an integer scalar M=2 and integer scalar Q=2. The 
% output is a single scalar proportion correct PC
%
% Introduced: Palamedes version 1.8.0 (FK)

function [PC]=PAL_SDT_AS_2uneqSLtoPC(x,r,g,p,M,Q)

%Checking the input arguments
if nargin ~= 6
    message = 'Wrong number of input arguments: requires a scalar of first stimulus level x, a scalar ratio of the second to first stimulus levels r, a N=2 vector of corresponding stimulus level scaling factors g, a N=2 vector of corresponding transducer exponents p, a scalar number of task alternatives M, and a scalar number of monitored channels Q';
    error('PALAMEDES:nargin',message);
end

x(2)=x(1).*r;
x(x<0.0)=0.0;
dP=(g.*x).^p;

sumDP=sum(dP)./(Q.^0.5);

PC=PAL_SDT_MAFC_DPtoPC(sumDP,M);