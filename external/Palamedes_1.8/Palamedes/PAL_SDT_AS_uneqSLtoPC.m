%
% PAL_SDT_AS_uneqSLtoPC(x,g,p,M,Q)
%
% uses numerical integration to calculate proportion correct PC based 
% on the ADDITIVE SUMMATION of unequal stimulus levels x in Q monitored 
% channels, with each stimulus level x subject to a scaling factor g 
% and a transducer exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task under the assumptions of SDT and 
% assuming an unbiased observer
%
% Syntax: [PC]=PAL_SDT_AS_uneqSLtoPC(x,g,p,M,Q);
% 
% returns a scalar proportion correct ('PC') for a vector of x each in 
% the range 0<x<inf, a vector of g (0<g<inf), a vector of p (0<p<inf) 
% such that d'=(gx)'^p, a scalar integer M (2<M<inf), scalar integer 
% Q (n=<Q<inf) and scalar integer n (0<n=<Q).
%
% Syntax: [PC]=PAL_SDT_AS_uneqSLtoPC(x,g,p,M,Q);
% 
% returns a scalar, vector or matrix of proportion correct ('PC') for a 
% vector of x (0<x<inf), a vector of corresponding g (0<g<inf), a vector 
% of corresponding p (0<p<inf) and a scalar integer M (2<M<inf).
%
% For two stimuli, with g=1 and p=1 (i.e. d'=x), M=2 and Q=2,the function 
% implements equation B10 in Shimozaki et al. (2003) 
% J. Opt. Soc. Amer., 20, 2197-2215.
%
% Example:
%
% [PC]=PAL_SDT_AS_uneqSLtoPC([0.25 0.5 1.0],[1 1 1],[1.5 1.5 1.5],2,6)
%
% returns:
% 
% PC =
%
%     0.6652
%
% The example input consists of an N=3 vector of x with an N=3 vector of 
% g set to 1, an N=3 vector of p set to 1.5, M to 2 and Q to 6.  The output 
% is a scalar of proportion correct
% 
%
% Introduced: Palamedes version 1.7.0 (FK&NP)


function [PC]=PAL_SDT_AS_uneqSLtoPC(x,g,p,M,Q)

%Checking the input arguments
if nargin ~= 5
    message = 'Wrong number of input arguments: input requires a vector of timulus level x, a vector of corresponding stimulus level scaling factors g, a vector of corresponding transducer exponents p, a scalar number of task alternatives M, and a scalar number of monitored channels Q';
    error('PALAMEDES:nargin',message);
end

x(x<0.0)=0.0;
dP=(g.*x).^p;
sumDP=sum(dP)./(Q.^0.5);

PC=PAL_SDT_MAFC_DPtoPC(sumDP,M);