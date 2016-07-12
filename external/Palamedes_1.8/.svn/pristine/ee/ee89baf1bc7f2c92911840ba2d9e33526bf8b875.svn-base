%
% PAL_SDT_PS_MonteCarlo_uneqSLtoPC(x,g,p,M,Q)
% 
% uses Monte Carlo simulation to calculate proportion correct based 
% on the PROBABILITY SUMMATION of unequal stimulus levels x in Q monitored 
% channels, with each stimulus level x subject to a scaling factor g 
% and a transducer exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task, under the assumptions of SDT and 
% assuming an unbiased observer
%
% Syntax: [PC]=PAL_SDT_PS_MonteCarlo_uneqSLtoPC(x,g,p,M,Q);
% 
% returns a scalar proportion correct ('PC') for a vector of x each in 
% the range 0<x<inf, a vector of g (0<g<inf), a vector of p (0<p<inf) 
% such that d'=(gx)'^p, a scalar integer M (2<M<inf), scalar integer 
% Q (n=<Q<inf) and scalar integer n (0<n=<Q).
%
% Note: Normally one should use PAL_SDT_PS_uneqSLtoPC which implements the
% same function but using numerical integration which is much faster.
% This routine is included for verification and educational purposes
%
% Syntax: [PC]=PAL_SDT_PS_MonteCarlo_uneqSPtoPC(sP,p,M,Q);
% 
% returns a scalar proportion correct ('PC') for a vector of s's each in 
% the ranges 0<s'<inf, a vector of p (0<p<inf), scalar integer M (2<M<inf), 
% scalar integer Q (n=<Q<inf) and scalar integer n (0<n=<Q).
%
% Example:
%
% [PC]=PAL_SDT_PS_MonteCarlo_uneqSLtoPC([0.65 0.5 1.0],[1 1 1],[1 1 1],2,6)
%
% returns something like:
% 
% PC = 
%
%   0.6957
%
% The input consists of an N=3 vector of x, with an N=3 vector of g, an 
% N=3 vector of p, and integer scalars M=2 and Q=6.The output is proportion 
% correct
%
% Introduced: Palamedes version 1.7.0 (FK&NP)

function PC = PAL_SDT_PS_MonteCarlo_uneqSLtoPC(x,g,p,M,Q)

%Checking the input arguments
if nargin ~= 5
    message = 'Wrong number of input arguments: input requires a vector of different stimulus levels x, a vector of corresponding stimulus level scaling factors g, a vector of corresponding transducer exponents p, a scalar number of task alternatives M, and a scalar number of monitored channels Q';
    error('PALAMEDES:nargin',message);
end

numReps=1000000;
n=length(x);

R = randn(Q,M,numReps);

for i=1:n
    dP=(g(i).*x(i)).^p(i);
    R(i,1,:) = R(i,1,:) + dP;
end

[trash, maxIndex]=max(max(R));
PC = length(maxIndex(maxIndex==1))./numReps;