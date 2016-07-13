%
% PAL_SDT_PS_uneqSLtoPC(x,g,p,M,Q)
%
% uses numerical integration to calculate proportion correct based 
% on the PROBABILITY SUMMATION of unequal stimulus levels in Q monitored 
% channels, with each stimulus level x subject to a scaling factor g 
% and a transducer with exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task, under the assumptions of SDT and 
% assuming an unbiased observer
%
% Syntax: [PC]=PAL_SDT_PS_uneqSLtoPC(x,g,p,M,Q);
% 
% returns a scalar proportion correct ('PC') for a vector of x each in 
% the range 0<x'<inf, a vector of g (0<g<inf), a vector of p (0<p<inf) 
% such that d'=(gx)'^p, a scalar integer M (2<M<inf), scalar integer 
% Q (n=<Q<inf) and scalar integer n (0<n=<Q).
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
% [PC]=PAL_SDT_PS_uneqSLtoPC([0.65 0.5 1.0],[1 1 1],[1 1 1],2,6)
%
% returns:
% 
% PC = 
%
%   0.6955
%
% The input consists of an N=3 vector of x with an N=3 vector of
% g, an an N=3 vector of p, an integer scalar M=2 and Q=6 monitored 
% channels. The output is proportion correct PC
%
% Introduced: Palamedes version 1.7.0 (FK&NP)                                

function [PC]=PAL_SDT_PS_uneqSLtoPC(x,g,p,M,Q)

%Checking the input arguments
if nargin ~= 5
    message = 'Wrong number of input arguments: requires a vector of stimulus level x, a vector of corresponding stimulus level scaling factors g, a vector of corresponding transducer exponents p, a scalar number of task alternatives M, and a scalar number of monitored channels Q';
    error('PALAMEDES:nargin',message);
end

n=length(x);
x(x<0.0)=0.0;
dP=(g.*x).^p;

PC=0.0;

if exist('quadgk.m','file')==2
    
    limit = Inf;
    if exist('OCTAVE_VERSION')
        limit = 100;    %Octave does not like Inf as limit
    end

    for i = 1:n
        PC=PC+quadgk(@(X)PAL_SDT_PS_uneqDPtoPCpartFunc(X,dP,i,M,Q,n),-limit,limit);
    end
        PC=PC+quadgk(@(X)PAL_SDT_PS_uneqDPtoPCpartFunc2(X,dP,M,Q,n),-limit,limit);
    
else
    if (max(dP)>8)
        warning('d-prime values >8 will not yield accurate results');
    end
    for i = 1:n
        PC=PC+quadl(@(X)PAL_SDT_PS_uneqDPtoPCpartFunc(X,dP,i,M,Q,n),-12,12); 
    end
        PC=PC+quadl(@(X)PAL_SDT_PS_uneqDPtoPCpartFunc2(X,dP,M,Q,n),-12,12);
    
end