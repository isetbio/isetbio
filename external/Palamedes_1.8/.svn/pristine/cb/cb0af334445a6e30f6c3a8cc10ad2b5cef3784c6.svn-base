%
% PAL_SDT_PS_2uneqSLtoPC(x,r,g,p,M,Q)
%
% uses numerical integration to calculate proportion correct based 
% on the PROBABILITY SUMMATION of two unequal intensity stimuli, one 
% with an intensity of x the other r*x, in Q monitored channels, with each 
% x subject to a scaling factor g and transducer exponent p 
% such that d'= (gx)^p, in a M-AFC (M-alternative-forced-choice) task, 
% under the assumptions of SDT and assuming an unbiased observer
%
% Syntax: [PC]=PAL_SDT_PS_2uneqSLtoPC(x,r,g,p,M,Q);
% 
% returns a scalar proportion correct 'PC' for a scalar 'x' (0<x'<inf) a 
% scalar 'r' (0<r<inf), a N=2 vector of 'g' (0<g<inf), a N=2 vector of 'p' 
% (0<p<inf), a scalar integer 'M' (2<M<inf) and a scalar integer 
% 'Q' (n=<Q<inf).
%
% Because there is only one value of x as an input argument, the other
% given by r*x, the routine is invertible, i.e. one can use 
% PAL_SDT_PS_PCto2uneqSL(PC,r,g,p,M,Q) to convert PC to x and r*x
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
% [PC]=PAL_SDT_PS_2uneqSLtoPC(2.5,1.5,[0.3 0.3],[1.5 1.5],2,2)
%
% returns:
% 
% PC = 
%
%   0.7936
%
% The input consists of a scalar x=2.5, a scalar r=1.5, a N=2 vector of
% g, N=2 vector of p, an integer scalar M=2 and integer scalar Q=2. The 
% output is a single scalar proportion correct PC
%
% Introduced: Palamedes version 1.8.0 (FK)                                

function [PC]=PAL_SDT_PS_2uneqSLtoPC(x,r,g,p,M,Q)

%Checking the input arguments
if nargin ~= 6
    message = 'Wrong number of input arguments: requires a scalar of first stimulus level x, a scalar ratio of the second to first stimulus levels r, a N=2 vector of corresponding stimulus level scaling factors g, a N=2 vector of corresponding transducer exponents p, a scalar number of task alternatives M, and a scalar number of monitored channels Q';
    error('PALAMEDES:nargin',message);
end

x(2)=r.*x(1);
x(x<0.0)=0.0;
dP=(g.*x).^p;

PC=0.0;
n=2;

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