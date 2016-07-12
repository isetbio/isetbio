%
% PAL_SDT_PS_SLtoPC(x,g,p,M,Q,n)
%
% uses numerical integration to calculate proportion correct PC based 
% on the PROBABILITY SUMMATION of n identical stimuli in Q monitored 
% channels, with stimulus level x subject to a scaling factor g and 
% transducer exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task, under the assumptions of 
% signal-detection-theory and assuming an unbiased observer
%
% Syntax: [PC]=PAL_SDT_PS_SLtoPC(x,g,p,M,Q,n);
% 
% returns a scalar, vector or matrix of proportion correct ('PC') for a 
% scalar, vector or matrix of x in the ranges 0<x<inf, a scalar 
% g (0<g<inf), a scalar p (0<p<inf), scalar integer M (2<M<inf), scalar 
% integer Q (n=<Q<inf) and scalar integer n (0<n=<Q).
%
% Note that the Matlab function quadgk, which performs a numerical 
% integration between the bounds -Inf and Inf, is used as the default in 
% the routine.  For versions of Matlab that do not have quadgk, quadl is 
% used instead with the bounds set to -12 and 12.  The results using quadl 
% are accurate except in extreme cases where d' is very high (>8)
%
% Example:
%
% [PC]=PAL_SDT_PS_SLtoPC([0.25 0.5 0.75 1.0],1,1,2,6,3)
%
% returns:
% 
% PC =
%
%   0.5602    0.6275    0.6975    0.7655
%
% The example input consists of an N=4 vector of x with scalars g set to 1,
% p set to 1, M to 2, Q to 6 and n to 3.  The output is the corresponding 
% N=4 vector of proportion correct
%
% Introduced: Palamedes version 1.7.0 (FK&NP)

function [PC]=PAL_SDT_PS_SLtoPC(x,g,p,M,Q,n)

%Checking the input arguments
if nargin ~= 6
    message = 'Wrong number of input arguments: requires a scalar, vector or matrix of stimulus level x, a scalar stimulus level scaling factor g,  a scalar transducer exponent p, a scalar number of task alternatives M, a scalar number of monitored channels Q, and a scalar number of stimuli/signals n';
    error('PALAMEDES:nargin',message);
end

[rows, cols]=size(x);
x(x<0.0)=0.0;
dP=(g.*x).^p; 

PC = zeros(rows,cols);

if exist('quadgk.m','file')==2
    
    limit = Inf;
    if exist('OCTAVE_VERSION')
        limit = 100;    %Octave does not like Inf as limit
    end

    for r = 1:rows
        for c = 1:cols
            PC(r,c)=quadgk(@(X)PAL_SDT_PS_DPtoPCpartFunc(X,dP(r,c),M,Q,n),-limit,limit);
        end
    end
else
    if (max(dP)>8)
        warning('d-prime values >8 will not yield accurate results');
    end
    for r = 1:rows
        for c = 1:cols
            PC(r,c)=quadl(@(X)PAL_SDT_PS_DPtoPCpartFunc(X,dP(r,c),M,Q,n),-12,12);
        end
    end
    
end