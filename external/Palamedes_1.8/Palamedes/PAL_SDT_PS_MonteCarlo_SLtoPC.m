%
% PAL_SDT_PS_MonteCarlo_SLtoPC(x,g,p,M,Q,n)
% 
% uses numerical integration to calculate proportion correct PC based 
% on the PROBABILITY SUMMATION of n identical stimuli in Q monitored 
% channels, with signal level x subject to a scaling factor g and 
% transducer exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task under the assumptions of 
% signal-detection-theory and assuming an unbiased observer
%
% Syntax: [PC]=PAL_SDT_PS_MonteCarlo_SLtoPC(x,g,p,M,Q,n);
% 
% returns a scalar, vector or matrix of proportion correct ('PC') for a 
% scalar, vector or matrix of x in the ranges 0<x<inf, a scalar 
% g (0<g<inf), a scalar p (0<p<inf), scalar integer M (2<M<inf), scalar 
% integer Q (n=<Q<inf) and scalar integer n (0<n=<Q).
%
% Note: Normally one should use instead PAL_SDT_PS_SLtoPC(x,g,p,M,Q,n), 
% which performs the same function but using numerical integration 
% which is much faster than Monte Carlo (approx. x 50).  The Monte Carlo 
% routine is included for verification and pedagogical purposes
%
%
% Example:
%
% [PC]=PAL_SDT_PS_MonteCarlo_SLtoPC([0.25 0.5 0.75 1.0],1,1,2,6,3)
%
% returns something like:
% 
% PC =
%
%   0.5598    0.6276    0.6984    0.7653
%
% The example input consists of an N=4 vector of x with scalars g set to 1,
% p set to 1, M to 2, Q to 6 and n to 3.  The output is the corresponding 
% N=4 vector of proportion correct
%
% Introduced: Palamedes version 1.7.0 (FK&NP)

function PC = PAL_SDT_PS_MonteCarlo_SLtoPC(x,g,p,M,Q,n)

%Checking the input arguments
if nargin ~= 6
    message = 'Wrong number of input arguments: input requires a scalar, vector or matrix of stimulus level x, a scalar stimulus level scaling factor g, a scalar transducer exponent p, a scalar number of task alternatives M, a scalar number of monitored channels, Q and a scalar number of stimuli/signals n';
    error('PALAMEDES:nargin',message);
end

numReps=1000000;

[rows, cols]=size(x);  
PC = zeros(rows,cols);

for r=1:rows
    for c=1:cols
        dP=(g.*x(r,c)).^p;
        R = randn(Q,M,numReps);
        R(1:n,1,:) = R(1:n,1,:)+dP;
        [trash, maxIndex]=max(max(R));
        PC(r,c) = length(maxIndex(maxIndex==1))./numReps;
    end
end