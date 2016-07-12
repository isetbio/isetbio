%
% PAL_SDT_AS_PCtoSL 
%
% uses an iterative search procedure of PAL_SDT_AS_SLtoPC to calculate the 
% stimulus level x from proportion correct PC, a PC that would result 
% from the ADDITIVE SUMMATION of n identical stimuli in Q monitored 
% channels, if x was subject to a scaling factor g and a transducer 
% exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task under the assumptions of 
% signal-detection-theory, assuming an unbiased observer.
%
% Syntax: [x]=PAL_SDT_AS_PCtoSL(PC,g,p,M,Q,n);
% 
% returns a scalar, vector or matrix of stimulus level 'x' for a scalar, 
% vector or matrix of proportion correct ('PC') in the ranges (1/M)<PC<1.0, 
% a scalar g (0<g<inf), scalar p (0<p<inf), scalar integer M (2<M<inf), 
% scalar integer Q (n=<Q<inf) and scalar integer n (0<n=<Q).
%
% Example:
%
% [x]=PAL_SDT_AS_PCtoSL([0.55 0.65 0.75 0.95],1,1,2,6,3)
%
% returns:
% 
% x =
%
%   0.1451    0.4449    0.7788    1.8993
%
% The example input consists of an N=4 vector of PCs with scalars g set 
% to 1, p set to 1, M to 2, Q to 6 and n to 3.  The output is the 
% corresponding N=4 vector of x
%
% Introduced: Palamedes version 1.7.0 (FK&NP)

function x=PAL_SDT_AS_PCtoSL(PC,g,p,M,Q,n)

%Checking the input arguments
if nargin ~= 6
    message = 'Wrong number of input arguments: input requires a scalar, vector or matrix of proportion correct PC, a scalar stimulus intensity scaling factor g, a scalar transducer exponent p, a scalar number of task alternatives M, a scalar number of monitored channels Q, and a scalar number of stimuli/signals n';
    error('PALAMEDES:nargin',message);
end

[rows, cols]=size(PC);
x = zeros(rows,cols);

func=@PAL_SDT_AS_SLtoPC;

for r=1:rows
    for c=1:cols
        x(r,c)=PAL_minimize(@PAL_sqDistanceYfuncX,1,[],PC(r,c),func,g,p,M,Q,n);
        if PC(r,c)==1.0
            x(r,c)=1/0;
        end
    end
end