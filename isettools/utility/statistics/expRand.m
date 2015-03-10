function val = expRand(mn, S)
% Exponentially distributed random number generator
%
%   val = expRand(mn, S)
%
% mean is the distribution mean and S is the size of the returned matrix.
%
% It computes the formula
%   
%   val = -mn.*log(rand(S))
%   This generates a random variable with 
%    1/mn * exp(-x/mn) as probability density.
%
% Example:
%   S = [100,100];
%   mn = 3;
%   val = expRand(mn,S); hist(val(:),50);  mean(val(:))
%
% (c) 2010 Stanford Synapse Team 

%% Check input
if notDefined('mn'), mn = 1; end;
if notDefined('S'), S = 1; end;

val = -mn .* log(rand(S));

end