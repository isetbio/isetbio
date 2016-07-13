%
%PAL_inverseWeibull
%
%PAL_inverseWeibull is no longer functional.
%
%Use x = PAL_Weibull(params,y,'inverse'); instead of 
%   x = PAL_inverseWeibull(params,y);
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.2.0, 1.4.0 (see History.m)

%Fixed error in help comments

function out = PAL_inverseWeibull(params, y)

message = 'PAL_inverseWeibull is no longer functional. Use ';
message = [message 'x = PAL_Weibull(params,y,''inverse'') instead '];
message = [message 'of x = PAL_inverseWeibull(params,y).'];
error(message);
out = [];