function bool = isodd(x)
% Always nice to know if there is something odd going on ;).
%
% bool = isodd(x)
%
%   Perhaps there should be an 'iseven' routine
%
% Example:
%   if isodd(3), disp('hello world'), end;
%
% Copyright ImagEval Consultants, LLC, 2005.

bool = mod(x,2);

end