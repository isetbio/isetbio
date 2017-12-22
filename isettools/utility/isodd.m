function cal = isodd(x)
% Always nice to know if there is something odd going on ;).
%
% Syntax:
%   bool = isodd(x)
%
% Description:
%    A function to determine if the passed argument is an odd value
%
% Inputs:
%    x - The input variable
%
% Outputs:
%    cal - The calculated boolean value (true for odd, false for even)
%
% Notes:
%    * [Note: XXX - Perhaps there should be an 'iseven' routine]
%    * [Note: JNM - Changed the output variable from bool to cal in order
%      to avoid clobbering the function bool(b)]
%    * [Note: JNM - This function only returns reasonable answers when
%      provided with integer inputs.]

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/21/17  jnm  Formatting & variable renaming

% Examples:
%{
    if isodd(3), disp('hello world'), end;
%}

cal = mod(x, 2);

end