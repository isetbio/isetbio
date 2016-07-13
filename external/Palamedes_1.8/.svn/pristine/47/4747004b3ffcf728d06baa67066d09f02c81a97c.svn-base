%
%PAL_whatIs  Identifies variable type
%
%   syntax: type = PAL_whatIs(arg)
%
%   returns:
%       0 in case 'arg' is empty (i.e., arg = [])
%       1 in case 'arg' is a numeric array (e.g., arg = [1 2]) 
%       2 in case 'arg' is a string (e.g., arg = 'hello')
%       3 in case 'arg' is a function handle (e.g., arg = @PAL_Gumbel)
%       4 in case 'arg' is a structure (e.g., arg.num = [1 2], arg.str ...
%           = 'hello')
%      -1 in case 'arg' is none of the above
%
%   Example: PAL_whatIs('hello') returns 2
%    
%Introduced: Palamedes version 1.1.0 (NP)

function type = PAL_whatIs(arg)

type = -1;
if isnumeric(arg)
    type = 1;
    if isempty(arg)
        type = 0;
    end
end
if ischar(arg)
    type = 2;
end
if isa(arg, 'function_handle')
    type = 3;
end
if isstruct(arg)
    type = 4;
end