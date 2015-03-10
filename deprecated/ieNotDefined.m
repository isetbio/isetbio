function notDefined = ieNotDefined( varString )
%Test whether a variable (usually a function argument) is defined
%
% notDefined = ieNotDefined( varString )
%
% Purpose:
%    This routine is used to determine if a variable is defined in the
%    calling function's workspace.  A variable is defined if (a) it exists
%    and (b) it is not empty. This routine is used throughout the ISET code
%    to test whether arguments have been passed to the routine or a default
%    should be assigned.
%
% notDefined: 1 (true) if the variable is not defined in the calling workspace 
%             0 (false) if the variable is defined in the calling workspace
%
%  Defined means the variable exists and is not empty in the function that
%  called this function.  
%
%  This routine replaced many calls of the form
%    if ~exist('varname','var') | isempty(xxx), ... end
%
%    with 
%
%    if ieNotDefined('varname'), ... end
%
% Copyright ImagEval Consultants, LLC, 2003.

% warning('Use notDefined instead');
if (~ischar(varString)), error('Varible name must be a string'); end

notDefined = 0;  % Assume the variable is defined

str = sprintf('''%s''',varString);   
cmd1 = ['~exist(' str ',''var'') == 1'];
cmd2 = ['isempty(',varString ') == 1'];
cmd = [cmd1, ' | ',cmd2];

% If either of these conditions holds, then not defined is true.
notDefined = evalin('caller',cmd1);     % Check that the variable exists in the caller space
if notDefined, return;                  % If it does not, return with a status of 0
else 
    notDefined = evalin('caller',cmd2); % Check if the variable is empty in the caller space
    if notDefined
        return;
    end
end

return;
