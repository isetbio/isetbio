function ndef = notDefined(varString)
% Test whether a variable (usually a function argument) is defined
%
% Syntax:
%	ndef = notDefined(varString)
%
% Description:
%    This routine is used to determine if a variable is defined in the
%    calling function's workspace. A variable is defined if (a) it exists
%    and (b) it is not empty. This routine is used throughout the ISET code
%    to test whether arguments have been passed to the routine or a default
%    should be assigned.
%
%    This routine replaced many calls of the form
%       'if ~exist('varname', 'var') || isempty(xxx), ... end'
%    with
%       'if notDefined('varname'), ... end'.
%
% Inputs:
%    varString - The string containing the name of the variable whose
%                existence you wish to verify.
%
% Outputs:
%    ndef      - A boolean indicating whether or not the variable is
%                defined. The return will be 1 (true) if the variable is
%                not defined in the parent workspace, and 0 (false) if the
%                variable was already defined. Defined meaning the variable
%                exists and is not empty in the calling function.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05  bw   summer 05, imported into mrVista 2.0
%    10/xx/05  ras  changed variable names to avoid a recursion error.
%    01/xx/06  ras  imported back into mrVista 1.0
%    01/xx/10  N    Nikhil, support for checking structure variables added
%    12/14/17  jnm  Formatting
%    01/19/18  jnm  Formatting update to match Wiki.

[rootVarString, fieldString] = strtok(varString, '.');
str = sprintf('''%s''', rootVarString);
cmd1 = ['~exist(' str ', ''var'') == 1'];
cmd2 = ['isempty(', rootVarString ') == 1'];

% create cmd3 if this is a structure
if ~isempty(fieldString)
    fieldString = fieldString(2:end);
    fieldString = strrep(fieldString, '.', ''', ''');
    cmd3 = ['~checkfields(', rootVarString, ', ''', fieldString ''')'];
end
cmd = [cmd1, ' || ', cmd2];

% If either of these conditions holds, then not defined is true
ndef = evalin('caller', cmd); % Check if variables exist in caller space
if ~ndef && ~isempty(fieldString)
    ndef = evalin('caller', cmd3); % Check if field exists in structure
end

end