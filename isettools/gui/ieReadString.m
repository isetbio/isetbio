function ieString = ieReadString(prompt, defString)
% Query the user for a string.
%
% Syntax:
%   ieString = ieReadString(prompt, defString)
%
% Description:
%    Query the user for an input string. If the user cancels, the returned
%    string is empty.
%
%    The code below contains examples of function usage. To access, type
%    'edit ieReadString.m' into the Command Window.
%
% Inputs:
%    prompt    - (Optional) String. The prompt to display to the user.
%                Default 'Enter'.
%    defString - (Optional). String. Default String to return. Default ''.
%
% Outputs:
%    ieString  - String. The return string.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    03/01/18  jnm  Formatting

% Example:
%{
    newName = ieReadString('Enter filter name: ', 'filter')
%}

if notDefined('prompt'), prompt = 'Enter'; end
if notDefined('defString'), defString = ''; end

ieString = [];

def = {defString};
dlgTitle = sprintf('IE Read String');
lineNo = 1;
answer = inputdlg(prompt, dlgTitle, lineNo, def);
if ~isempty(answer), ieString = answer{1}; end

end