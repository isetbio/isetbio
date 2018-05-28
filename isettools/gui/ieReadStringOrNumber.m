function ieInput = ieReadStringOrNumber(prompt, defString)
% Query the user for a string or a number.
%
% Syntax:
%   newFilterName = ieReadString(prompt, defString)
%
% Description:
%    If the user cancels, the returned string is empty.
%
%    The function contains examples of usage below. To access these
%    examples, type 'edit ieReadStringOrNumber.m' into the Command Window.
%
% Inputs:
%    prompt    - String. Prompt string requesting user input.
%    defString - String. The Default string associated with the prompt.
%
% Outputs:
%    ieInput   - VARIES. The entered value, be it Numerical or String.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/02/18  jnm  Formatting

% Example:
%{
    newName = ieReadString('Enter filter name: ', 'filter')
%}

if notDefined('prompt'), prompt = 'Enter'; end
if notDefined('defString'), defString = ''; end

ieInput = [];

def = {defString};
dlgTitle = sprintf('IE Read String');
lineNo = 1;
answer = inputdlg(prompt, dlgTitle, lineNo, def);
if ~isempty(answer)
    ieInput = answer{1};
    % if number - convert
    if ~isempty(str2num(ieInput))
        ieInput = str2num(ieInput);
    end
end

return;