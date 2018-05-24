function val = ieReadNumber(str, defaultValue, fmt)
% Graphical query for a number
%
% Syntax:
%   val = ieReadNumber([str], [defaultValue = 1], [fmt])
%
% Description:
%    Graphical wrapper to query the user for a number.  You can also enter
%    a vector, say by typing in 1:10, and the return values will be (1:10);
%
%    The code below contains examples of function usage. To access, type
%    'edit ieReadNumber.m' into the Command Window.
%
% Inputs:
%    str          - (Optional) String. User prompt. Default 'Enter value'.
%    defaultValue - (Optional) Numerical. Default value if nothing entered.
%                   Default 1.
%    fmt          - (Optional) String. Numerical format. Default '  %.2e'
%
% Outputs:
%    val          - Numerical. The user entered value in specified format.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/01/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - Requires user input.
    val = ieReadNumber('Enter column number')
    val = ieReadNumber('Enter column number', 1 / 17, ' %.3f')
%}

if notDefined('str'), str = 'Enter value'; end
if notDefined('defaultValue'), defaultValue = 1; end
if notDefined('fmt'), fmt = '  %.2e'; end

prompt = {str};
def = {num2str(defaultValue, fmt)};
dlgTitle = sprintf('ISET read number');
lineNo = 1;
answer = inputdlg(prompt, dlgTitle, lineNo, def);

if isempty(answer), val = []; else, val = eval(answer{1}); end

return