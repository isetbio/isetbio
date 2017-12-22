function txt = addText(txt, str)
% Add text to an existing text string
%
% Syntax:
%   txt = addText(txt, str)
%
% Description:
%    Utility for combining strings
%
% Inputs:
%    txt - Existing string
%    str - String to add to the existing string
%
% Outputs:
%    txt - Modified string to return
%
% Notes:
%    * [Note: JNM - mrMessage does not seem to be included by default?]
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/12/17  jnm  Formatting

% Example:
%{
    txt = 'Hello World! ';
    txt = addText(txt, 'What a beautiful day!');
    fprintf("%s\n", txt);
%}

txt = [txt, sprintf(str)];

end
