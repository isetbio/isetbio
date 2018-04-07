function sformatted = ieParamFormat(s)
% Converts s to a standard ISET parameter format  (lower case, no spaces)
%
% Syntax:
%   sformatted = ieParamFormat(s)
%
% Description:
%    The string is sent to lower case and spaces are removed.
%  
%    If the input argument is a cell array, then ieParamFormat is
%    called on the odd entries of the array.  This allows conversion
%    of a varargin that contains only key/value pairs to a form where
%    only the keys are translated to standard ISET parameter format.
%
% Inputs:
%    s          - The entity that you wish to parameter-ize. If a singular
%                 string, convert the entirety to the parameter format (all
%                 lower case, no spaces). If an array, perform the
%                 conversion on odd entries only.
%
% Outputs:
%    sFormatted - The modified s (string or array).
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/10       Copyright ImagEval Consultants, LLC, 2010
%	 12/05/17  dhb  Handle cell arrays.
%    12/12/17  jnm  Formatting
%    01/19/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    ieParamFormat('Exposure Time')
%}
%{
    keyValuePairs{1} = 'Exposure Time';
    keyValuePairs{2} = 1;
    keyValuePairs{3} = 'iWasCamelCase';
    keyValuePairs{4} = 'Do Not Convert Me';
    keyValuePairs = ieParamFormat(keyValuePairs)
%}

if (~ischar(s) && ~iscell(s))
    error('s has to be a string or cell array');
end

% Lower case
if (ischar(s))
    % To lower and remove spaces
    sformatted = lower(s);
    sformatted = strrep(sformatted, ' ', '');
else
    if (iscell(s))
        sformatted = s;
        for ii = 1:2:length(s)
            sformatted{ii} = ieParamFormat(s{ii});
        end
    end
end
