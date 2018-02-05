function m = ieMemoryCheck(unit, level)
% Calculate the memory use in the calling routine or base
%
% Syntax:
%   m = ieMemoryCheck(unit, level)
%
% Description:
%    Calculate the memory use in the calling routine or base.
%
%    The output can be presented in the formats, megabytes ('MB'), 
%    kilobytes ('KB'), or bytes ('B')
%
% Inputs:
%    unit  - (Optional) The unit in which to present the output. Options
%            are 'mb', 'kb', and 'b'. Default is 'b' (bytes).
%    level - (Optional) The level of memory you are querying. Options are
%            'base' and 'caller'. Default is 'caller'.
%
% Outputs:
%    m     - The memory in use, in the specified units.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval, LLC, 2005
%    12/12/17  jnm  Formatting
%    01/19/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    m = ieMemoryCheck('KB', 'base')
    m = ieMemoryCheck('MB', 'caller')
%}

if notDefined('unit'), unit = 'b'; end
if notDefined('level'), level = 'caller'; end

m = 0; 
t = evalin(level, 'whos');  
for ss = 1:length(t), m = m + (t(ss).bytes); end

switch lower(unit)
    case 'mb'
        s = 1e6;
    case 'kb'
        s = 1e3;
    otherwise
        s = 1;
end

if nargout == 0, fprintf('Memory %f %s\n', m/s, unit); end

end
