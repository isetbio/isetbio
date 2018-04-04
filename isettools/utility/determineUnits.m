function units = determineUnits(defaultUnit, varargin)
% Function determines if unit other than default specified
%
% Syntax:
%   units = defaultUnits(defaultUnit, [varargin])
%
% Description:
%    Function to simplify unit designation outside of within every case
%    statement in every single function.
%
%    There are examples contained in the code. To access, type 'edit
%    determineUnits.m' in the Command Window.
%
% Inputs:
%    defaultUnit - String. Default unit for case.
%    varargin    - (Optional). parent varargin.
%
% Outputs:
%    units       - String. Units to be used for case.
%
% Optional key/value pairs:
%    None.
%

% History:
%    03/08/18  jnm  Created.

% Examples:
%{
    units = determineUnits('meters');
    units = determineUnits('um','cyclesPerDeg');
%}

if isempty(varargin)
    units = defaultUnit;
else
    units = varargin{1};
end

end