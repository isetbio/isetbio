function displayClose(varargin)
% GUI callback invoked when display window closes.
%
% Syntax:
%   displayClose([varargin])
%
% Description:
%    This is GUI callback function which is invoked when the display window
%    is to be closed.
%
% Inputs:
%    None required.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    **NEEDS TO BE FILLED OUT**
%

% History:
%    05/xx/14  HJ   Created May, 2014
%    05/18/18  jnm  Formatting

% Close figure
global vcSESSION;

if checkfields(vcSESSION, 'GUI', 'vcDisplayWindow')
    vcSESSION.GUI = rmfield(vcSESSION.GUI,'vcDisplayWindow');
end

closereq;

end