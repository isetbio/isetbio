function oiClose
% Close optical image window
%
% Syntax:
%   oiClose
%
% Description:
%    Close window function for optical image and remove figure handle from
%    vcSESSION structure.
%
% Inputs:
%    None.
%
%  Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/19/18  jnm  Formatting

global vcSESSION;

if checkfields(vcSESSION, 'GUI', 'vcOptImgWindow')
    vcSESSION.GUI = rmfield(vcSESSION.GUI, 'vcOptImgWindow');
end

closereq;

return;