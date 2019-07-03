function handles = oiImageInitCustomStrings(handles)
% Initialize the strings in the OI window custom popup menu
%
% Syntax:
%   handles = oiInitCustomStrings(handles);
%
% Description:
%    The popup strings in the OI window are initialized by data in the the
%    .fig file. Then, the strings are modified using data contained in the
%    global variable, vcSESSION.CUSTOM.
%
%    The popup menu strings can also maniuplated using the Add Custom line.
%
% Inputs:
%    handles - Handle. The handle to the optical image's window.
%
% Outputs:
%    handles - Handle. The handle to the optical image's window.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    03/07/18  jnm  Formatting

defaultOICompute = get(handles.popCustom,'String');
customOI = ieSessionGet('oicomputelist');
list = [defaultOICompute, customOI];
set(handles.popCustom,'String',list);

end