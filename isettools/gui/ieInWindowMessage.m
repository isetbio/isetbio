function ieInWindowMessage(str, handles, duration)
% Place a message in a text box within an ISET window
%
% Syntax:
%   ieInWindowMessage(str, handles, duration)
%
% Description:
%    The message placed in the window is kept visible for duration seconds.
%    If duration is not sent, the default is to leave the message in place.
%    To clear the box, set str = [];
%
%    The code below contains examples of function usage. To access, type
%    'edit ieInWindowMessage.m' into the Command Window.
%
% Inputs:
%    str      - (Optional). String. The string to display. Default []
%    handles  - (Optional). Handle. Handle to window. Default display in CW
%    duration - (Optional). Integer. Duration to display. Default empty.
%               Only removes string from display if handles is present,
%               else the display in the Command Window does not change.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    03/01/18  jnm  Formatting

% Examples:
%{
   handles = ieSessionGet('opticalimagehandle');
   ieInWindowMessage('Hello World', handles, []);
   ieInWindowMessage('', handles);
%}

if notDefined('duration'), duration = []; end
if notDefined('str'), str = []; end
if notDefined('handles'), disp(str); return; end

% Place the string in the message area.
set(handles.txtMessage, 'String', str);

% If duration is set, replace the string with blank after duration seconds.
if ~isempty(duration)
    pause(duration);
    set(handles.txtMessage, 'String', '');
end

end
