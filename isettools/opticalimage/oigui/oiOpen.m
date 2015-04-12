function bool = oiOpen(hObject,eventdata,handles)
% Initialize oiWindow
% 
%    bool = oiOpen(hObject,eventdata,handles)
%
% Copyright ImagEval Consultants, LLC, 2005.

%
bool = true;

% Choose default command line output for microLensWindow
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

vcSetFigureHandles('OI',hObject,eventdata,handles);

%  Check the preferences for ISET and adjust the font size.
ieFontInit(hObject);

return;