function vcSetFigureHandles(figType,hObject,eventdata,handles);
% Set figure handle information at opening window
%
%  vcSetFigureHandles(figType,hObject,eventdata,handles);
%
% Purpose:
%    Maintain the handles and main figure object representation
%
%  vcSetFigureHandles('ISA',hObject,eventdata,handles);
%
%  Shortly, this routine will go away and be replaced by ieSessionSet
%  entirely.
%
% Copyright ImagEval Consultants, LLC, 2005.

switch lower(figType)
    case 'main'
        ieSessionSet('mainwindow',hObject,eventdata,handles);
                
    case 'scene'
        ieSessionSet('scenewindow',hObject,eventdata,handles);
        
    case {'oi','opticalimage'}
        ieSessionSet('oiwindow',hObject,eventdata,handles);
        
    case {'isa','sensor'}
        ieSessionSet('sensorwindow',hObject,eventdata,handles);
        
    case {'conemosaic'}
        ieSessionSet('conemosaicwindow',hObject,eventdata,handles);
        
    otherwise
        error('Unknown figure type');
end


return;
