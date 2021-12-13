function vcSetFigureHandles(figType, hObject, eventdata, handles)
% Set figure handle information at opening window
%
% Syntax:
%   vcSetFigureHandles(figType, hObject, eventdata, handles);
%
% Description:
%    Maintain the handles and main figure object representation
%
%    vcSetFigureHandles('ISA', hObject, eventdata, handles);
%
%    Shortly, this routine will go away and be replaced by ieSessionSet
%    entirely. I do not know if this has been accomplished yet.
%
% Inputs:
%    figType   - String. String describing the figure type. Options are
%                'main', 'scene', 'oi', 'opticalImage', 'isa', 'sensor',
%                and 'coneMosaic'.
%    hObject   - Object. The figure object. Type should match string in
%                figType for desired functionality.
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/09/18  jnm  Formatting

switch lower(figType)
    case 'main'
        ieSessionSet('mainwindow', hObject, eventdata, handles);

    case 'scene'
        ieSessionSet('scenewindow', hObject, eventdata, handles);

    case {'oi', 'opticalimage'}
        ieSessionSet('oiwindow', hObject, eventdata, handles);
        
    case {'conemosaic'}
        ieSessionSet('conemosaicwindow', hObject, eventdata, handles);

    otherwise
        error('Unknown figure type');
end

return;
