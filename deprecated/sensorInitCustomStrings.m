function sensorInitCustomStrings(handles)
% Deprecated
% Initialize popup menu strings that select  sensor compute function.
%
%    sensorInitCustomStrings(handles);
%
% The handles are from the sensorImageWindow
%
% The routine reads the default list in the .fig file, and  then adding
% the list from vcSESSION.CUSTOM as well as the Add Custom line. 
%
%
% Copyright ImagEval Consultants, LLC, 2005.

error('Should not be called.')

% Read the default custom compute list, and merge it with the additional
% compute functions saved in the session structure.
% defaultSensorCompute = get(handles.popCustomCompute,'String');
customSensorCompute =  ieSessionGet('sensorComputeList');
list = cellMerge(defaultSensorCompute,customSensorCompute);
set(handles.popCustomCompute,'String',list);

fontChange = ieSessionGet('increasefontsize');
if ~isempty(fontChange), 
    ieFontChangeSize(ieSessionGet('sensorFigure'),fontChange);
end

return;