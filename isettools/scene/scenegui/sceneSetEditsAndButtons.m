function sceneSetEditsAndButtons(handles)
% Fill scene window fields based on the current scene information
%
%    sceneSetEditsAndButtons(handles)
%
% Fill the fields  the current scene information, including editdistance,
% editluminance, FOV, etc.  
%
% Display of the image data is handled separately by sceneShowImage.
%
% Copyright ImagEval Consultants, LLC, 2003.


%% Use scene data to set boxes in window
scene = vcGetObject('SCENE');

if isempty(scene)
    % No scene, so set empty
    str = [];
    set(handles.editDistance,'String',str);
    set(handles.editLuminance,'String',str);
    set(handles.editHorFOV,'String',str);
    
    % Select scene popup contents
    set(handles.popupSelectScene,...
        'String','No Scenes',...
        'Value',1);
else
    % Text boxes on right: we should reduce the fields in SCENE.
    set(handles.editDistance,'String',num2str(sceneGet(scene,'distance')));
    meanL = sceneGet(scene,'mean luminance');
   
    set(handles.editLuminance,'String',sprintf('%.1f',meanL));
    set(handles.editHorFOV,'String',sprintf('%.2f',sceneGet(scene,'fov')));
    
    % Select scene popup contents
    set(handles.popupSelectScene,...
        'String',vcGetObjectNames('SCENE'),...
        'Value',vcGetSelectedObject('SCENE'));    
end

%% Description box on upper right
set(handles.txtSceneDescription,'String',sceneDescription(scene));

%% Set the gamma and displayFlag from the scene window.
figNum = vcSelectFigure('SCENE'); 
figure(figNum);
displayFlag = get(handles.popupDisplay,'Value');

gam = str2double(get(handles.editGamma,'String'));
%sceneShowImage(scene,displayFlag,gam);
sceneGet(scene, 'rgb', gam, displayFlag);

return;
