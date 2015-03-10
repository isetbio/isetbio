function varargout = sceneWindow(varargin)
% Graphical user interface to manage the ISET SCENE properties.
%
%     varargout = sceneWindow(varargin)
%
%      SCENEWINDOW, by itself, creates a new SCENEWINDOW or raises the existing
%      singleton.
%
%      H = SCENEWINDOW returns the handle to a new SCENEWINDOW or the handle to
%      the existing singleton*.
%
%      SCENEWINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCENEWINDOW.M with the given input arguments.
%
%      SCENEWINDOW('Property','Value',...) creates a new SCENEWINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sceneWindow_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sceneWindow_OpeningFcn via varargin.
%
% Copyright ImagEval Consultants, LLC, 2003.

% Last Modified by GUIDE v2.5 28-Sep-2013 23:13:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @sceneWindow_OpeningFcn, ...
    'gui_OutputFcn',  @sceneWindow_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before sceneWindow is made visible.
function sceneWindow_OpeningFcn(hObject, eventdata, handles, varargin)

% sceneOpen.p performs essential functions, but only if the key is
% verified. Hence, commenting out this test will cause the sceneWindow to
% fail.
if ~sceneOpen(hObject,eventdata,handles), 
    error('Key or license failure');   
else sceneRefresh(hObject, eventdata, handles); 
end

return

% --- Outputs from this function are returned to the command line.
function varargout = sceneWindow_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;

return

% --- Executes during object creation, after setting all properties.
function editDistance_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
return

% --------------------------------------------------------------------
function editDistance_Callback(hObject, eventdata, handles)

% Should be set(SCENE,'editDistance',value);
[val,scene] = vcGetSelectedObject('SCENE');
if ~isempty(scene)
    scene.distance = str2double(get(hObject,'String'));
    scene.consistency = 0;
    vcReplaceObject(scene,val);
end
sceneRefresh(hObject, eventdata, handles);
return;

% --- Executes during object creation, after setting all properties.
function editLuminance_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

return;

% --------------------------------------------------------------------
function editLuminance_Callback(hObject, eventdata, handles)

[val,scene] = vcGetSelectedObject('SCENE');
if ~isempty(scene)
    meanL = str2double(get(hObject,'String'));
    scene = sceneAdjustLuminance(scene,meanL);
    vcReplaceObject(scene,val);
end
sceneRefresh(hObject, eventdata, handles);

return;


% --- Executes during object creation, after setting all properties.
function editHorFOV_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
return;

% --------------------------------------------------------------------
function editHorFOV_Callback(hObject, eventdata, handles)

[val,scene] = vcGetSelectedObject('SCENE');
if ~isempty(scene)
    scene = sceneSet(scene,'fov',str2double(get(hObject,'String')));
    vcReplaceObject(scene,val);
end
sceneRefresh(hObject, eventdata, handles);

return;

% --- Executes during object creation, after setting all properties.
function popupSelectScene_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
return;

% --- Executes on selection change in popupSelectScene.
function popupSelectScene_Callback(hObject, eventdata, handles)

contents = get(hObject,'String');
if strcmp(contents,'No Scene')
    return;
else
    val = get(hObject,'Value');
    vcSetSelectedObject('SCENE',val);
end
sceneRefresh(hObject, eventdata, handles);

return;

% --- Executes during object creation, after setting all properties.
function editRow_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
return;

% --------------------------------------------------------------------
function editRow_Callback(hObject, eventdata, handles)

[val,scene] = vcGetSelectedObject('SCENE');
scene.nRows = str2double(get(hObject,'String'));
vcReplaceObject(scene,val);

return;


% --- Executes during object creation, after setting all properties.
function editCol_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
return;


% --------------------------------------------------------------------
function editGamma_Callback(hObject, eventdata, handles)

% We store the display gamma value in the window without ever putting it
% into the scene structure.  We could ... but we don't.  Maybe we should.
sceneRefresh(hObject,eventdata,handles);


return;

% --- Executes on button press in btnInterpolate.
function btnInterpolate_Callback(hObject, eventdata, handles)
% Call back from the 'Interp' button
%
% Read the data in the row and col edit fields.  Re-sample the current data
% in the photons field so that it has the desired number of rows and
% columns.
[val,scene] = vcGetSelectedObject('SCENE');
sz = sceneGet(scene,'size');
r0 = sz(1); c0 = sz(2);

% Call GUI to set new row and col dimensions.
rc = sceneSetRowCol;
figure(gcbf);  % Return control to this figure.
newRow = rc(1); newCol = rc(2);
% These are the desired values determined by sceneSetRowCol

% Find the appropriate scale factor.
sFactor = [newRow/r0, newCol/c0];
scene = sceneInterpolate(scene,sFactor);
vcReplaceObject(scene,val);
sceneRefresh(hObject, eventdata, handles);

return;


%%%%%%%%%%%%%%%%%%%% Menus are controlled below here %%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
return

% --------------------------------------------------------------------
function EditMenu_Callback(hObject, eventdata, handles)
return

% --------------------------------------------------------------------
function PlotMenu_Callback(hObject, eventdata, handles)
return

% --------------------------------------------------------------------
function PlotLuminance_Callback(hObject, eventdata, handles)

[val,scene] = vcGetSelectedObject('SCENE');

% Check that the luminance field has been calculated
if ~checkfields(scene,'data','luminance')
    [lum, meanL] = sceneCalculateLuminance(scene);
    scene = sceneSet(scene,'luminance',lum);
    scene = sceneSet(scene,'meanLuminance',meanL);
    vcReplaceAndSelectObject(scene,val);
end

% Plots log10 or linear luminance
% Should be replaced by plotScene call
plotScene(scene,'luminance mesh log');

return;

% --------------------------------------------------------------------
function menPlotLumLin_Callback(hObject, eventdata, handles)

[val,scene] = vcGetSelectedObject('SCENE');

if ~checkfields(scene,'data','luminance')
    [scene.data.luminance, scene.data.meanL] = sceneCalculateLuminance(scene);
    vcReplaceAndSelectObject(scene,val);
end

% Plots log10 or linear luminance as a mesh.
% Should be replaced by plotScene call
plotScene(scene,'luminance mesh linear');

return;

% --------------------------------------------------------------------
function editGamma_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
return;

% --------------------------------------------------------------------
function popupDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% This is what Matlab does
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% I changed it to be like the others.  Not sure what's right.  Check on
% Ubuntu some day.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

return

% --------------------------------------------------------------------
function popupDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to popupDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupDisplay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupDisplay

% We read the state of this popup inside of sceneSetEditsAndButtons.  The
% value it has (contents{get(hObject,'Value')}) determines the displayFlag
% value.  1 means Standard RGB (the default) and we will be adding other
% display modes for NIR/SWIR and so forth over time.
sceneRefresh(hObject, eventdata, handles);
return


% --------------------------------------------------------------------
function plotRadiance_Callback(hObject, eventdata, handles)
% Plot | Radiance (Quanta)
% plotSceneRadiance('photons');
plotScene(vcGetObject('scene'),'radiance photons roi');
return

% --------------------------------------------------------------------
function menuPlotRadianceE_Callback(hObject, eventdata, handles)
% Plot | Radiance (Energy)
plotScene(vcGetObject('scene'),'radiance energy roi');
return

% --------------------------------------------------------------------
function menuPlotReflectance_Callback(hObject, eventdata, handles)
% Plot | Reflectance
plotScene(vcGetObject('scene'),'reflectance');
return

% --------------------------------------------------------------------
function menuPlotIlluminant_Callback(hObject, eventdata, handles)
% Plot | Illuminant (energy)
plotScene(vcGetObject('scene'),'illuminant energy roi');
return

% --------------------------------------------------------------------
function menuPlotIllumPhotons_Callback(hObject, eventdata, handles)
% Plot | Illuminant (photons)
plotScene(vcGetObject('scene'),'illuminant photons roi');
return

% --------------------------------------------------------------------
function menuPlotIlluminantImage_Callback(hObject, eventdata, handles)
%  Plot | Illuminant image
plotScene(vcGetObject('scene'),'illuminant image');
return

% --------------------------------------------------------------------
function menuPlotIlluminantComment_Callback(hObject, eventdata, handles)
% Deprecated - Plot | Illuminant comment
scene = vcGetObject('SCENE');
disp(sceneGet(scene,'illuminantComment'));
return

% --------------------------------------------------------------------
function menuPlotRadImGrid_Callback(hObject, eventdata, handles)
% Plot | Radiance image (grid)
plotScene(vcGetObject('SCENE'), 'radianceimagewithgrid');
return

% --------------------------------------------------------------------
function menuPlotWavebandImage_Callback(hObject, eventdata, handles)
% Plot | Waveband image
plotScene(vcGetObject('SCENE'), 'radiance waveband image');
return

% --------------------------------------------------------------------
function menuPlotImTrueSize_Callback(hObject, eventdata, handles)
% Plot | True Size Image
% Shows the image in the window in a separate window, with just the image
% The spatial scale is 1 to 1 with the data (true size)

% Get the data
scene = vcGetObject('scene');
spd = sceneGet(scene,'photons');
w   = sceneGet(scene,'wave');
gam = str2double(get(handles.editGamma,'String'));

% Create the image
RGB = imageSPD(spd,w,gam);

% Save the current graph window figure;
noNewGraphWin = 1;
figNumSave = vcSelectFigure('GRAPHWIN',noNewGraphWin);

% Create a new figure for the truesize image
figNumTRUESIZE = figure;
imshow(RGB);
truesize
% Turn off the menu and so forth
set(figNumTRUESIZE,'menubar','none','name','ISET');

% Make sure the saved image is reset to the current graph window
if ~isempty(figNumSave), ieSessionSet('graphWinFigure',figNumSave); end

return


% --------------------------------------------------------------------
function menuPlotMultipleRGB_Callback(hObject, eventdata, handles)
%Plot | Multiple image (RGB)
% Multiple figures with images of the session scene data
imageMultiview('scene');
return;

% --------------------------------------------------------------------
function menuPlotNewGraphWindow_Callback(hObject, eventdata, handles)
% Plot | New Graph Window
vcNewGraphWin;
return

% --------------------------------------------------------------------
function menuFileSave_Callback(hObject, eventdata, handles)
[val,scene] = vcGetSelectedObject('SCENE');
vcSaveObject(scene);
return

% --------------------------------------------------------------------
function menuFileLoad_Callback(hObject, eventdata, handles)
vcImportObject('SCENE');
sceneRefresh(hObject, eventdata, handles);
return

% --------------------------------------------------------------------
function menuSaveImage_Callback(hObject, eventdata, handles)
% File | Save (.tif)
gam = str2double(get(handles.editGamma,'String'));
[val, scene] = vcGetSelectedObject('SCENE');
sceneSaveImage(scene,[],gam);
return;

% --------------------------------------------------------------------
function menuFileClose_Callback(hObject, eventdata, handles)
sceneClose;
return;

% --------------------------------------------------------------------
function menuCopyScene_Callback(hObject, eventdata, handles)
[val,scene] = vcGetSelectedObject('scene');

newName = ieReadString('New scene name','new-scene');
if isempty(newName),  return;
else                  scene = sceneSet(scene,'name',newName);
end

vcAddAndSelectObject('scene',scene);
sceneRefresh(hObject, eventdata, handles);  
return;

% --------------------------------------------------------------------
function menuAn_Callback(hObject, eventdata, handles)
% Menu->Analyze callback
return;

% --------------------------------------------------------------------
function menuAnalyzeLine_Callback(hObject, eventdata, handles)
return;

% --------------------------------------------------------------------
function menuAnalyzeLineH_Callback(hObject, eventdata, handles)
plotScene(vcGetObject('SCENE'),'luminance hline');
return;

% --------------------------------------------------------------------
function menuAnalyzeLineV_Callback(hObject, eventdata, handles)
plotScene(vcGetObject('SCENE'),'luminance vline');
return;

% --------------------------------------------------------------------
function menuAnalyzeLFFTv_Callback(hObject, eventdata, handles)
plotScene(vcGetObject('SCENE'),'luminance fft vline');
return;

% --------------------------------------------------------------------
function menuAnalyzeLFFTH_Callback(hObject, eventdata, handles)
plotScene(vcGetObject('SCENE'),'luminance fft hline');
return;

% --------------------------------------------------------------------
function menuAnalyzeLineWave_Callback(hObject, eventdata, handles)
return;

% --------------------------------------------------------------------
function menuAnalyzeLWH_Callback(hObject, eventdata, handles)
scene = vcGetObject('SCENE'); 
plotScene(scene,'radiance hline');
return;

% --------------------------------------------------------------------
function menuAnalyzeLWV_Callback(hObject, eventdata, handles)
scene = vcGetObject('SCENE'); 
plotScene(scene,'radiance vline');
return;

% --------------------------------------------------------------------
function menuAnalyzeROI_Callback(hObject, eventdata, handles)
% Analyze | ROI
return;

% --------------------------------------------------------------------
function menuLuminance_Callback(hObject, eventdata, handles)
% Analyze | ROI Summary | Luminance
scene = vcGetObject('scene');
plotScene(scene,'luminance roi');
return;


% --------------------------------------------------------------------
function menuAnIlluminantCCT_Callback(hObject, eventdata, handles)
% Analyze | Illuminant CCT

scene = vcGetObject('scene');
wave = sceneGet(scene,'wave');
spd = sceneGet(scene,'illuminant energy');

% This size makes the title visible
str = sprintf('       ---------  Correlated color temp %.0f  ---------',spd2cct(wave,spd));
msgbox(str,'Illuminant');

return

% --------------------------------------------------------------------
function menuPlotCIE_Callback(hObject, eventdata, handles)
% Analyze | ROI Summary | Chromaticity
scene = vcGetObject('SCENE');
plotScene(scene,'chromaticity roi');
return;

% --------------------------------------------------------------------
function menuPlotDepth_Callback(hObject, eventdata, handles)
% Plot | Depth Map

scene = vcGetObject('scene');

if isempty(sceneGet(scene,'depth map'))
    handles = ieSessionGet('sceneimagehandle');
    ieInWindowMessage('No depth map data.',handles,3);
else
    plotScene(scene,'depth map');
end
return

% --------------------------------------------------------------------
function menuPlotDepthContour_Callback(hObject, eventdata, handles)
% Plot | Depth Contour

scene = vcGetObject('scene');
if isempty(oiGet(scene,'depth map'))
    handles = ieSessionGet('scene handle');
    ieInWindowMessage('No depth data.',handles,3);
else
    plotScene(scene,'depth map contour');
end

return

% --- Executes during object creation, after setting all properties.
function popupImScale_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
return;

% --- Executes on selection change in popupImScale.
function popupImScale_Callback(hObject, eventdata, handles)
% Call back for the Adjust scene popup that scales the scene size.
%

contents = get(hObject,'String');
str = contents{get(hObject,'Value')};
switch lower(str)
    case 'x 1'
        return;
    case 'x 2'
        sFactor = 2;
    case 'x 4'
        sFactor = 4;
    case 'x 1/2'
        sFactor = 1/2;
    case 'x 1/4'
        sFactor = 1/4;
    otherwise
        error('Unknown scale factor');
end

[val,scene] = vcGetSelectedObject('SCENE');
scene = sceneInterpolate(scene,sFactor);
vcReplaceObject(scene,val);
sceneRefresh(hObject, eventdata, handles);

return;


% --------------------------------------------------------------------
function menuFileRef_Callback(hObject, eventdata, handles)
sceneRefresh(hObject, eventdata, handles);
return;


% --------------------------------------------------------------------
function sceneRefresh(hObject, eventdata, handles)
% Refresh callback.
sceneSetEditsAndButtons(handles)
return;

% --------------------------------------------------------------------
function menuScene_Callback(hObject, eventdata, handles)
return;

% --------------------------------------------------------------------
function menuSceneMacbeth_Callback(hObject, eventdata, handles)
return;

% --------------------------------------------------------------------
function menuSceneMacbethC_Callback(hObject, eventdata, handles)
val = vcNewObjectValue('SCENE');
scene =  sceneCreate('macbethC');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuSceneMacbethTungsten_Callback(hObject, eventdata, handles)

val = vcNewObjectValue('SCENE');
scene =  sceneCreate('macbethTungsten');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuSceneMacbethD50_Callback(hObject, eventdata, handles)

val = vcNewObjectValue('SCENE');
scene =  sceneCreate('macbethD50');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuSceneMacbethFluorescent_Callback(hObject, eventdata, handles)

val = vcNewObjectValue('SCENE');
scene =  sceneCreate('macbethFluorescent');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuSceneMacbethD65_Callback(hObject, eventdata, handles)

val = vcNewObjectValue('SCENE');
scene =  sceneCreate('macbethD65');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuSceneMacbethVisIR_Callback(hObject, eventdata, handles)
% Scene | Macbeth Charts | Visible-InfraRed

spectrum.wave = ieReadNumber('Enter waves','380:4:1068','%.0f');
scene =  sceneCreate('macbethEE_IR',[],spectrum);
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuSceneLstar_Callback(hObject, eventdata, handles)
% Create vertical bars with equal L* steps

scene =  sceneCreate('lstar');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return

% --------------------------------------------------------------------
function menuSceneFile_Callback(hObject, eventdata, handles)
return;

% --------------------------------------------------------------------
function menuSceneMultiSpec_Callback(hObject, eventdata, handles)
% Scene | Choose | Multispectral 
%
scene = sceneFromFile([],'multispectral');
if isempty(scene), return; end
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuSceneChooseRGB_Callback(hObject, eventdata, handles)
% Scene | Choose | RGB
%
scene = sceneFromFile([],'rgb');
if isempty(scene), return; end
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;


% --------------------------------------------------------------------
function menuFileChooseFileMono_Callback(hObject, eventdata, handles)
% Scene | Choose | Monochrome
%
scene = sceneFromFile([],'monochrome');
if isempty(scene), return; end
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;


% --------------------------------------------------------------------
function menuScenesTest_Callback(hObject, eventdata, handles)
return;

% --------------------------------------------------------------------
function menuSceneMackay_Callback(hObject, eventdata, handles)

val = vcNewObjectValue('SCENE');
scene = sceneCreate('mackay');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuScenesSweep_Callback(hObject, eventdata, handles)
val = vcNewObjectValue('SCENE');
scene = sceneCreate('sweep');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuSceneSlantedBar_Callback(hObject, eventdata, handles)
scene = sceneCreate('slantedBar');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

function menuSceneZonePlate_Callback(hObject, eventdata, handles)
scene = sceneCreate('zonePlate');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuSceneFreqOrient_Callback(hObject, eventdata, handles)
scene = sceneCreate('freqorientpattern');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuSceneCheckerboard_Callback(hObject, eventdata, handles)
scene = sceneCreate('checkerBoard');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuScenePointArray_Callback(hObject, eventdata, handles)
% Scene | Patterns | PointArray (D65)

scene = sceneCreate('pointarray',[],[],'d65');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;


% --------------------------------------------------------------------
function menuSceneGridLines_Callback(hObject, eventdata, handles)
% Scene | Patterns | PointArray (D65)

scene = sceneCreate('GridLines',[],[],'d65');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

function menuSceneRadialLines_Callback(hObject, eventdata, handles)
scene = sceneCreate('radialLines');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuSceneNoise_Callback(hObject, eventdata, handles)

val = vcNewObjectValue('SCENE');
scene = sceneCreate('noise');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuHarmonic_Callback(hObject, eventdata, handles)

val = vcNewObjectValue('SCENE');
scene = sceneCreate('harmonic');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuSceneTestLine_Callback(hObject, eventdata, handles)
% Scene | Patterns | Line (D65 vert) 
% Impulse with a D65 spectral curve
%  

[val,scene] = vcGetSelectedObject('SCENE');
scene = sceneCreate('impulse1dd65',256);
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject,eventdata,handles);

return;

% --------------------------------------------------------------------
function menuScenesRamp_Callback(hObject, eventdata, handles)

val = vcNewObjectValue('SCENE');
scene = sceneCreate('ramp');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuUniform_Callback(hObject, eventdata, handles)
return;

% --------------------------------------------------------------------
function menuSceneUniformPhoton_Callback(hObject, eventdata, handles)
% Scene | Uniform | Equal Photon
val = vcNewObjectValue('SCENE');
scene = sceneCreate('uniformequalphoton');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;


% --------------------------------------------------------------------
function menuSceneUniformEE_Callback(hObject, eventdata, handles)
% Scene | Uniform | Equal Energy
val = vcNewObjectValue('SCENE');
scene = sceneCreate('uniformee');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuSceneUniformEESpecify_Callback(hObject, eventdata, handles)
% Scene | Uniform | Equal Energy (specify)
sz = 32;  % Spatial samples
wavelength = ieReadNumber('Enter waves','380:4:1068','%.0f');
if isempty(wavelength), disp('User canceled.'); return; end

scene = sceneCreate('uniformeeSpecify',sz,wavelength);
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuSceneUniformBBspecify_Callback(hObject, eventdata, handles)
% Scene | Uniform | Blackbody (specify)
% sz = 32;  % Spatial samples

% Read the size 
prompt = {'Image size','Color Temp','Wave (nm)'};
def = {'32','5000','400:700'};
answer = inputdlg(prompt,'Uniform blackbody',1,def);
if isempty(answer), disp('User canceled'); return;
else
    sz = str2double(answer{1});
    cTemp = str2double(answer{2});
    wave = eval(answer{3});
end

scene = sceneCreate('uniformbb',sz,cTemp,wave);
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuUniformD65_Callback(hObject, eventdata, handles)
% Scene | Uniform | D65
val = vcNewObjectValue('SCENE');
scene = sceneCreate('uniformD65');
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuEdit_Callback(hObject, eventdata, handles)
return;

% --------------------------------------------------------------------
function menuEditNewScene_Callback(hObject, eventdata, handles)
scene = sceneCreate;
vcAddAndSelectObject('SCENE',scene);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuEditClearWindow_Callback(hObject, eventdata, handles)
ieInWindowMessage('',ieSessionGet('sceneWindowHandles'),[]);
return

% --------------------------------------------------------------------
function menuEditFontSize_Callback(hObject, eventdata, handles)
ieFontChangeSize(handles.figure1);
return;

% --------------------------------------------------------------------
function editCrop_Callback(hObject, eventdata, handles)

[val,scene] = vcGetSelectedObject('SCENE');
[scene,rect] = sceneCrop(scene);

illF        = sceneGet(scene,'illuminant format');
switch illF
    case 'spatial spectral'
        % Get the illuminant photons out and put them in the main photon
        % slot.
        illuminantSPD = sceneGet(scene,'illuminant photons');
        sceneI = sceneSet(scene,'photons',illuminantSPD);
        
        % Crop them
        sceneI = sceneCrop(sceneI,rect);
        illuminantSPD = sceneGet(sceneI,'photons');
        
        % Stuff them back into the illuminant slot
        scene = sceneSet(scene,'illuminant photons',illuminantSPD);
    otherwise
end

vcReplaceObject(scene,val);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuEditSceneName_Callback(hObject, eventdata, handles)

[val,scene] = vcGetSelectedObject('SCENE');

newName = ieReadString('New scene name','new-scene');
if isempty(newName),  return;
else    scene = sceneSet(scene,'name',newName);
end

vcReplaceAndSelectObject(scene,val)
sceneRefresh(hObject,eventdata,handles);

return;

% --------------------------------------------------------------------
function menuEditDelete_Callback(hObject, eventdata, handles)
% Edit | Delete Current Scene
vcDeleteSelectedObject('SCENE');
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuEditDeleteSome_Callback(hObject, eventdata, handles)
% Edit | Delete Some Scenes
vcDeleteSomeObjects('scene');
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuEditTransform_Callback(hObject, eventdata, handles)
% hObject    handle to menuEditTransform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
return;

% --------------------------------------------------------------------
function menuEditTranspose_Callback(hObject, eventdata, handles)
% Edit | Transform | Transpose
[val,scene] = vcGetSelectedObject('SCENE');
photons     = imageTranspose(sceneGet(scene,'photons'));
scene       = sceneSet(scene,'cphotons',photons);
illF        = sceneGet(scene,'illuminant format');
switch illF
    case 'spatial spectral'
        photons     = imageTranspose(sceneGet(scene,'illuminant photons'),'leftRight');
        scene       = sceneSet(scene,'illuminant photons',photons);
    otherwise
end
vcReplaceObject(scene,val); sceneWindow();
return;

% --------------------------------------------------------------------
function menuEditRotate_Callback(hObject, eventdata, handles)
% Edit | Transform | Rotate
return;

% --------------------------------------------------------------------
function menuEditFlip_Callback(hObject, eventdata, handles)
% Edit | Transform | Flip
return;

% --------------------------------------------------------------------
function menuEditFlipHorizontal_Callback(hObject, eventdata, handles)
% Edit | Transform | Flip | Horizontal
[val,scene] = vcGetSelectedObject('SCENE');
photons     = imageFlip(sceneGet(scene,'photons'),'leftRight');
scene       = sceneSet(scene,'photons',photons);
illF        = sceneGet(scene,'illuminant format');
switch illF
    case 'spatial spectral'
        photons     = imageFlip(sceneGet(scene,'illuminant photons'),'leftRight');
        scene       = sceneSet(scene,'illuminant photons',photons);
    otherwise
end

vcReplaceObject(scene,val); sceneWindow();
return;

% --------------------------------------------------------------------
function menuEditFlipVertical_Callback(hObject, eventdata, handles)
% Edit | Transform | Flip | Vertical
[val,scene] = vcGetSelectedObject('SCENE');
photons     = imageFlip(sceneGet(scene,'photons'),'upDown');
scene       = sceneSet(scene,'photons',photons);

illF        = sceneGet(scene,'illuminant format');
switch illF
    case 'spatial spectral'
        photons     = imageFlip(sceneGet(scene,'illuminant photons'),'upDown');
        scene       = sceneSet(scene,'illuminant photons',photons);
    otherwise
end
vcReplaceObject(scene,val); sceneWindow();
return;

% --------------------------------------------------------------------
function menuEditRotCW_Callback(hObject, eventdata, handles)
% Edit | Transform | Rotate | ClockWise
[val,scene] = vcGetSelectedObject('SCENE');
photons     = imageRotate(sceneGet(scene,'photons'),'cw');
scene       = sceneSet(scene,'photons',photons);

illF        = sceneGet(scene,'illuminant format');
switch illF
    case 'spatial spectral'
        photons     = imageRotate(sceneGet(scene,'illuminant photons'),'cw');
        scene       = sceneSet(scene,'illuminant photons',photons);
    otherwise
end
vcReplaceObject(scene,val); sceneWindow();
return;

% --------------------------------------------------------------------
function menuEditRotCCW_Callback(hObject, eventdata, handles)
% Edit | Transform | Rotate | CounterClockWise
[val,scene] = vcGetSelectedObject('SCENE');
photons     = imageRotate(sceneGet(scene,'photons'),'ccw');
scene       = sceneSet(scene,'cphotons',photons);

illF        = sceneGet(scene,'illuminant format');
switch illF
    case 'spatial spectral'
        photons     = imageRotate(sceneGet(scene,'illuminant photons'),'ccw');
        scene       = sceneSet(scene,'illuminant photons',photons);
    otherwise
end
vcReplaceObject(scene,val); sceneWindow();
return;

% --------------------------------------------------------------------
function menuEditZoom_Callback(hObject, eventdata, handles)
zoom
return;

% --------------------------------------------------------------------
function menuEditViewer_Callback(hObject, eventdata, handles)
scene = vcGetObject('scene');
img = sceneGet(scene,'photons');
rgb = imageSPD(img,sceneGet(scene,'wavelength'));
ieViewer(rgb);

return;

% --------------------------------------------------------------------
function menuSpectral_Callback(hObject, eventdata, handles)
return;


% --------------------------------------------------------------------
function menuResampleWave_Callback(hObject, eventdata, handles)
[val,scene] = vcGetSelectedObject('SCENE');
scene = sceneInterpolateW(scene,[]);
vcReplaceObject(scene,val);
sceneRefresh(hObject,eventdata,handles);
return;


% --------------------------------------------------------------------
function menuEditAdjustMonochrome_Callback(hObject, eventdata, handles)
% Edit | Adjust SPD | Adjust Monochrome Wavelength

s = vcGetObject('scene');
w = sceneGet(s,'wavelength');
newWave = ieReadNumber('Enter new wavelength',w,'%.0f');
s = sceneSet(s,'wave',newWave);

vcReplaceAndSelectObject(s);
sceneRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuSPDdiv_Callback(hObject, eventdata, handles)

[val,scene] = vcGetSelectedObject('SCENE');
scene = sceneSPDScale(scene,[],'divide');
vcReplaceObject(scene,val);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuSPDMult_Callback(hObject, eventdata, handles)
% Edit | Adjust SPD | Multiply
[val,scene] = vcGetSelectedObject('SCENE');
scene = sceneSPDScale(scene,[],'multiply');
vcReplaceObject(scene,val);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuEditSetIlluminant_Callback(hObject, eventdata, handles)
% Edit | Adjust SPD | Change Illuminant

scene = vcGetObject('scene');
scene = sceneAdjustIlluminant(scene);

% Replace and refresh
vcReplaceObject(scene);
sceneRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuHelp_Callback(hObject, eventdata, handles)
return

% --------------------------------------------------------------------
% function menuHelpISETmanual_Callback(hObject, eventdata, handles)
% % Help | Iset manual (pdf)
% ieManualViewer('pdf','ISET_Manual');
% return

% --------------------------------------------------------------------
function menuHelpAppNotes_Callback(hObject, eventdata, handles)
% Help | Documentation (web)
web('http://imageval.com/documentation/','-browser');
return

% --------------------------------------------------------------------
function menuHelpSceneProgrammers_Callback(hObject, eventdata, handles)
% Help | Scene Programmers (online)
web('http://www.imageval.com/public/ISET-Functions/ISET/scene/index.html','-browser');
return

% --------------------------------------------------------------------
function menuHelpProgGuide_Callback(hObject, eventdata, handles)
% Help | Iset Programmers (online)
ieManualViewer('manual');
return


