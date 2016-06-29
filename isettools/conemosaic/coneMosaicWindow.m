function varargout = coneMosaicWindow(varargin)
%Sensor image coneMosaicWindow interface
%
%   varargout = coneMosaicWindow(varargin)
%   CONEMOSAICWINDOW M-file for coneMosaicWindow.fig
%
%  Graphical user interface to manage the Image Sensor Array (ISA) properties.
%
%  CONEMOSAICWINDOW, by itself, creates a new CONEMOSAICWINDOW or raises the existing
%  singleton*.
%
%  H = CONEMOSAICWINDOW returns the handle to a new CONEMOSAICWINDOW or the handle to
%  the existing singleton*.
%
%  CONEMOSAICWINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%  function named CALLBACK in CONEMOSAICWINDOW.M with the given input arguments.
%
%  CONEMOSAICWINDOW('Property','Value',...) creates a new CONEMOSAICWINDOW or raises the
%  existing singleton*.  Starting from the left, property value pairs are
%  applied to the GUI before sensorImageWindow_OpeningFunction gets called.  An
%  unrecognized property name or invalid value makes property application
%  stop.  All inputs are passed to coneMosaicWindow_OpeningFcn via varargin.
%
%  *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% Copyright ImagEval Consultants, LLC, 2005.

% Last Modified by GUIDE v2.5 28-Jun-2016 20:20:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @coneMosaicWindow_OpeningFcn, ...
    'gui_OutputFcn',  @coneMosaicWindow_OutputFcn, ...
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
end

% --- Executes just before coneMosaicWindow is made visible.
function coneMosaicWindow_OpeningFcn(hObject, eventdata, handles, varargin)
%#ok<*DEFNU>
%#ok<*INUSD>

% check inputs
if isempty(varargin) || ~isa(varargin{1}, 'coneMosaic')
    error('cone mosaic object required');
end

% Choose default command line output for coneMosaicWindow
handles.output = hObject;
handles.cMosaic = varargin{1};

% Update handles structure
guidata(hObject, handles);
handles.cMosaic.window = hObject;

figure(hObject);
ieFontInit(hObject);

coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --- Outputs from this function are returned to the command line.
function varargout = coneMosaicWindow_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
end

function btnComputeImage_Callback(hObject, eventdata, handles)
% Button press computes the image from the optics data
[~, oi] = vcGetSelectedObject('OI');
if isempty(oi) || isempty(oiGet(oi, 'photons'))
    error('No optical image photon data available');
end
handles.cMosaic.compute(oi);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

function menuAnComputeFromOI_Callback(hObject, eventdata, handles)
btnComputeImage_Callback(hObject, eventdata, handles);
end

% Edit box - adjust number of rows
function editRows_Callback(hObject, eventdata, handles)
handles.cMosaic.rows = str2double(get(hObject, 'String'));
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% Edit box - adjust number of columns
function editCols_Callback(hObject, eventdata, handles)
handles.cMosaic.cols = str2double(get(hObject, 'String'));
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% Edit box - adjust integration time
function editExpTime_Callback(hObject, eventdata, handles)
handles.cMosaic.integrationTime = 1e-3*str2double(get(hObject, 'String'));
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% GUI object create functions
function editRows_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
end

function editCols_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
end

function editExpTime_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
end

% menu call back functions
function menuFile_Callback(hObject, eventdata, handles)
end

function menuFileClose_Callback(~, ~, handles)
close(handles.coneMosaicWindow);
end

function menuEdit_Callback(hObject, eventdata, handles)
end

function menuEditName_Callback(hObject, eventdata, handles)
str = ieReadString('New name', handles.cMosaic.name);
if ~isempty(str), handles.cMosaic.name = str; end
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

function sensorEditClearData_Callback(hObject, eventdata, handles)
handles.cMosaic.clearData();
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function editGam_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
end

function editGam_Callback(hObject, eventdata, handles)
coneMosaicGUIRefresh(hObject,eventdata,handles);
end

function coneMosaicGUIRefresh(~, ~, handles)
%Update the cone mosaic gui window interface
%
%   coneMosaciGUIRefresh(handles)
%
% HJ/BW, ISETBIO TEAM, 2016

% get coneMosaic object
cm = handles.cMosaic;

% set row and cols
set(handles.editRows, 'string', num2str(cm.rows));
set(handles.editCols, 'string', num2str(cm.cols));

% set integration time
set(handles.editExpTime, 'string', sprintf('%.1f',cm.integrationTime*1e3));

% set KLMS ratio
str = sprintf('[%.1f, %.1f, %.1f, %.1f]', cm.spatialDensity(1), ...
    cm.spatialDensity(2), cm.spatialDensity(3), cm.spatialDensity(4));
set(handles.editKLMS, 'string', str);

% set description strings
str = cm.description('skipMacular', true, 'skipPigment', true);
set(handles.txtMosaic, 'string', str);
set(handles.txtConeProperties, 'string', cm.pigment.description);

% set photopigment properties
set(handles.editConeWidth, 'string', num2str(cm.pigment.width*1e6));
set(handles.editConeHeight, 'string', num2str(cm.pigment.height*1e6));

str = sprintf('[%.1f, %.1f, %.1f]', cm.pigment.opticalDensity(1), ...
    cm.pigment.opticalDensity(2), cm.pigment.opticalDensity(3));
set(handles.editConeOpticalDensity, 'string', str);

str = sprintf('[%.2f, %.2f, %.2f]', cm.pigment.peakEfficiency(1), ...
    cm.pigment.peakEfficiency(2), cm.pigment.peakEfficiency(3));
set(handles.editConePeakEfficiency, 'string', str);

% set macular density
set(handles.editMacularDensity, 'string', num2str(cm.macular.density));

% set image content to axes
cm.plot('cone mosaic', 'hf', handles.axes2);

% check if absorptions and current are available
if isempty(cm.absorptions)
    set(handles.menuPlotMosaicMeanAbsorptions, 'Enable', 'off');
else
    set(handles.menuPlotMosaicMeanAbsorptions, 'Enable', 'on');
end

if isempty(cm.current)
    set(handles.menuPlotMosaicMeanCurrent, 'Enable', 'off');
else
    set(handles.menuPlotMosaicMeanCurrent, 'Enable', 'on');
end

end

% --------------------------------------------------------------------
function menuAnSNR_Callback(hObject, eventdata, handles)
% Analyze SNR menu
end

% --------------------------------------------------------------------
function menuAnHistogram_Callback(hObject, eventdata, handles)
% Analyze Histogram Menu
end

% --------------------------------------------------------------------
function menuISOSat_Callback(hObject, eventdata, handles)
%
speed = isoSpeed('saturation');
str = sprintf('ISO speed (saturation):    %.0f\n\n',speed);
str = [str,sprintf('Measured for a uniform D65 optical image.\n\n')];
str = [str,sprintf('Larger means saturates at lower lux-sec level.\n\n')];

msgbox(str);

end

% --------------------------------------------------------------------
function menuAnExposureValue_Callback(hObject, eventdata, handles)

EV = exposureValue(vcGetObject('OPTICS'),vcGetObject('ISA'));

str = sprintf('Exposure value (log2(f/#^2 / T)):    %.2f',EV);
ieInWindowMessage(str,ieSessionGet('sensorwindowhandles'));

end

%-------------Photometric Exposure Value (lux-sec)
function menuAnPhotExp_Callback(hObject, eventdata, handles)
% Analyze | SNR | Photometric Exp

str = sprintf('Photometric exposure (lux-sec): %.2f',...
    photometricExposure(vcGetObject('OI'),vcGetObject('ISA')));

% Display in coneMosaicWindow message
ieInWindowMessage(str,ieSessionGet('sensorwindowhandles'));

end

function menuPlot_Callback(hObject, eventdata, handles)
% Menu Plot
end

function menuAnalyze_Callback(hObject, eventdata, handles)
% Analyze menu
end

function menuAnLine_Callback(hObject, eventdata, handles)
% Analyze->Line menu
end

function menuHorizontal_Callback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
function menuVertical_Callback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
function menuAnLineV_Callback(hObject, eventdata, handles)
% Analyze | Line | Vertical | Electrons
sensorPlot(vcGetObject('sensor'),'electrons vline');
%OLD:  sensorPlotLine(vcGetObject('ISA'),'v','volts','space');
end

% --------------------------------------------------------------------
function menuAnLineH_Callback(hObject, eventdata, handles)
% Analyze | Line | Horizontal | Volts
sensorPlot(vcGetObject('sensor'),'volts hline');
end

% --------------------------------------------------------------------
function menuHorLineE_Callback(hObject, eventdata, handles)
% Analyze | Line | Horizontal | Electrons
sensorPlot(vcGetObject('sensor'),'electrons hline');
end

% --------------------------------------------------------------------
function menuVertLineE_Callback(hObject, eventdata, handles)
% Analyze | Line | Vertical | Electrons
sensorPlot(vcGetObject('sensor'),'electrons vline');
end

% --------------------------------------------------------------------
function menuHorLineDV_Callback(hObject, eventdata, handles)
% sensorPlotLine(vcGetObject('sensor'),'h','dv','space');
sensorPlot(vcGetObject('sensor'),'dv hline');
end

% --------------------------------------------------------------------
function menuVertLineDV_Callback(hObject, eventdata, handles)
% sensorPlotLine(vcGetObject('sensor'),'v','dv','space');
sensorPlot(vcGetObject('sensor'),'dv vline');
end

% --------------------------------------------------------------------
function menuFFThor_Callback(hObject, eventdata, handles)
sensorPlotLine(vcGetObject('sensor'),'h','volts','fft');
end

% --------------------------------------------------------------------
function menuFFTVert_Callback(hObject, eventdata, handles)
sensorPlotLine(vcGetObject('sensor'),'v','volts','fft');
end

% --------------------------------------------------------------------
function menuAnPixHistQ_Callback(hObject, eventdata, handles)
sensorPlotHistogram('e');
end

% --------------------------------------------------------------------
function menuAnPixHistV_Callback(hObject, eventdata, handles)
% Analyze | Line | Vertical | Volts
sensorPlot(vcGetObject('sensor'),'volts hist');
end

% --------------------------------------------------------------------
function menuAnSensorSNR_Callback(hObject, eventdata, handles)
%
sensorPlotSNR;
end

% --------------------------------------------------------------------
function menuAnROIStats_Callback(hObject, eventdata, handles)
% Analysis->ROI statistics
end

% --------------------------------------------------------------------
function menuAnBasicV_Callback(hObject, eventdata, handles)
%
sensorStats([],'basic','volts');
end

% --------------------------------------------------------------------
function menuAnBasicE_Callback(hObject, eventdata, handles)
sensorStats([],'basic','electrons');
end

% --------------------------------------------------------------------
function menuEditFontSize_Callback(~, ~, handles)
ieFontChangeSize(handles.coneMosaicWindow);
end

% --------------------------------------------------------------------
function menuHelp_Callback(hObject, eventdata, handles)
end

function menuAppNotes_Callback(hObject, eventdata, handles)
% Help | Documentation (web)
web('https://github.com/isetbio/isetbio/wiki','-browser');
end

function editKLMS_Callback(hObject, eventdata, handles)
% hObject    handle to editKLMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editKLMS as text
%        str2double(get(hObject,'String')) returns contents of editKLMS as a double
end

% --- Executes during object creation, after setting all properties.
function editKLMS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editKLMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes during object creation, after setting all properties.
function txtMosaic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMosaic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end



function editConeWidth_Callback(hObject, eventdata, handles)
% hObject    handle to editConeWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.cMosaic.pigment.width = 1e-6 * str2double(get(hObject, 'String'));
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --- Executes during object creation, after setting all properties.
function editConeWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editConeWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function editConeHeight_Callback(hObject, eventdata, handles)
% hObject    handle to editConeHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.pigment.height = 1e-6 * str2double(get(hObject, 'String'));
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function editConeHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editConeHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function editConeOpticalDensity_Callback(hObject, eventdata, handles)
% hObject    handle to editConeOpticalDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2double(get(hObject, 'String'));
assert(numel(val) == 3, 'invalid input for optical density');
handles.cMosaic.pigment.opticalDensity = val;
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function editConeOpticalDensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editConeOpticalDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function editMacularDensity_Callback(hObject, eventdata, handles)
% hObject    handle to editMacularDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.macular.density = str2double(get(hObject, 'String'));
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function editMacularDensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMacularDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function editConePeakEfficiency_Callback(hObject, eventdata, handles)
% hObject    handle to editConePeakEfficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2double(get(hObject, 'String'));
assert(numel(val) == 3, 'invalid input for peak efficiency');
handles.cMosaic.pigment.peakEfficiency = val;
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function editConePeakEfficiency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editConePeakEfficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --------------------------------------------------------------------
function menuPlotMacular_Callback(hObject, eventdata, handles)
% hObject    handle to menuPlotMacular (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

function menuPlotCone_Callback(hObject, eventdata, handles)
% hObject    handle to menuPlotCone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

function menuPlotMosaic_Callback(hObject, eventdata, handles)
% hObject    handle to menuPlotMosaic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

function menuPlotMosaicConeMosaic_Callback(~, ~, handles)
% hObject    handle to menuPlotMosaicConeMosaic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.plot('cone mosaic');
end

function menuPlotMosaicMeanAbsorptions_Callback(~, ~, handles)
% hObject    handle to menuPlotMosaicMeanAbsorptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.plot('mean absorptions');
end

function menuPlotMosaicMeanCurrent_Callback(~, ~, handles)
% hObject    handle to menuPlotMosaicMeanCurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.plot('mean current');
end

function menuPlotConeAbsorptance_Callback(~, ~, handles)
% hObject    handle to menuPlotConeAbsorptance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.plot('cone fundamentals');
end

function menuPlotMacularTransmittance_Callback(~, ~, handles)
% hObject    handle to menuPlotMacularTransmittance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.plot('macular transmittance');
end

function menuPlotMacularAbsorptance_Callback(~, ~, handles)
% hObject    handle to menuPlotMacularAbsorptance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.plot('macular absorptance');
end

function menuPlotMacularAbsorbance_Callback(~, ~, handles)
% hObject    handle to menuPlotMacularAbsorbance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.plot('macular absorbance');

end


% --------------------------------------------------------------------
function menuFileRefresh_Callback(hObject, eventdata, handles)
% hObject    handle to menuFileRefresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
coneMosaicGUIRefresh(hObject, eventdata, handles);
end


% --------------------------------------------------------------------
function menuPlotEMPath_Callback(~, ~, handles)
% hObject    handle to menuPlotEMPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.plot('eye movement path');
end


% --------------------------------------------------------------------
function menuEditGenerateEM_Callback(hObject, eventdata, handles)
% hObject    handle to menuEditGenerateEM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = ieReadString('Number of frames', '5000');
if ~isempty(str)
    handles.cMosaic.emGenSequence(str2double(str));
    handles.cMosaic.clearData();
    coneMosaicGUIRefresh(hObject, eventdata, handles);
end

end