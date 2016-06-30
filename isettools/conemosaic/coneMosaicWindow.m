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

% Last Modified by GUIDE v2.5 29-Jun-2016 16:13:08

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
%#ok<*ST2NM>

% check inputs
if isempty(varargin) || ~isa(varargin{1}, 'coneMosaic')
    error('cone mosaic object required');
end

% Choose default command line output for coneMosaicWindow
handles.output = hObject;
handles.cMosaic = varargin{1};
handles.mov = [];  % absorption movie
handles.curMov = [];  % photocurrent movie

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
set(handles.popupImageType, 'Value', 2); % mean absorptions
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
set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editCols_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editExpTime_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
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

function menuEditClearData_Callback(hObject, eventdata, handles)
handles.cMosaic.clearData();
handles.mov = [];
handles.curMov = [];
guidata(hObject, handles);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function editGam_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor', get(0,'defaultUicontrolBackgroundColor'));
end

function editGam_Callback(hObject, eventdata, handles)
handles.mov = [];
handles.curMov = [];
coneMosaicGUIRefresh(hObject,eventdata,handles);
end

function coneMosaicGUIRefresh(hObject, eventdata, handles)
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

% popup menu content
str = {'Cone mosaic'};
if ~isempty(cm.absorptions)
    str = [str {'Mean absorptions', 'Absorption movie'}];
end

if ~isempty(cm.current)
    str = [str {'Mean photocurrent', 'Photocurrent movie'}];
end

index = get(handles.popupImageType, 'Value');
if index > length(str), index = 1; end
plotType = str{index};
set(handles.popupImageType, 'Value', index);
set(handles.popupImageType, 'String', str);


switch plotType
    case 'Cone mosaic'
        % cone mosaic image
        resetMovieControl(handles);
        cm.plot('cone mosaic', 'hf', handles.axes2);
    case 'Mean absorptions'
        % mean cone absorptions
        resetMovieControl(handles);
        cm.plot('mean absorptions', 'hf', handles.axes2);
        
        % set up right click menu (context menu)
        c = uicontextmenu;
        handles.axes2.Children.UIContextMenu = c;
        uimenu(c, 'Label', 'hLine response', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'vLine response', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'hLine LMS', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'vLine LMS', 'Callback', @contextMenuPlot);
    case 'Absorption movie'
        set(handles.btnPlayPause, 'Visible', 'on');
        set(handles.sliderMovieProgress, 'Visible', 'on');
        if isempty(handles.mov)
            % generate movie
            [~, handles.mov] = cm.plot('absorptions', 'hf', 'none', ...
                'gamma', str2double(get(handles.editGam, 'String')));
            guidata(hObject, handles);
        end
        
        % set up right click menu (context menu)
        c = uicontextmenu;
        handles.axes2.Children.UIContextMenu = c;
        uimenu(c, 'Label', 'hLine response', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'vLine response', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'hLine LMS', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'vLine LMS', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'time series', 'Callback', @contextMenuPlot);
        
        % play movie
        btnPlayPause_Callback(hObject, eventdata, handles);
    case 'Mean photocurrent'
        resetMovieControl(handles);
        cm.plot('mean current', 'hf', handles.axes2);
        
        % set up right click menu (context menu)
        c = uicontextmenu;
        handles.axes2.Children.UIContextMenu = c;
        uimenu(c, 'Label', 'hLine response', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'vLine response', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'hLine LMS', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'vLine LMS', 'Callback', @contextMenuPlot);
    case 'Photocurrent movie'
        set(handles.btnPlayPause, 'Visible', 'on');
        set(handles.sliderMovieProgress, 'Visible', 'on');
        if isempty(handles.curMov) % generate movie for photocurrent
            [~, handles.curMov] = cm.plot('current', 'hf', 'none', ...
                'gamma', str2double(get(handles.editGam, 'String')));
            guidata(hObject, handles);
        end
        
        % set up right click menu (context menu)
        c = uicontextmenu;
        handles.axes2.Children.UIContextMenu = c;
        uimenu(c, 'Label', 'hLine response', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'vLine response', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'hLine LMS', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'vLine LMS', 'Callback', @contextMenuPlot);
        uimenu(c, 'Label', 'time series', 'Callback', @contextMenuPlot);
        
        % play movie
        btnPlayPause_Callback(hObject, eventdata, handles);
    otherwise
        error('Unknown plot type');
end

end

function contextMenuPlot(source, callbackdata)
% call back function for context menu
handles = guidata(source);

% determine which data to use (absorption or current)
contents = get(handles.popupImageType, 'String');
index = get(handles.popupImageType, 'Value');
if index > length(contents), index = 1; end
plotType = contents{index};

% get a point
[x, y] = ginput(1);

switch plotType
    case 'Mean absorptions'
        data = mean(handles.cMosaic.absorptions, 3);
        yStr = 'Absorptions';
    case 'Absorption movie'
        cnt = round(get(handles.sliderMovieProgress, 'Value'));
        if strcmp(source.Label, 'time series')
            data = handles.cMosaic.absorptions;
        else
            data = handles.cMosaic.absorptions(:, :, cnt);
        end
        
        % map x, y to cone positions
        x = x / size(handles.mov, 2) * size(data, 2);
        y = y / size(handles.mov, 1) * size(data, 1);
        yStr = 'Absorptions';
    case 'Mean photocurrent'
        data = mean(handles.cMosaic.current, 3);
        yStr = 'Photocurrent (pA)';
    case 'Photocurrent movie'
        cnt = round(get(handles.sliderMovieProgress, 'Value'));
        if strcmp(source.Label, 'time series')
            data = handles.cMosaic.current;
        else
            data = handles.cMosaic.current(:, :, cnt);
        end
        
        % map x, y to cone positions
        x = x / size(handles.curMov, 2) * size(data, 2);
        y = y / size(handles.curMov, 1) * size(data, 1);
        yStr = 'Photocurrent (pA)';
end
x = ieClip(round(x), 1, size(data, 2));
y = ieClip(round(y), 1, size(data, 1));

switch source.Label
    case 'hLine response'
        vcNewGraphWin; plot(data(y, :), 'LineWidth', 2); grid on;
        xlabel('Horionzal position (cones)'); ylabel(yStr);
    case 'vLine response'
        vcNewGraphWin; plot(data(:, x), 'LineWidth', 2); grid on;
        xlabel('Vertical position (cones)'); ylabel(yStr);
    case 'hLine LMS'
        vcNewGraphWin; names = 'LMS';
        for ii = 2 : 4 % L, M, S
            subplot(3, 1, ii-1);
            pos = find(handles.cMosaic.pattern(y, :) == ii);
            plot(pos, data(y, pos), '.-', 'LineWidth', 2); grid on;
            xlabel('Horizontal Position (cones');
            ylabel([names(ii-1) ' ' lower(yStr)]);
        end
    case 'vLine LMS'
        vcNewGraphWin; names = 'LMS';
        for ii = 2 : 4 % L, M, S
            subplot(3, 1, ii-1);
            pos = find(handles.cMosaic.pattern(:, x) == ii);
            plot(pos, data(pos, x), '.-', 'LineWidth', 2); grid on;
            xlabel('Vertical Position (cones');
            ylabel([names(ii-1) ' ' lower(yStr)]);
        end
    case 'time series'
        vcNewGraphWin;
        t = (1:size(data, 3)) * handles.cMosaic.sampleTime * 1e3;
        plot(t, squeeze(data(x, y, :)), 'LineWidth', 2);
        grid on; xlabel('Time (ms)'); ylabel(yStr);
    otherwise
        error('Unknown label type');
end

end

function resetMovieControl(handles)
% reset movie control
set(handles.btnPlayPause, 'Visible', 'off');
set(handles.sliderMovieProgress, 'Visible', 'off');

set(handles.btnPlayPause, 'Value', 1);
set(handles.sliderMovieProgress, 'Value', 1);
end

function menuPlot_Callback(hObject, eventdata, handles)
% Menu Plot
end

function menuEditFontSize_Callback(~, ~, handles)
ieFontChangeSize(handles.coneMosaicWindow);
end

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
density = str2num(get(handles.editKLMS, 'String'));
assert(numel(density) == 4, 'invalid input');

density = density / sum(density);
handles.cMosaic.spatialDensity = density;
menuEditClearData_Callback(hObject, eventdata, handles);
end

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

function txtMosaic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMosaic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end


function editConeWidth_Callback(hObject, eventdata, handles)
% hObject    handle to editConeWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

newWidth = 1e-6 * str2double(get(hObject, 'String'));

if handles.cMosaic.pigment.width ~= newWidth
    handles.cMosaic.pigment.width = newWidth;
    menuEditClearData_Callback(hObject, eventdata, handles);
end


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
newHeight = 1e-6 * str2double(get(hObject, 'String'));
if handles.cMosaic.pigment.height ~= newHeight
    handles.cMosaic.pigment.height = newHeight;
    menuEditClearData_Callback(hObject, eventdata, handles);
end
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
val = str2num(get(hObject, 'String'));
assert(numel(val) == 3, 'invalid input for optical density');

if any(handles.cMosaic.pigment.opticalDensity ~= val)
    handles.cMosaic.pigment.opticalDensity = val;
    menuEditClearData_Callback(hObject, eventdata, handles);
end
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
val = str2double(get(hObject, 'String'));
if handles.cMosaic.macular.density ~= val
    handles.cMosaic.macular.density = val;
    menuEditClearData_Callback(hObject, eventdata, handles);
end
end

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
val = str2num(get(hObject, 'String'));
assert(numel(val) == 3, 'invalid input for peak efficiency');
if any(handles.cMosaic.pigment.peakEfficiency ~= val)
    handles.cMosaic.pigment.peakEfficiency = val;
    menuEditClearData_Callback(hObject, eventdata, handles);
end
end

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

function menuFileRefresh_Callback(hObject, eventdata, handles)
% hObject    handle to menuFileRefresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

function menuPlotEMPath_Callback(~, ~, handles)
% hObject    handle to menuPlotEMPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.plot('eye movement path');
end

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

function popupImageType_Callback(hObject, eventdata, handles)
% hObject    handle to popupImageType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

function popupImageType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupImageType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function sliderMovieProgress_Callback(~, ~, handles)
% hObject    handle to sliderMovieProgress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index = get(handles.popupImageType, 'Value');
if index == 3  % absorption movie
    mov = handles.mov;
elseif index == 5  % photocurrent movie
    mov = handles.curMov;
end

cnt = round(get(handles.sliderMovieProgress, 'Value'));
assert(cnt <= size(mov, 4), 'slider choice out of range');
axes(handles.axes2); imshow(mov(:, :, :, cnt)); drawnow;

% register right click menu
c = uicontextmenu;
handles.axes2.Children.UIContextMenu = c;
uimenu(c, 'Label', 'hLine response', 'Callback', @contextMenuPlot);
uimenu(c, 'Label', 'vLine response', 'Callback', @contextMenuPlot);
uimenu(c, 'Label', 'hLine LMS', 'Callback', @contextMenuPlot);
uimenu(c, 'Label', 'vLine LMS', 'Callback', @contextMenuPlot);
uimenu(c, 'Label', 'time series', 'Callback', @contextMenuPlot);

end

function sliderMovieProgress_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderMovieProgress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function btnPlayPause_Callback(~, ~, handles)
% hObject    handle to btnPlayPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index = get(handles.popupImageType, 'Value');
if index == 3  % absorption movie
    mov = handles.mov;
elseif index == 5  % photocurrent movie
    mov = handles.curMov;
end

nFrames = size(mov, 4);
set(handles.sliderMovieProgress, 'Min', 1);
set(handles.sliderMovieProgress, 'Max', nFrames);
set(handles.sliderMovieProgress, 'SliderStep', [1/nFrames, 10/nFrames]);

if get(handles.btnPlayPause, 'Value')
    % play the video
    set(handles.btnPlayPause, 'String', 'Pause');
    cnt = round(get(handles.sliderMovieProgress, 'Value'));
    axes(handles.axes2);
    while get(handles.btnPlayPause, 'Value')
        imshow(mov(:, :, :, cnt));
        set(handles.sliderMovieProgress, 'Value', cnt);
        
        drawnow; cnt = cnt + 1;
        if cnt > nFrames, cnt = 1; end
    end
else
    % pause video
    set(handles.btnPlayPause, 'String', 'Play');
    
    % register right click menu
    c = uicontextmenu;
    handles.axes2.Children.UIContextMenu = c;
    uimenu(c, 'Label', 'hLine response', 'Callback', @contextMenuPlot);
    uimenu(c, 'Label', 'vLine response', 'Callback', @contextMenuPlot);
    uimenu(c, 'Label', 'hLine LMS', 'Callback', @contextMenuPlot);
    uimenu(c, 'Label', 'vLine LMS', 'Callback', @contextMenuPlot);
    uimenu(c, 'Label', 'time series', 'Callback', @contextMenuPlot);
end

end
