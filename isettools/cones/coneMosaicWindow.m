function varargout = coneMosaicWindow(varargin)
%Cone image coneMosaicWindow interface
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

% Last Modified by GUIDE v2.5 08-Jun-2017 21:21:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
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
handles.cMosaic.hdl = hObject;

% Adjust the database and bring this figure to the front
vcSetFigureHandles('conemosaic',hObject,eventdata,handles);
figure(hObject);

% Set the popup default image selection to mean absorptions when the window
% opens.
str = get(handles.popupImageType, 'String');
if iscell(str) && length(str) > 1
    % This is mean absorptions
    set(handles.popupImageType, 'Value',2);  
end

% Refresh and move on
coneMosaicGUIRefresh(hObject, eventdata, handles);

% Set the font size based on the ISETBIO preferences
ieFontInit(hObject);

% Very important for good rendering speed
set(hObject, 'Renderer', 'OpenGL')

end

% --- Outputs from this function are returned to the command line.
function varargout = coneMosaicWindow_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
end

function btnComputeImage_Callback(hObject, eventdata, handles)
% Computes the image from the optics data - button at the bottom
oi = vcGetObject('OI');
if isempty(oi) || isempty(oiGet(oi, 'photons'))
    warning('No optical image.  Use ieAddObject(oi) to store.');
    return;
end

handles.cMosaic.compute(oi);
handles.cMosaic.name = oiGet(oi,'name');
set(handles.popupImageType, 'Value', 2); % mean absorptions

coneMosaicGUIRefresh(hObject, eventdata, handles);

end

function menuAnComputeFromOI_Callback(hObject, eventdata, handles)
% Cones | Compute absorptions 
% Computes from an OI in the database
btnComputeImage_Callback(hObject, eventdata, handles);
end

% Edit box - adjust number of rows
function editRows_Callback(hObject, eventdata, handles)
% Columns text box
handles.cMosaic.rows = str2double(get(hObject, 'String'));
menuEditClearData_Callback(hObject, eventdata, handles)
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% Edit box - adjust number of columns
function editCols_Callback(hObject, eventdata, handles)
% Columns text box
handles.cMosaic.cols = str2double(get(hObject, 'String'));
menuEditClearData_Callback(hObject, eventdata, handles)
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
delete(handles.coneMosaicWindow);
end

function menuEdit_Callback(hObject, eventdata, handles)
end

function menuEditName_Callback(hObject, eventdata, handles)
% Edit | Rename
str = ieReadString('New name', handles.cMosaic.name);
if ~isempty(str), handles.cMosaic.name = str; end
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

function menuEditClearData_Callback(hObject, eventdata, handles)
% Edit | Clear data
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
% Text book Gamm

set(handles.editGam,'value',str2double(get(handles.editGam,'string')));
coneMosaicGUIRefresh(hObject,eventdata,handles);

end

function coneMosaicGUIRefresh(hObject, eventdata, handles)
% Update the cone mosaic window interface - main pulldown
%
%   coneMosaciGUIRefresh(handles)
%
% HJ/BW, ISETBIO TEAM, 2016

% get coneMosaic object
cm = handles.cMosaic;

% Place name in text string box
set(handles.txtName,'string',sprintf('%s',cm.name));

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

% set outersegment description
str = cm.descriptionOS;
set(handles.txtOS,'string',str);

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
if index > length(str), index = 1; end  % Mean absorptions
plotType = str{index};
set(handles.popupImageType, 'Value', index);
set(handles.popupImageType, 'String', str);

%% Here are the different window options

switch plotType
    case 'Cone mosaic'
        % cone mosaic image
        % TODO:  For large mosaics, the computation is slow.  We should
        % compute it once and store it.
        cm.plot('cone mosaic', 'hf', handles.axes2);

        resetMovieControl(handles,'off');
        
        % enable plot options in menu
        enable.hLine = 'off';    enable.vLine = 'off';
        enable.hLineLMS = 'off'; enable.vLineLMS = 'off';
        enable.timeSeries = 'off';
        contextMenuEnable(handles,enable);


    case 'Mean absorptions'
        % mean cone absorptions
        cm.plot('mean absorptions', 'hf', handles.axes2);
        axis image
        resetMovieControl(handles,'off');

        % Why isn't axis image in the plot() routine?
        
        % enable plot options in menu
        enable.hLine = 'on';    enable.vLine = 'on';
        enable.hLineLMS = 'on'; enable.vLineLMS = 'on';
        enable.timeSeries = 'off';
        contextMenuInit(handles);
        contextMenuEnable(handles,enable);
        
    case 'Absorption movie'
        resetMovieControl(handles,'on');  %Turn off the movie controller

        enable.hLine = 'on';    enable.vLine = 'on';
        enable.hLineLMS = 'on'; enable.vLineLMS = 'on';
        enable.timeSeries = 'on';
        contextMenuEnable(handles,enable);
        
        % Initiate movie display GUI elements.
        set(handles.btnPlayPause, 'Visible', 'on');
        set(handles.btnPlayPause, 'Value', 1);  % Auto start the movie
        set(handles.sliderMovieProgress, 'Visible', 'off');
        
        % play movie if more than one frame
        btnPlayPause_Callback(hObject, eventdata, handles);
        
    case 'Mean photocurrent'
        cm.plot('mean current', 'hf', handles.axes2);
        axis image
        
        resetMovieControl(handles,'off');  %Turn off the movie controller
        
        enable.hLine = 'on';    enable.vLine = 'on';
        enable.hLineLMS = 'on'; enable.vLineLMS = 'on';
        enable.timeSeries = 'off';
        contextMenuInit(handles);
        contextMenuEnable(handles,enable);

    case 'Photocurrent movie'
        resetMovieControl(handles,'on');  %Turn off the movie controller

        % Graphics elements
        set(handles.btnPlayPause, 'Visible', 'on');
        set(handles.btnPlayPause, 'Value', 1);  % Auto start the movie
        set(handles.sliderMovieProgress, 'Visible', 'off');
        
        % enable plot options in menu
        enable.hLine = 'on';    enable.vLine = 'on';
        enable.hLineLMS = 'on'; enable.vLineLMS = 'on';
        enable.timeSeries = 'on';
        contextMenuInit(handles);
        contextMenuEnable(handles,enable);
        
        % play movie
        btnPlayPause_Callback(hObject, eventdata, handles);
    otherwise
        error('Unknown plot type');
end


end

function c = contextMenuInit(handles)
% Set up right click menu (context menu)
% 
%  Typical sequence is 
%    Set up enable.XXX
%    c = contextMenuInit(handles,enable)
%    contextMenuEnable(enable)
%

c = uicontextmenu;
if ~isempty(handles.axes2.Children)
    for ichild = 1:size(handles.axes2.Children,1)
        handles.axes2.Children(ichild).UIContextMenu = c;
    end
    uimenu(c, 'Label', 'hLine response', 'Callback', @contextMenuPlot);
    uimenu(c, 'Label', 'vLine response', 'Callback', @contextMenuPlot);
    uimenu(c, 'Label', 'hLine LMS', 'Callback', @contextMenuPlot);
    uimenu(c, 'Label', 'vLine LMS', 'Callback', @contextMenuPlot);
end

end

function contextMenuEnable(handles,enable)
% enable plot options in menu
%
%    c = contextMenuInit(handles,enable)
%    % Set up enable.XXX
%    contextMenuEnable(enable)
%
% enable is a structure of 'on' and 'off' values
%
%   enable.hLine
%   enable.hLineLMS
%   enable.vLine
%   enable.vLineLMS
%   enable.timeSeries

set(handles.menuPlotHLine, 'Enable', enable.hLine);
set(handles.menuPlotVLine, 'Enable', enable.vLine);
set(handles.menuPlotHLineLMS, 'Enable',enable.hLineLMS);
set(handles.menuPlotVLineLMS, 'Enable', enable.vLineLMS);
set(handles.menuPlotTimeSeries, 'Enable', enable.timeSeries);

end

function contextMenuPlot(source, callbackdata)
% call back function for context menu
%
% Source is a Menu object 
% The guidata of source contains all the gui objects
handles = guidata(source);

% determine which data to use (absorption or current)
contents = get(handles.popupImageType, 'String');
index    = get(handles.popupImageType, 'Value');
if index > length(contents), index = 1; end
plotType = contents{index};

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
        yStr = 'Photocurrent (pA)';
end

% The plots below are with respect to a point.
% Get the point
[x, y] = ginput(1); % Rounded and clipped to the data
x = ieClip(round(x), 1, size(data, 2));
y = ieClip(round(y), 1, size(data, 1));

% Draw a circle around the selected point.
viscircles([x,y],0.7);

switch source.Label
    case 'hLine response'
        vcNewGraphWin; plot(data(y, :), 'LineWidth', 2); 
        grid on; xlabel('Horizontal position (cones)'); ylabel(yStr);
        set(gca,'userdata',data(y,:));
        
    case 'vLine response'
        vcNewGraphWin; plot(data(:, x), 'LineWidth', 2); 
        grid on; xlabel('Vertical position (cones)'); ylabel(yStr);
        set(gca,'userdata',data(:,x));
        
    case 'hLine LMS'
        % Save the work more completely in the window, please!
        vcNewGraphWin([],'tall'); names = 'LMS';
        c = {'ro-','go-','bo-'};
        for ii = 2 : 4 % L, M, S
            subplot(3, 1, ii-1);
            pos = find(handles.cMosaic.pattern(y, :) == ii);
            plot(pos, data(y, pos), c{ii-1}, 'LineWidth', 2); grid on;
            uData.pos{ii-1} = pos; uData.data{ii-1}=data(y,pos);
            xlabel('Horizontal Position (cones');
            ylabel([names(ii-1) ' ' yStr]);
            set(gca,'xlim',[1 size(data,2)]);
        end
        set(gca,'userdata',uData);
        
    case 'vLine LMS'
        % Save the work more completely in the window, please!
        vcNewGraphWin([],'tall'); names = 'LMS';
        c = {'ro-','go-','bo-'};
        for ii = 2 : 4 % L, M, S
            subplot(3, 1, ii-1);
            pos = find(handles.cMosaic.pattern(:, x) == ii);
            plot(pos, data(pos, x), c{ii-1}, 'LineWidth', 2); grid on;
            uData.pos{ii-1} = pos; uData.data{ii-1}=data(pos,x);
            xlabel('Vertical Position (cones');
            ylabel([names(ii-1) ' ' yStr]);
            set(gca,'xlim',[1 size(data,1)]);
        end
        set(gca,'userdata',uData);
        
    case 'time series'
        % Time series is enabled for the absorption and current movie modes
        vcNewGraphWin;
        if index == 3  %  absorption movie
            % Show the absorption time series
            mx = max(handles.cMosaic.absorptions(:));
            mn = min(handles.cMosaic.absorptions(:));
            t = (1:size(data, 3)) * handles.cMosaic.integrationTime * 1e3;
        elseif index == 5  % photocurrent movie
            mx = max(handles.cMosaic.current(:));
            mn = min(handles.cMosaic.current(:));
            t = (1:size(data, 3)) * handles.cMosaic.integrationTime * 1e3;
        end
        plot(t, squeeze(data(y, x, :)), 'LineWidth', 2);
        uData.timerseries = t;
        uData.x = x; uData.y = y;
        grid on; xlabel('Time (ms)'); ylabel(yStr);
        set(gca,'ylim',[mn mx]);
        set(gca,'userdata',uData);
        
    otherwise
        error('Unknown label type');
end

end

function resetMovieControl(handles,status)
% reset movie controls for playing a movie

switch status
    case 'off'
        set(handles.btnPlayPause, 'Visible', 'off');
        set(handles.btnPlayPause, 'Value', 0);  % Pause the movie

        set(handles.sliderMovieProgress, 'Visible', 'off');
        set(handles.sliderMovieProgress, 'Value', 1);
    case 'on'
        set(handles.btnPlayPause, 'Visible', 'on');
        set(handles.btnPlayPause, 'Value', 1);  % Auto start the movie

        set(handles.sliderMovieProgress, 'Visible', 'off');
        set(handles.sliderMovieProgress, 'Value', 1);

end
        
end

function menuPlot_Callback(hObject, eventdata, handles)
% Menu Plot
end

function menuEditFontSize_Callback(~, ~, handles)
ieFontSizeSet(handles.coneMosaicWindow);
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

if any(handles.cMosaic.pigment.opticalDensity(:) ~= val(:))
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
if any(handles.cMosaic.pigment.peakEfficiency(:) ~= val(:))
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
% Top  level 
%   Mosaic 
% hObject    handle to menuPlotMosaic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

function menuPlotMosaicConeMosaic_Callback(~, ~, handles)
% Plot | Mosaic | Cone Mosaic
%
% hObject    handle to menuPlotMosaicConeMosaic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.plot('cone mosaic',...
    'showCorrespondingRectangularMosaicInstead',false);
end

function menuPlotMosaicMeanAbsorptions_Callback(~, ~, handles)
% Plot | Mosaic | Mean absorptions
%
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
% Plot | Macular | Transmittance
%
% hObject    handle to menuPlotMacularTransmittance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cMosaic.plot('macular transmittance');
end

function menuPlotMacularAbsorptance_Callback(~, ~, handles)
% Plot | Macular | Absorptance
%
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


% --------------------------------------------------------------------
function menuCones_Callback(hObject, eventdata, handles)
% Cones - Main Pull down for computing.
end

function menuConesGenerateEM_Callback(hObject, eventdata, handles)
% Cones | Generate eye movements
%
str = ieReadString('Number of frames', '500');
if ~isempty(str)
    handles.cMosaic.emGenSequence(str2double(str));
    menuEditClearData_Callback(hObject, eventdata, handles);
    set(handles.popupImageType, 'Value', 2); % mean absorptions
    coneMosaicGUIRefresh(hObject, eventdata, handles);
end

end

% --------------------------------------------------------------------
function menuConesAbsorptions_Callback(hObject, eventdata, handles)
% Cones | Compute absorptions
%
% Loads current oi to compute the absorptions.  If no oi is selected, it
% complains.

oi = vcGetObject('OI');
if isempty(oi) || isempty(oiGet(oi, 'photons'))
    warning('No optical image.  Use ieAddObject(oi) to store.');
    return;
end

fprintf('Calculating with optical image %s\n',oiGet(oi,'name'));
handles.cMosaic.compute(oi);
handles.cMosaic.name = oiGet(oi,'name');
set(handles.popupImageType, 'Value', 2); % mean absorptions
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuConesPhotocurrent_Callback(hObject, eventdata, handles)
% Cones | Compute photocurrent
%
handles.cMosaic.computeCurrent;
set(handles.popupImageType, 'Value', 4); % mean current
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuConePhotocurrentNoise_Callback(hObject, eventdata, handles)
% Cones | Toggle photocurrent noise
% Also executes computeCurrent

set(handles.btnPlayPause,'Value',0);  % Turn off any movie.

% Flip from whatever state to the other
switch handles.cMosaic.os.noiseFlag
    case 'random'
        handles.cMosaic.os.noiseFlag = 'none';
    case 'frozen'
        handles.cMosaic.os.noiseFlag = 'random';
    case 'none'
        handles.cMosaic.os.noiseFlag = 'frozen';
end

% We used to use this.  But now, I think it should always be off.
handles.menuConePhotocurrentNoise.Checked = 'off';
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

%---------------

function popupImageType_Callback(hObject, eventdata, handles)
% hObject    handle to popupImageType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Stop any movies.
set(handles.btnPlayPause,'Value',0);  % Turn off the movie.

% Refresh.
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
if index == 3,     mov = handles.cMosaic.absorptions;
elseif index == 5, mov = handles.cMosaic.current;
end

gam = get(handles.editGam,'value');

mov = ieScale(mov,0,1) .^ gam;
mind = min(mov(:)); maxd = max(mov(:));

cnt = round(get(handles.sliderMovieProgress, 'Value'));
assert(cnt <= size(mov, ndims(mov)), 'slider choice out of range');
axes(handles.axes2); 
imagesc(mov(:,:,cnt)); 
axis image; set(gca,'xticklabel','','yticklabel',''); caxis([mind maxd]); 
drawnow;
set(handles.txtMovieFrame,'string',cnt);

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
% Play/Pause button near slider

% We don't want the numbers during the movie, I think.
axis off;

% Which type of data, absorptions or current
index = get(handles.popupImageType, 'Value');
if index == 3,       mov = handles.cMosaic.absorptions;
elseif index == 5,   mov = handles.cMosaic.current;
end

% Gamma for display
gam = get(handles.editGam,'value');

% Check that there are some frames for a movie
if ismatrix(mov), nFrames = 1;
else              nFrames = size(mov, ndims(mov));
end
if nFrames == 1, 
    str = sprintf('Only one frame. No movie to show.'); 
    ieInWindowMessage(str,handles,3);
    return; 
end

% Display up the slider
set(handles.sliderMovieProgress, 'Min', 1);
set(handles.sliderMovieProgress, 'Max', nFrames);
set(handles.sliderMovieProgress, 'SliderStep', [1, round(max(nFrames/4,1))]);  %Minor step and major step

% If play, then play.
if get(handles.btnPlayPause, 'Value')
    % play video if  value is not zero
    set(handles.btnPlayPause, 'String', 'Pause');
    set(handles.sliderMovieProgress, 'Visible', 'off');
    set(handles.txtMovieFrame,'Visible','off');

    gData = guidata(handles.coneMosaicWindow);
    axes(gData.axes2);
    while get(handles.btnPlayPause, 'Value')
        % Plays the movie until the pause button is pushed.
        % I removed the slider and counter for now.  We could pass the pause button into
        % the ieMovie routine.  Say, ('hdlCounter',XXX,'hdlPause',yyy)
        ieMovie(mov,'gamma',gam); 
        pause(0.2);  % The brief pause identifies the end of the movie.
    end
    
else
    % pause video when the value is false (zero)
    set(handles.btnPlayPause, 'String', 'Play');
    
    % Turn on the slider and frame counter
    set(handles.txtMovieFrame,'Visible','on');
    val = str2double(get(handles.txtMovieFrame,'String'));
    val = max(1,val);
    set(handles.sliderMovieProgress, 'Value', val);
    set(handles.sliderMovieProgress, 'Visible', 'on');
    
    % We should show the right image frame here
    % TODO
    
end

% register right click menu
c = uicontextmenu;
for ichild = 1:size(handles.axes2.Children,1)
    handles.axes2.Children(ichild).UIContextMenu = c;
end
uimenu(c, 'Label', 'hLine response', 'Callback', @contextMenuPlot);
uimenu(c, 'Label', 'vLine response', 'Callback', @contextMenuPlot);
uimenu(c, 'Label', 'hLine LMS', 'Callback', @contextMenuPlot);
uimenu(c, 'Label', 'vLine LMS', 'Callback', @contextMenuPlot);
uimenu(c, 'Label', 'time series', 'Callback', @contextMenuPlot);

end

function editEccentricity_Callback(hObject, eventdata, handles)
% hObject    handle to editEccentricity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEccentricity as text
%        str2double(get(hObject,'String')) returns contents of editEccentricity as a double
end

% --- Executes during object creation, after setting all properties.
function editEccentricity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEccentricity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function menuPlotHLine_Callback(hObject, eventdata, handles)
contextMenuPlot(hObject, []);
end

function menuPlotVLine_Callback(hObject, eventdata, handles)
contextMenuPlot(hObject, []);
end

function menuPlotHLineLMS_Callback(hObject, eventdata, handles)
contextMenuPlot(hObject, []);
end

function menuPlotVLineLMS_Callback(hObject, eventdata, handles)
contextMenuPlot(hObject, []);
end

function menuPlotTimeSeries_Callback(hObject, ~, handles)
set(handles.btnPlayPause, 'Value', 0);  % Pause the movie
contextMenuPlot(hObject, []);
end
