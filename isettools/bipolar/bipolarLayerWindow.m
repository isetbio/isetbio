function varargout = bipolarLayerWindow(varargin)
% BIPOLARLAYERWINDOW MATLAB code for bipolarlayerwindow.fig
%      BIPOLARLAYERWINDOW, by itself, creates a new BIPOLARLAYERWINDOW or raises the existing
%      singleton*.
%
%      H = BIPOLARLAYERWINDOW returns the handle to a new BIPOLARLAYERWINDOW or the handle to
%      the existing singleton*.
%
%      BIPOLARLAYERWINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BIPOLARLAYERWINDOW.M with the given input arguments.
%
%      BIPOLARLAYERWINDOW('Property','Value',...) creates a new BIPOLARLAYERWINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bipolarLayerWindow_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bipolarLayerWindow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bipolarlayerwindow

% Last Modified by GUIDE v2.5 13-Jun-2017 11:35:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bipolarLayerWindow_OpeningFcn, ...
                   'gui_OutputFcn',  @bipolarLayerWindow_OutputFcn, ...
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

% --- Executes just before bipolarlayerwindow is made visible.
function bipolarLayerWindow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bipolarlayerwindow (see VARARGIN)

p = inputParser;
vFunc = @(x)(isequal(class(x),'bipolarLayer'));
p.addRequired('bipolar',vFunc);

% Check that we have the bipolar layer
p.parse(varargin{:});

bp = varargin{1};
bp.fig = hObject;   % Store this figure handle

% Choose default command line output for bipolarlayerwindow
handles.output = hObject;

% This is the bipolar layer object
handles.bipolar = bp;

% Store the cell types in the list box.  Maybe we should have a 'name' slot
% in the mosaic and store that, where the default name is the cell type.
mosaicNames = cell(1,length(bp.mosaic));
for ii=1:length(bp.mosaic)
    mosaicNames{ii} = bp.mosaic{ii}.cellType;  
end
set(handles.listMosaics,'String',mosaicNames);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bipolarlayerwindow wait for user response (see UIRESUME)
% uiwait(handles.bipolarlayerwindow);

% Refresh/Initialize window information
bipolarWindowRefresh(handles);

% Very important for good rendering speed
set(hObject, 'Renderer', 'OpenGL')

end

% --- Outputs from this function are returned to the command line.
function varargout = bipolarLayerWindow_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)
% hObject    handle to menuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function menuEdit_Callback(hObject, eventdata, handles)
% hObject    handle to menuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.linearMov = [];
% handles.psthMov = [];
end

% --------------------------------------------------------------------
function menuPlot_Callback(hObject, eventdata, handles)
% hObject    handle to menuPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function menuMosaic_Callback(hObject, eventdata, handles)
% hObject    handle to menuMosaic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function menuAnalyze_Callback(hObject, eventdata, handles)
% hObject    handle to menuAnalyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on selection change in popupResponseSelect.
function popupResponseSelect_Callback(hObject, eventdata, handles)
% Popup over main response window
%
% hObject    handle to popupResponseSelect (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupResponseSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupResponseSelect

bipolarWindowRefresh(handles);

end

% --- Executes during object creation, after setting all properties.
function popupResponseSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupResponseSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --------------------------------------------------------------------
function menuPlotCTS_Callback(hObject, eventdata, handles)
% Plot | Current time series
%

% We should make this a function that selects a point and plots it on the
% window axis.

% The current time series is at a point. Get the point
[x, y] = ginput(1); % Rounded and clipped to the bipolar mosaic size

sz = size(handles.bipolar.responseCenter);
xlim = get(gca,'xlim'); ylim = get(gca,'ylim');
pos(1) = ieClip(round(x - xlim(1)), 1, sz(2));
pos(2) = ieClip(round(y - ylim(1)), 1, sz(1));

% Draw a circle around the selected point.
viscircles([x,y],0.7);

% pos is in units of the matrix of data, not microns on the surface.
% Deal with this.
handles.bipolar.plot('response time series','pos',pos);

end

% --------------------------------------------------------------------
function menuPlotSpatialRFMosaic_Callback(hObject, eventdata, handles)
% Plot | Spatial RF mosaic
%
% hObject    handle to menuPlotSpatialRFMosaic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.bipolar.plot('mosaic');

end


% --------------------------------------------------------------------
function menuPlotSpatialRF_Callback(hObject, eventdata, handles)
%
% hObject    handle to menuPlotSpatialRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Single spatial RF')

end

% --------------------------------------------------------------------
function menuPlotTemporalIR_Callback(hObject, eventdata, handles)
% hObject    handle to menuPlotTemporalIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Temporal IRF')

end

% --------------------------------------------------------------------
function menFileSave_Callback(hObject, eventdata, handles)
% File | Save
%
% hObject    handle to menFileSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Save.  NYI')
end

% --------------------------------------------------------------------
function menuFileRefresh_Callback(hObject, eventdata, handles)
% hObject    handle to menuFileRefresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bipolarWindowRefresh(handles)
end

% --------------------------------------------------------------------
function menuFileClose_Callback(hObject, eventdata, handles)
% File | Close
% hObject    handle to menuFileClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.bp.fig = [];
delete(handles.bipolarWindow);
end

%% Internal functions

function bipolarWindowRefresh(handles)
% Update all the text fields and such with the data in the mosaic

bp    = handles.bipolar;
fig   = figure(bp.fig);
gdata = guidata(fig);

% Show the appropriate response axis plot
axes(gdata.axisResponse);
cla(gdata.axisResponse,'reset');

% Selected string in the popup
contents = cellstr(get(gdata.popupResponseSelect,'String'));
str = contents{get(gdata.popupResponseSelect,'Value')};

g = str2double(get(handles.editGamma,'string'));

% Get the selected mosaic from the listbox properly ...
nMosaic = get(gdata.listMosaics,'Value');
% Need to deal with interface for layer, not mosaic

% Update the main axis window with the relevant plot
switch(str)
    
    case 'Bipolar mosaic'
        ieInWindowMessage('Building mosaic',handles);
        gdata.bipolar.plot('mosaic','nMosaic',nMosaic);
        ieInWindowMessage('',handles);

    case 'Bipolar mean (image)'
        gdata.bipolar.plot('response image','gamma',g,'nMosaic',nMosaic);
        colorbar;
        
    case 'Bipolar movie'
        ieInWindowMessage('Showing movie',handles);
        gdata.bipolar.plot('Response movie','gamma',g,'nMosaic',nMosaic);
        ieInWindowMessage('',handles);

    otherwise
        error('Unknown plot type %s\n',str);
end

% Text description
str = bp.describe;
set(gdata.txtBipolarProperties,'string',str);

end


% --- Executes on button press in btnPlayPause.
function btnPlayPause_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to btnPlayPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of btnPlayPause
end

% --- Executes on slider movement.
function sliderMovieProgress_Callback(hObject, eventdata, handles)
% hObject    handle to sliderMovieProgress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
end

% --- Executes during object creation, after setting all properties.
function sliderMovieProgress_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderMovieProgress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes during object creation, after setting all properties.
function editGamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
% --- Edit gamma for display.
function editGamma_Callback(hObject, eventdata, handles)
% hObject    handle to editGamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGamma as text
%        str2double(get(hObject,'String')) returns contents of editGamma as a double
bipolarWindowRefresh(handles)
%
end

%% Layer management

% --- Executes on selection change in listMosaics.
function listMosaics_Callback(hObject, eventdata, handles)
% Select a mosaic from the listbox.

% User selected a new mosaic.  Refresh.
bipolarWindowRefresh(handles);

end

% --- Executes during object creation, after setting all properties.
function listMosaics_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listMosaics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
