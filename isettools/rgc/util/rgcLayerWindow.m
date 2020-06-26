function varargout = rgcLayerWindow(varargin)
% MATLAB code for rgcLayerWindow.fig
%
% Syntax:
%   [varargout] = rgcLayerWindow([varargin])
%
% Description:
%    RGCLAYERWINDOW, by itself, creates a new RGCLAYERWINDOW or raises the
%    existing singleton*.
%
%    H = RGCLAYERWINDOW returns the handle to a new RGCLAYERWINDOW or the
%    handle to the existing singleton*.
%
%    RGCLAYERWINDOW('CALLBACK', hObject, eventData, handles, ...) calls the
%    local function named CALLBACK in RGCLAYERWINDOW.M using the given
%    input arguments.
%
%    RGCLAYERWINDOW('Property', 'Value', ...) creates a new RGCLAYERWINDOW
%    or raises the existing singleton*. Starting from the left, property
%    value pairs are applied to the GUI before rgcLayerWindow_OpeningFcn
%    gets called. An unrecognized property name or invalid value makes
%    property application stop. All inputs are passed to
%    rgcLayerWindow_OpeningFcn via varargin.
%
%    * See GUI Options on GUIDE's Tools menu. Choose "GUI allows only one
%      instance to run (singleton)".
%
% Inputs:
%    None.
%
% Outputs:
%    None required.
%
% Optional key/value pairs:
%    Needs to be completed
%
% Notes:
%    * TODO: Edit the above text to modify the response to help better fit
%      to rgcLayerWindow.
%
% See Also:
%   GUIDE, GUIDATA, GUIHANDLES
%

% History:
%    07/15/17  GUIDE  Last Modified by GUIDE v2.5 15-Jul-2017 21:05:28
%    06/03/19  JNM    Documentation pass

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton', gui_Singleton, ...
    'gui_OpeningFcn', @rgcLayerWindow_OpeningFcn, ...
    'gui_OutputFcn', @rgcLayerWindow_OutputFcn, ...
    'gui_LayoutFcn', [] , 'gui_Callback', []);
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

% --- Executes just before rgcLayerWindow is made visible.
function rgcLayerWindow_OpeningFcn(hObject, eventdata, handles, varargin)
% Open the layer window (prior to it being made visible)
%
% Syntax:
%   rgcLayerWindow_OpeningFcn(hObject, eventdata, handles, [varargin])
%
% Description:
%    Open the rgc layer window prior to making it visible.
%
%    This function has no output args, see OutputFcn.
%
% Inputs:
%    hObject   - Handle. The handle to the figure.
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. A structure with handles & user data (see GUIDATA)
%    varargin  - (Optional) VARIES. The command line arguments to
%                rgcLayerWindow (see VARARGIN).
%

p = inputParser;
vFunc = @(x)(isequal(class(x), 'rgcLayer'));
p.addRequired('rgc', vFunc);

% Check that we have the bipolar layer
p.parse(varargin{:});

rgcL = varargin{1};
rgcL.fig = hObject;   % Store this figure handle

% Choose default command line output for bipolarlayerwindow
handles.output = hObject;
handles.rgcLayer = varargin{1};
handles.spikesMovie = [];  % spike movie

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rgcLayerWindow wait for user response (see UIRESUME)
% uiwait(handles.rgcLayerWindow);

mosaicNames = cell(1, length(rgcL.mosaic));
for ii = 1:length(rgcL.mosaic)
    mosaicNames{ii} = rgcL.mosaic{ii}.cellType;  
end
set(handles.listMosaics, 'String', mosaicNames);

% Refresh/Initialize window information
rgcLayerWindowRefresh(handles);

% Very important for good rendering speed
set(hObject, 'Renderer', 'OpenGL')

handles.linearMov = [];
handles.psthMov = [];

end

% --- Outputs from this function are returned to the command line.
function varargout = rgcLayerWindow_OutputFcn(hObject, eventdata, handles)
% Output function for the rgc layer window.
%
% Syntax:
%   [varargout] = rgcLayerWindow_OutputFcn(hObject, eventdata, handles)
%
% Description:
%    Output function for the rgc layer window.
%
% Inputs:
%    hObject   - Handle. The handle to the figure.
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    varargout - Cell. cell array for returning output args (see VARARGOUT)
%
% Optional key/value pairs:
%    None.
%

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)
% (File) Call function to the menu file object.
%
% Syntax:
%   menuFile_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function to the File object on the menus.
%
% Inputs:
%    hObject   - Handle. The handle to menuFile (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

% --------------------------------------------------------------------
function menuEdit_Callback(hObject, eventdata, handles)
% (Edit) Call function to the menu edit object.
%
% Syntax:
%   menuEdit_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function to the Edit object on the menus.
%
% Inputs:
%    hObject   - Handle. The handle to menuEdit (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

handles.linearMov = [];
handles.psthMov = [];
end

% --------------------------------------------------------------------
function menuPlot_Callback(hObject, eventdata, handles)
% (Plot) Call function to the menu Plot object.
%
% Syntax:
%   menuPlot_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function to the Plot object on the menus.
%
% Inputs:
%    hObject   - Handle. The handle to menuPlot (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

% --------------------------------------------------------------------
function menuMosaic_Callback(hObject, eventdata, handles)
% (Mosaic) Call function to the menu Mosaic object.
%
% Syntax:
%   menuMosaic_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function to the Mosaic menu object.
%
% Inputs:
%    hObject   - Handle. The handle to menuMosaic (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

end

% --------------------------------------------------------------------
function menuAnalyze_Callback(hObject, eventdata, handles)
% (Analyze) Call function to the menu Analyze object.
%
% Syntax:
%   menuAnalyze_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function to the Analyze menu object.
%
% Inputs:
%    hObject   - Handle. The handle to menuAnalyze (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

end

% --- Executes on selection change in popupResponseSelect.
function popupResponseSelect_Callback(hObject, eventdata, handles)
% Popup over main response window
%
% Syntax:
%   popupResponseSelect_Callback(hObject, eventdata, handles)
%
% Description:
%    A popup over the main response window, which executes on a selection
%    change in the window.
%
% Inputs:
%    hObject   - Handle. The handle to popupResponseSelect (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hints:
%        contents = cellstr(get(hObject, 'String'))  - returns
%           popupResponseSelect contents as a cell array
%        contents{get(hObject, 'Value')} - returns the selected item from
%           popupResponseSelect object
%

rgcLayerWindowRefresh(handles);

end

% --- Executes during object creation, after setting all properties.
function popupResponseSelect_CreateFcn(hObject, eventdata, handles)
% Popup object creation using all defined properties.
%
% Syntax:
%   popupResponseSelect_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Create the object for popup response select using defined properties.
%
% Inputs
%    hObject   - Handle. The handle to popupResponseSelect (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty. Handles are not created until all of the CreateFcns
%                have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

end

% Plot Menu
% --------------------------------------------------------------------
function menuPlotPSTH_Callback(hObject, eventdata, handles)
% (Plot | PSTH) Call the PSTH menu option under Plot
%
% Syntax:
%   menuPlotPSTH_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function for the PSTH menu option under Plot.
%
% Inputs:
%    hObject   - Handle. The handle to the menuPlotPSTH (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
    
nMosaic = get(handles.listMosaics, 'Value');
rgcMosaic = handles.rgcLayer.mosaic{nMosaic}; 

% Try to force this into a new window with a flag.
rgcMosaic.plot('psth');

end

% --------------------------------------------------------------------
function menuLinearPreSpike_Callback(hObject, eventdata, handles)
% (Plot | Linear) Call the linear pre-spike menu option under Plot
%
% Syntax:
%   menuLinearPreSpike_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function for the Linear (pre-spike) option under Plot.
%
% Inputs:
%    hObject   - Handle. The handle to menuLinearPreSpike (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Plots the psth of all the cells combined. Kind of weird.
responseLinear = handles.rgcMosaic.get('responseLinear');

vcNewGraphWin;
plot(RGB2XWFormat(responseLinear)');       
end

% File Menu
% --------------------------------------------------------------------
function menuFileSave_Callback(hObject, eventdata, handles)
% (File | Save Image) Call Save Image under File
%
% Syntax:
%   menuFileSave_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function for Save Image under File.
%
% Inputs:
%    hObject   - Handle. The handle to menuFileSave (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
disp('Save NYI')
end

% --------------------------------------------------------------------
function menuFileRefresh_Callback(hObject, eventdata, handles)
% (File | Refresh) Call Refresh under File
%
% Syntax:
%   menuFileRefresh_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function for the Refresh menu object under File.
%
% Inputs:
%    hObject   - Handle. The handle to menuFileRefresh (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

rgcLayerWindowRefresh(handles)
end

% --------------------------------------------------------------------
function menuFileClose_Callback(hObject, eventdata, handles)
% (File | Close) Call Close under File
%
% Syntax:
%   menuFileClose_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function for the Close menu object under file.
%
% Inputs:
%    hObject   - Handle. The handle to the menuFileClose (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

handles.rgcM.fig = [];
delete(handles.mosaicWindow);
end

% --- Executes on selection change in listMosaics.
function listMosaics_Callback(hObject, eventdata, handles)
% The object containing the list of mosaics in the figure
%
% Syntax:
%   listMosaics_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function for the list object of the mosaics in the figure.
%
% Inputs:
%    hObject   - Handle. The handle to the listMosaics (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hints: To return the listMosaics contents as a cell array, use:
%           contents = cellstr(get(hObject, 'String'))
%      To return a selected item from listMosaics, use:
%           contents{get(hObject, 'Value')}
%

rgcLayerWindowRefresh(handles)
end

% --- Executes during object creation, after setting all properties.
function listMosaics_CreateFcn(hObject, eventdata, handles)
% The create function call for the listMosaic object
%
% Syntax:
%   listMosaics_CreateFcn(hObjet, eventdata, handles)
%
% Description:
%   The create function call for the listMosaic object
%
% Inputs:
%    hObject   - Handle. The handle to listMosaics (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty - handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: listbox controls usually have a white background on Windows.
%      See ISPC and COMPUTER.
%

if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end
end

function editGamma_Callback(hObject, eventdata, handles)
% editGamma box. Sets display gamma.
%
% Syntax:
%   editGamma_Callback(hObject, eventdata, handles)
%
% Description:
%    The call function to the editGamma object.
%
% Inputs:
%    hObject   - Handle. The handle to the editGamma (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Refresh to update for the new gamma value
rgcLayerWindowRefresh(handles)
end

% --- Executes during object creation, after setting all properties.
function editGamma_CreateFcn(hObject, eventdata, handles)
% The create function call for the editGamma object
%
% Syntax:
%   editGamma_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    The create function call for the editGamma object.
%
% Inputs:
%    hObject   - Handle. The handle to the menuPlotPSTH (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty - handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: edit controls usually have a white background on Windows.
%      See ISPC and COMPUTER.
%

if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end
end

%% Internal functions
function rgcLayerWindowRefresh(handles)
% Update all the text fields and such with the data in the mosaic
%
% Syntax:
%   rgcLayerWindowRefresh(handles)
%
% Description:
%    This function is used to update all of the text fields (etc) with the
%    data from the mosaic.
%
% Inputs:
%    handles - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

rgcL = handles.rgcLayer;
fig = figure(rgcL.fig);
gdata = guidata(fig);

% Show the appropriate response axis plot
axes(gdata.axisResponse);
cla(gdata.axisResponse, 'reset');

% Selected string in the popup
contents = cellstr(get(gdata.popupResponseSelect, 'String'));
str = contents{get(gdata.popupResponseSelect, 'Value')};

nMosaic = get(gdata.listMosaics, 'Value');
rgcL = gdata.rgcLayer;
rgcL.mosaic{nMosaic}.fig = rgcL.fig;
gam = str2double(get(gdata.editGamma', 'String'));

switch(str)
    case 'Receptive field mosaic'
        ieInWindowMessage('Building mosaic', handles);
        rgcL.mosaic{nMosaic}.plot('mosaic fill', 'gam', gam);
        ieInWindowMessage('', handles);

    case 'Spike mean (image)'
        rgcL.mosaic{nMosaic}.plot('spike mean image', 'gam', gam);
    case 'PSTH mean (image)'
        rgcL.mosaic{nMosaic}.plot('psth mean image', 'gam', gam);
    case 'Linear movie'
        ieInWindowMessage('Showing movie', handles, []);
        rgcL.mosaic{nMosaic}.plot('linear movie', 'gam', gam); 
        ieInWindowMessage('', handles, []);
    case 'Spike movie'
        ieInWindowMessage('Showing movie', handles, []);
        rgcL.mosaic{nMosaic}.plot('spike movie', 'gam', gam);
        ieInWindowMessage('', handles, []);
    case 'PSTH movie'
        ieInWindowMessage('Showing movie', handles, []);
        rgcL.mosaic{nMosaic}.plot('spike movie', 'gam', gam);
        ieInWindowMessage('', handles, []);
    otherwise
        error('Unknown plot type %s\n', str);
end

% Make a button for rfOverlay. ALways false, for now.
% rfOverlay = false;
% if rfOverlay, rgcL.mosaic{nMosaic}.plot('mosaic'); end

% Text description - implemented in rgcMosaic base class.
set(gdata.rgcProperties, 'string', rgcL.describe);

end
