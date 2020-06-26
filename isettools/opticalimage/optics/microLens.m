function varargout = microLens(varargin)
% M-file for microLens.fig
%
% Syntax:
%   varargout = microLens([varargin])
%
% Description:
%    Calling microlens, by itself, creates a new MICROLENS or raises the
%    existing singleton*.
%
%    H = MICROLENS returns the handle to a new MICROLENS or the handle to
%    the existing singleton*.
%
%    MICROLENS('CALLBACK', hObject, eventData, handles, ...) will call the
%    local function with the name CALLBACK in MICROLENS.M with the given
%    input arguments.
%
%    MICROLENS('Property', 'Value', ...) creates a new MICROLENS or raises
%    the existing singleton*. Starting from the left, property value pairs
%    are applied to the GUI before microLens_OpeningFunction gets called.
%    An unrecognized property name or invalid value makes property
%    application stop. All inputs are passed to microLens_OpeningFcn via
%    varargin.
%
%    *See GUI Options on GUIDE's Tools menu. Choose "GUI allows only one
%    instance to run (singleton)".
%
% Inputs:
%    None required.
%
% Outputs:
%    None required.
%
% Optional key/value pairs:
%    **NEEDS TO BE FILLED IN**
%
% Notes:
%    * Edit the above text to modify the response to help microLens
%
% See Also:
%    GUIDE, GUIDATA, GUIHANDLES
%

% History:
%    03/22/04       Last Modified by GUIDE v2.5 22-Mar-2004 17:59:57
%    03/08/18  jnm  Formatting

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton', gui_Singleton, ...
    'gui_OpeningFcn', @microLens_OpeningFcn, ...
    'gui_OutputFcn', @microLens_OutputFcn, ...
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

% --- Executes just before microLens is made visible.
function microLens_OpeningFcn(hObject, eventdata, handles, varargin)
% (microLens | Open) Opening function for microLens executes before visible
%
% Syntax:
%   microLens_OpeningFcn(hObject, eventdata, handles, [varargin])
%
% Description:
%    Function call to open the microLens object. Executes before the
%    microLens object is visible.
%
% Inputs:
%    hObject   - Handle. Handle to microLens
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%    varargin  - (Optional) VARIES. Optional Command Line arguments to
%                microLens. See VARARGIN.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Choose default command line output for microLens
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes microLens wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = microLens_OutputFcn(hObject, eventdata, handles)
% (microLens | Output) Command Line Output for the microLens
%
% Syntax:
%   varargout = microLens_OutputFcn(hObject, eventdata, handles)
%
% Description:
%    Function to output from the microLens to the Command Line.
%
% Inputs:
%    hObject   - Handle. Handle to microLens
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    varargout - (Optional) Cell Array for returning optional output
%                arguments. See VARARGOUT.
%
% Optional key/value pairs:
%    None.
%

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in btnCompute.
function btnCompute_Callback(hObject, eventdata, handles)
% (btnCompute) Object call to the Compute Button
%
% Syntax:
%   btnCompute_Callback(hObject, eventdata, handles)
%
% Description:
%    Compute Button call
%
% Inputs:
%    hObject   - handle to btnCompute (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% (edit1 | Create) Object call to Create edit1
%
% Syntax:
%   edit1_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Create function for edit1
%
% Inputs:
%    hObject   - Handle. Handle to edit1 (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - empty - handles are not created until after all of the
%                CreateFcns have been called
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: edit controls usually have a white background on Windows.
%    * See ISPC and COMPUTER.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
end

function edit1_Callback(hObject, eventdata, handles)
% (edit1) Object call to edit1
%
% Syntax:
%   edit1_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit1 call
%
% Inputs:
%    hObject   - Handle. Handle to edit1 (see GCBO)
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
%       get(hObject, 'String')
%           returns contents of edit1 as text
%       str2double(get(hObject, 'String'))
%           returns contents of edit1 as double
%
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% (edit2 | Create) Object call to Create edit2
%
% Syntax:
%   edit2_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Create function for edit2
%
% Inputs:
%    hObject   - Handle. Handle to edit2 (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - empty - handles are not created until after all of the
%                CreateFcns have been called
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: edit controls usually have a white background on Windows.
%    * See ISPC and COMPUTER.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
end

function edit2_Callback(hObject, eventdata, handles)
% (edit2) Object call to edit2
%
% Syntax:
%   edit2_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit2 call
%
% Inputs:
%    hObject   - Handle. Handle to edit2 (see GCBO)
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
%       get(hObject, 'String')
%           returns contents of edit2 as text
%       str2double(get(hObject, 'String'))
%           returns contents of edit2 as double
%
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% (edit3 | Create) Object call to Create edit3
%
% Syntax:
%   edit3_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Create function for edit3
%
% Inputs:
%    hObject   - Handle. Handle to edit3 (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - empty - handles are not created until after all of the
%                CreateFcns have been called
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: edit controls usually have a white background on Windows.
%    * See ISPC and COMPUTER.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
end

function edit3_Callback(hObject, eventdata, handles)
% (edit3) Object call to edit3
%
% Syntax:
%   edit3_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit3 call
%
% Inputs:
%    hObject   - Handle. Handle to edit3 (see GCBO)
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
%       get(hObject, 'String')
%           returns contents of edit3 as text
%       str2double(get(hObject, 'String'))
%           returns contents of edit3 as double
%
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% (edit4 | Create) Object call to Create edit4
%
% Syntax:
%   edit4_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Create function for edit4
%
% Inputs:
%    hObject   - Handle. Handle to edit4 (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - empty - handles are not created until after all of the
%                CreateFcns have been called
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: edit controls usually have a white background on Windows.
%    * See ISPC and COMPUTER.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
end

function edit4_Callback(hObject, eventdata, handles)
% (edit4) Object call to edit4
%
% Syntax:
%   edit4_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit4 call
%
% Inputs:
%    hObject   - Handle. Handle to edit4 (see GCBO)
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
%       get(hObject, 'String')
%           returns contents of edit4 as text
%       str2double(get(hObject, 'String'))
%           returns contents of edit4 as double
%
end

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% (edit5 | Create) Object call to Create edit5
%
% Syntax:
%   edit5_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Create function for edit5
%
% Inputs:
%    hObject   - Handle. Handle to edit5 (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - empty - handles are not created until after all of the
%                CreateFcns have been called
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: edit controls usually have a white background on Windows.
%    * See ISPC and COMPUTER.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
end

function edit5_Callback(hObject, eventdata, handles)
% (edit5) Object call to edit5
%
% Syntax:
%   edit5_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit5 call
%
% Inputs:
%    hObject   - Handle. Handle to edit5 (see GCBO)
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
%       get(hObject, 'String')
%           returns contents of edit5 as text
%       str2double(get(hObject, 'String'))
%           returns contents of edit5 as double
%
end

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% (edit6 | Create) Object call to Create edit6
%
% Syntax:
%   edit6_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Create function for edit6
%
% Inputs:
%    hObject   - Handle. Handle to edit6 (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - empty - handles are not created until after all of the
%                CreateFcns have been called
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: edit controls usually have a white background on Windows.
%    * See ISPC and COMPUTER.
%

if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
end

function edit6_Callback(hObject, eventdata, handles)
% (edit6) Object call to edit6
%
% Syntax:
%   edit6_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit6 call
%
% Inputs:
%    hObject   - Handle. Handle to edit6 (see GCBO)
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
%       get(hObject, 'String')
%           returns contents of edit6 as text
%       str2double(get(hObject, 'String'))
%           returns contents of edit6 as double
%
end

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% (edit7 | Create) Object call to Create edit7
%
% Syntax:
%   edit7_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Create function for edit7
%
% Inputs:
%    hObject   - Handle. Handle to edit7 (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - empty - handles are not created until after all of the
%                CreateFcns have been called
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: edit controls usually have a white background on Windows.
%    * See ISPC and COMPUTER.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
end

function edit7_Callback(hObject, eventdata, handles)
% (edit7) Object call to edit7
%
% Syntax:
%   edit7_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit7 call
%
% Inputs:
%    hObject   - Handle. Handle to edit7 (see GCBO)
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
%       get(hObject, 'String')
%           returns contents of edit7 as text
%       str2double(get(hObject, 'String'))
%           returns contents of edit7 as double
%
end

% --------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)
% (Menu | File) Object call to the File menu option
%
% Syntax:
%   menuFile_Callback(hObject, eventdata, handles)
%
% Description:
%    Top level of File menu
%
% Inputs:
%    hObject   - Handle. Handle to menuFile (see GCBO)
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
% (Menu | Edit) Object call to the Edit menu option
%
% Syntax:
%   menuEdit_Callback(hObject, eventdata, handles)
%
% Description:
%    Top level of Edit menu
%
% Inputs:
%    hObject   - Handle. Handle to menuEdit (see GCBO)
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
function menuPlot_Callback(hObject, eventdata, handles)
% (Menu | Plot) Object call to the Plot menu option
%
% Syntax:
%   menuPlot_Callback(hObject, eventdata, handles)
%
% Description:
%    Top level of Plot menu
%
% Inputs:
%    hObject   - Handle. Handle to menuPlot (see GCBO)
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
% (Menu | Analyze) Object call to the Analyze menu option
%
% Syntax:
%   menuAnalyze_Callback(hObject, eventdata, handles)
%
% Description:
%    Top level of Analyze menu
%
% Inputs:
%    hObject   - Handle. Handle to menuAnalyze (see GCBO)
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
function menuHelp_Callback(hObject, eventdata, handles)
% (Menu | Help) Object call to the Help menu option
%
% Syntax:
%   menuHelp_Callback(hObject, eventdata, handles)
%
% Description:
%    Top level of Help menu
%
% Inputs:
%    hObject   - Handle. Handle to menuHelp (see GCBO)
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