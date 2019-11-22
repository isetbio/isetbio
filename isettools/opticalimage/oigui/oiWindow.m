function varargout = oiWindow(varargin)
% Optical image window
%
% Syntax:
%   [varargout] = oiWindow([varargin])
%
% Description:
%    Graphical user interface to manage the ISET OPTICALIMAGE properties.
%
%    OIWINDOW, by itself, creates a new OIWINDOW or raises the existing
%    singleton*.
%
%    H = OIWINDOW returns the handle to a new OIWINDOW or the handle to
%    the existing singleton*.
%
%    H = OIWINDOW(oi) adds the oi to the database and opens the window.
%    Equivalent to ieAddObject(oi); oiWindow;
%
%    OIWINDOW('CALLBACK', hObject, eventData, handles, ...) calls the local
%    function named CALLBACK in OIWINDOW.M with the given input arguments.
%
%    OIWINDOW('Property', 'Value', ...) creates a new OIWINDOW or raises
%    the existing singleton*. Starting from the left, property value pairs
%    are applied to the GUI before oiWindow_OpeningFunction gets called. An
%    unrecognized property name or invalid value makes property application
%    stop. All inputs are passed to oiWindow_OpeningFcn via varargin.
%
% Inputs:
%    None required.
%
% Outputs:
%    None required.
%
% Optional key/value pairs:
%    **Needs to be filled out**
%

% History:
%    xx/xx/03         Copyright ImagEval Consultants, LLC, 2003.
%    01/25/17  GUIDE  Last Modified by GUIDE v2.5 25-Jan-2017 15:25:04
%    03/19/18  jnm    Formatting (Fin 03/22)
%    06/18/19  JNM    Formatting update (also 07/01/19)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton', gui_Singleton, ...
    'gui_OpeningFcn', @oiWindow_OpeningFcn, ...
    'gui_OutputFcn', @oiWindow_OutputFcn, ...
    'gui_LayoutFcn', [] , ...
    'gui_Callback', []);
if nargin && isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before oiWindow is made visible.
function oiWindow_OpeningFcn(hObject, eventdata, handles, varargin)
% Set defaults for oiWindow
%
% Syntax:
%   oiWindow_OpeningFcn(hObject, eventdata, handles)
%
% Description:
%    Set up defaults for oiWindow just before it is made visible.
%
% Inputs:
%    hObject   - Handle. Handle to oiWindow (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Choose default command line output for microLensWindow
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
vcSetFigureHandles('OI', hObject, eventdata, handles);

%  Check the preferences for ISET and adjust the font size.
ieFontInit(hObject);

if ~isempty(varargin)
    oi = varargin{1};
    if strcmp(oi.type,'opticalimage'), ieAddObject(oi);
    else, warning('Unexpected variable input.\n');
    end
end
oiRefresh(hObject, eventdata, handles);

return

% --- Outputs from this function are returned to the command line.
function varargout = oiWindow_OutputFcn(hObject, eventdata, handles)
% Default command line output from handles structure
%
% Syntax:
%   varargout = oiWindow_OutputFcn(hObject, eventdata, handles)
%
% Description:
%    Function to retrieve the default command line output contained within
%    the handles structure.
%
% Inputs:
%    hObject   - Handle. Handle to oiWindow (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
varargout{1} = handles.output;
return;

% --- Executes on button press in btnDeleteOptImg.
function oiDelete(hObject, eventdata, handles)
% (Edit | Delete Current OI) Function called by menu selection
%
% Syntax:
%   oiDelete(hObject, eventdata, handles)
%
% Description:
%    Delete current OI as called from within menu selection function.
%
% Inputs:
%    hObject   - Handle. Handle to oiDelete (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
vcDeleteSelectedObject('OPTICALIMAGE');
[val, oi] = vcGetSelectedObject('OPTICALIMAGE');
if isempty(oi)
    oi = oiCreate;
    ieReplaceObject(oi, 1);
end

oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuEditDeleteSome_Callback(hObject, eventdata, handles)
% (Edit | Delete Some OIs) Option to delete SOME OI's (not singular)
%
% Syntax:
%   menuEditDeleteSome_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option under Edit to delete some OI's (not current oi, that is a
%    separate option).
%
% Inputs:
%    hObject   - Handle. Handle to menuEditDeleteSome (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
vcDeleteSomeObjects('oi');
oiRefresh(hObject, eventdata, handles);
return;

% --- Executes during object creation, after setting all properties.
function SelectOptImg_CreateFcn(hObject, eventdata, handles)
% (Create | OI) Create ability to select an OI according to defaults
%
% Syntax:
%   SelectOptImg_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Creates the option to select an optical image according to the
%    pre-defined default values.
%
% Inputs:
%    hObject   - Handle. Handle to SelectOptImg (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty - The handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in SelectOptImg.
function SelectOptImg_Callback(hObject, eventdata, handles)
% (Select | OI) Select an Optical Image
%
% Syntax:
%   SelectOptImg_Callback(hObject, eventdata, handles)
%
% Description:
%    Optical image selection change callback
%
% Inputs:
%    hObject   - Handle. Handle to SelectOptImg (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiNames = get(hObject, 'String');
thisName = oiNames{get(hObject, 'Value')};

switch lower(thisName)
    case 'new'
        oiNew(hObject, eventdata, handles);
    otherwise
        val = get(hObject, 'Value')-1;
        vcSetSelectedObject('OPTICALIMAGE', val);
        [val, oi] = vcGetSelectedObject('OPTICALIMAGE');
end

oiRefresh(hObject, eventdata, handles);
return;

function oiRefresh(hObject, eventdata, handles)
% (Refresh) Refresh function called by Menu selection
%
% Syntax:
%   oiRefresh(hObject, eventdata, handles)
%
% Description:
%    Refresh function called when menu option selected
%
% Inputs:
%    hObject   - Handle. Handle to oiRefresh (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiSetEditsAndButtons(handles);
return;

%%%%%%%%%%%%%%%%%%%% Menus are controlled below here %%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% (File) Top File Menu object
%
% Syntax:
%   FileMenu_Callback(hObject, eventdata, handles)
%
% Description:
%    The top level of the File Menu option
%
% Inputs:
%    hObject   - Handle. Handle to FileMenu (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --------------------------------------------------------------------
function menuSaveImage_Callback(hObject, eventdata, handles)
% (File | Save Image) Save OI as an image.
%
% Syntax:
%   menuSaveImage_Callback(hObject, eventdata, handles)
%
% Description:
%    Function to save the current OI as an image (.png).
%
% Inputs:
%    hObject   - Handle. Handle to menuSaveImage (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
gam = str2double(get(handles.editGamma, 'String'));
oi = vcGetObject('OPTICALIMAGE');
oiSaveImage(oi, [], gam);
return;

% --------------------------------------------------------------------
function menuFileClose_Callback(hObject, eventdata, handles)
% (File | Close) Close option under File
%
% Syntax:
%   menuFileClose_Callback(hObject, eventdata, handles)
%
% Description:
%    Calls the close function under the File menu heading.
%
% Inputs:
%    hObject   - Handle. Handle to menuFileClose (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiClose;
return;

% --------------------------------------------------------------------
function EditMenu_Callback(hObject, eventdata, handles)
% (Edit) Top level menu option for Edit
%
% Syntax:
%   EditMenu_Callback(hObject, eventdata, handles)
%
% Description:
%    Top level menu option for Edit
%
% Inputs:
%    hObject   - Handle. Handle to EditMenu (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --- Executes during object creation, after setting all properties.Ex
function editFnumber_CreateFcn(hObject, eventdata, handles)
% (Create | F-Number) Creates the F-number textbox following defaults
%
% Syntax:
%   editFnumber_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Creates the F-number textbox according to the pre-existing defaults.
%
% Inputs:
%    hObject   - Handle. Handle to editFnumber (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty - The handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

%---------------------------------
function editFnumber_Callback(hObject, eventdata, handles)
% (Edit | F-Number) F-number text edit call back (fnumber)
%
% Syntax:
%   editFnumber_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit the value located in the F-number text box.
%
% Inputs:
%    hObject   - Handle. Handle to editFnumber (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% The f-number is the ratio of the focal length divided by the aperture.
[val, oi] = vcGetSelectedObject('OPTICALIMAGE');
fNumber = str2double(get(hObject, 'String'));
oi = oiSet(oi, 'optics fnumber', fNumber);

vcReplaceObject(oi, val);
oiRefresh(hObject, eventdata, handles);
return;

% --- Executes during object creation, after setting all properties.
function editFocalLength_CreateFcn(hObject, eventdata, handles)
% (Create | Focal Length) Create Focal Length Textbox according to defaults
%
% Syntax:
%   editFocalLength_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Function to create the Focal Length textbox, according to the existing
%    default values.
%
% Inputs:
%    hObject   - Handle. Handle to editFocalLength (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

function editFocalLength_Callback(hObject, eventdata, handles)
% (Edit | Focal Lenght) Focal length textbox edit
%
% Syntax:
%   editFocalLength_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit the value in the Focal Length textbox.
%
% Inputs:
%    hObject   - Handle. Handle to editFocalLength (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Read the edit box
focalLength = str2double(get(hObject, 'String')) / 1000;

% Set the optics focal length
[val, oi] = vcGetSelectedObject('OPTICALIMAGE');
oi = oiSet(oi, 'optics focal length', focalLength);

vcReplaceObject(oi, val);
oiRefresh(hObject, eventdata, handles);
return;

% --- Executes during object creation, after setting all properties.
function editDefocus_CreateFcn(hObject, eventdata, handles)
% (Create | Defocus) Defocus creation function
%
% Syntax:
%   editDefocus_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Defocus creation function
%
% Inputs:
%    hObject   - Handle. Handle to editDefocus (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty - Handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

% --- Executes on button press in btnSimulate.
function btnSimulate_Callback(hObject, eventdata, handles)
% (Button | Press) Simulate triggered by Button Press
%
% Syntax:
%   btnSimulate_Callback(hObject, eventdata, handles)
%
% Description:
%    This function simulates a button press.
%
%    This callback reads the current scene and optics and then calculates
%    the optical image irradiance with the current parameters. We do not
%    calculate a new optical image. [Note: XXX - Probably, we should put
%    this in a separate function rather than keeping it in here.]
%
% Inputs:
%    hObject   - Handle. Handle to menuEditZoom (see GCBO)
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
%    * TODO: Determine if note in description section has merit (should we
%      move/remove the described functionality?
%
scene = vcGetObject('scene');
if isempty(scene)
    ieInWindowMessage('No scene data.', handles);
    beep;
    return;
else
    ieInWindowMessage('', handles);
end

oi = vcGetObject('OPTICALIMAGE');

% Check oiCompute whether the custom button is selected or not.
oi = oiCompute(scene, oi);
oi = oiSet(oi, 'consistency', 1);

% Save the OI in the vcSESSION as the selected optical image.
ieReplaceObject(oi);

% hObject = oiwindow;
oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuFileLoadOI_Callback(hObject, eventdata, handles)
% (File | Load OI Data) Load an OI into the current session
%
% Syntax:
%   menuFileLoadOI_Callback(hObject, eventdata, handles)
%
% Description:
%    Load OI function. Function loads OI data from a .mat file into the
%    current session.
%
% Inputs:
%    hObject   - Handle. Handle to menuFileLoadOI (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
newVal = vcImportObject('OPTICALIMAGE');
vcSetSelectedObject('OPTICALIMAGE', newVal);
oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuFileSaveOI_Callback(hObject, eventdata, handles)
% (File | Save OI Data) Save OI (.mat) function
%
% Syntax:
%   menuFileSaveOI_Callback(hObject, eventdata, handles)
%
% Description:
%    Save OI function. Function saves the current OI to a .mat file.
%
% Inputs:
%    hObject   - Handle. Handle to menuFileSaveOI (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, oi] = vcGetSelectedObject('OPTICALIMAGE');
fullName = vcSaveObject(oi);
return;

% --- Executes during object creation, after setting all properties.
function editGamma_CreateFcn(hObject, eventdata, handles)
% (Edit | Gamma Value) Create the Gamma Text Box according to defaults
%
% Syntax:
%   editGamma_Callback(hObject, eventdata, handles)
%
% Description:
%    Create the Gamma textbox and value according to defaults.
%
% Inputs:
%    hObject   - Handle. Handle to editGamma (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty - The handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

function editGamma_Callback(hObject, eventdata, handles)
% (Edit | Gamma) Gamma textbox alterations
%
% Syntax:
%   editGamma_Callback(hObject, eventdata, handles)
%
% Description:
%    Gamma textbox alterations.
%
%    When we refresh the GUI the value is read and the image is displayed
%    with the new gamma value.
%
% Inputs:
%    hObject   - Handle. Handle to editGamma (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiRefresh(hObject, eventdata, handles);
return;

% --- Executes on selection change in popupDisplay.
function popupDisplay_Callback(hObject, eventdata, handles)
% (PopUp | Display) Refresh.
%
% Syntax:
%   popupDisplay_Callback(hObject, eventdata, handles)
%
% Description:
%    When we refresh, the rendering method is read and the oiShowImage will
%    call the relevant rendering routine.
%
% Inputs:
%    hObject   - Handle. Handle to popupDisplay (see GCBO)
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
%      contents = get(hObject, 'String') - Will return popupDisplay
%                                          contents as a cell array
%      contents{get(hObject, 'Value')}   - Will return the selected item
%                                          from popupDisplay
oiRefresh(hObject, eventdata, handles);
return

% --- Executes during object creation, after setting all properties.
function popupDisplay_CreateFcn(hObject, eventdata, handles)
% (PopUp | Display | Create) Create the popUp following defaults.
%
% Syntax:
%   popupDisplay_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Create the popUp display.
%
% Inputs:
%    hObject   - Handle. Handle to popupDisplay (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty - The handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: popupmenu controls usually have a white background on Windows.
%      See ISPC and COMPUTER.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return

% --- Executes on button press in btnNew.
function oiNew(hObject, eventdata, handles)
% (Edit | Create new OI) Create a new OI
%
% Syntax:
%   oiNew(hObject, eventdata, handles)
%
% Description:
%    Create New OI under the Edit menu.
%
% Inputs:
%    hObject   - Handle. Handle to oiNew (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, newOI] = vcGetSelectedObject('opticalimage');
newVal = vcNewObjectValue('opticalimage');
newOI.name = vcNewObjectName('opticalimage');
newOI.type = 'opticalimage';
newOI = oiClearData(newOI);
vcAddAndSelectObject('OPTICALIMAGE', newOI);
oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuLens_Callback(hObject, eventdata, handles)
% (Optics | Lens) Lens sub-heading
%
% Syntax:
%   menuLens_Callback(hObject, eventdata, handles)
%
% Description:
%    Lens sub-heading under Optics
%
% Inputs:
%    hObject   - Handle. Handle to menuLens (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --------------------------------------------------------------------
function menuLensDensity_Callback(hObject, eventdata, handles)
% (Optics | Lens | Lens Density) Set the human lens density
%
% Syntax:
%   menuLensDensity_Callback(hObject, eventdata, handles)
%
% Description:
%    Set the human lens density.
%
% Inputs:
%    hObject   - Handle. Handle to menuLensDensity (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('oi');
d = oiGet(oi, 'lens density');
val = ieReadNumber('Enter human lens pigment density', d, '%.2f');
if isempty(val), disp('Canceled'); return; end
oi = oiSet(oi, 'lens density', val);
vcReplaceObject(oi);
return;

% --------------------------------------------------------------------
function menuOptTrans_Callback(hObject, eventdata, handles)
% (Optics | Lens | Transmittance)
%
% Syntax:
%   menuOptTrans_Callback(hObject, eventdata, handles)
%
% Description:
%    Read the optical transmittance in wavelength. This is useful for
%    diffraction cases.
%
%    We could use a function that multiplies the transmittance by another
%    function, such as a lens or macular pigment transmittance. As things
%    stand we load a transmittance, but we should probably have a function
%    that gets the existing one and multipllies it by another.
%
% Inputs:
%    hObject   - Handle. Handle to menuOptTrans (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, oi] = vcGetSelectedObject('OI');
optics = oiGet(oi, 'optics');
wave = opticsGet(oi, 'wave');
fullName = vcSelectDataFile('optics');
if isempty(fullName)
    return;
else
    optics = opticsSet(optics, 'transmittance', ...
        ieReadSpectra(fullName, wave));
end
oi = oiSet(oi, 'optics', optics);
vcReplaceObject(oi, val);
return;

% --------------------------------------------------------------------
function menuFileRefresh_Callback(hObject, eventdata, handles)
% (File | Refresh) Refresh function
%
% Syntax:
%   menuFileRefresh_Callback(hObject, eventdata, handles)
%
% Description:
%    Refresh function under File.
%
% Inputs:
%    hObject   - Handle. Handle to menuFileRefresh (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuEditScale_Callback(hObject, eventdata, handles)
% (Edit | Scale irradiance) Scale irradiance levels by s
%
% Syntax:
%   menuEditScale_Callback(hObject, eventdata, handles)
%
% Description:
%    Scale the irradiance levels by s - received from popUp.
%
% Inputs:
%    hObject   - Handle. Handle to menuEditScale (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
s = ieReadNumber('Enter scale factor', 1, ' %.2f');
if isempty(s), return; end

[val, oi] = vcGetSelectedObject('OPTICALIMAGE');
irrad = oiGet(oi, 'photons');
if isempty(irrad)
    handles = ieSessionGet('opticalimagehandle');
    ieInWindowMessage('Can not scale:  No irradiance data.', handles, []);
else
    handles = ieSessionGet('opticalimagehandle');
    ieInWindowMessage('', handles, []);
end
ill = oiGet(oi, 'illuminance');
meanIll = oiGet(oi, 'meanIlluminance');
oi = oiSet(oi, 'photons', irrad*s);
if ~isempty(ill), oi = oiSet(oi, 'illuminance', s * ill); end
if ~isempty(meanIll), oi = oiSet(oi, 'meanIlluminance', s * meanIll); end
ieReplaceObject(oi, val)
oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuEditFontSize_Callback(hObject, eventdata, handles)
% (Edit | Change Font Size) Font Size Alterations
%
% Syntax:
%   menuEditFontSize_Callback(hObject, eventdata, handles)
%
% Description:
%    Function to change the font size.
%
% Inputs:
%    hObject   - Handle. Handle to menuEditFontSize (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ieFontSizeSet(handles.figure1);
return;

% --------------------------------------------------------------------
function menuEditName_Callback(hObject, eventdata, handles)
% (Edit | Re-name OI) Re-naming function
%
% Syntax:
%   menuEditName_Callback(hObject, eventdata, handles)
%
% Description:
%    OI Re-naming function
%
% Inputs:
%    hObject   - Handle. Handle to menuEditName (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

[val, oi] = vcGetSelectedObject('OPTICALIMAGE');
newName = ieReadString('New optical image name', 'new-oi');

if isempty(newName), return;
else, oi = oiSet(oi, 'name', newName);
end

ieReplaceObject(oi, val)
oiRefresh(hObject, eventdata, handles);

return;
% --------------------------------------------------------------------
function menuCopyOI_Callback(hObject, eventdata, handles)
% (Edit | Copy OI) Edit's Copy OI function
%
% Syntax:
%   menuCopyOI_Callback(hObject, eventdata, handles)
%
% Description:
%    Copy OI function contained under Edit.
%
% Inputs:
%    hObject   - Handle. Handle to menuCopyOI (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, oi] = vcGetSelectedObject('OI');
newName = ieReadString('New optical image name', 'new-oi');
if isempty(newName), return; else, oi = oiSet(oi, 'name', newName); end
vcAddAndSelectObject('OPTICALIMAGE', oi);
oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuEditDelete_Callback(hObject, eventdata, handles)
% (Edit | Delete Current OI) Delete OI Function (current, not some)
%
% Syntax:
%   menuEditDelete_Callback(hObject, eventdata, handles)
%
% Description:
%    Delete current OI function (not delete some)
%
% Inputs:
%    hObject   - Handle. Handle to menuEditDelete (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiDelete(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuEditCreate_Callback(hObject, eventdata, handles)
% (Edit | Create new OI) Create New OI
%
% Syntax:
%   menuEditCreate_Callback(hObject, eventdata, handles)
%
% Description:
%    Create New OI.
%
% Inputs:
%    hObject   - Handle. Handle to menuEditCreate (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiNew(hObject, eventdata, handles);
return;

function menuEditClearMessage_Callback(hObject, eventdata, handles)
% (Edit | Clear Window Message) Clear Window Message Function
%
% Syntax:
%   menuEditClearMessage_Callback(hObject, eventdata, handles)
%
% Description:
%    Clear Window Message function.
%
% Inputs:
%    hObject   - Handle. Handle to menuEditClearMessage (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ieInWindowMessage('', ieSessionGet('opticalimagehandle'), []);
return;

% --------------------------------------------------------------------
function menuEditZoom_Callback(hObject, eventdata, handles)
% (Edit | Zoom) Zoom function
%
% Syntax:
%   menuEditZoom_Callback(hObject, eventdata, handles)
%
% Description:
%    Zoom function
%
% Inputs:
%    hObject   - Handle. Handle to menuEditZoom (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
zoom
return;

% --------------------------------------------------------------------
function menuEditViewer_Callback(hObject, eventdata, handles)
% (Edit | Viewer) Load the SI Data
%
% Syntax:
%   menuOpticsLoadSI_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot the irradiance in terms of photons.
%
% Inputs:
%    hObject   - Handle. Handle to menuEditViewer (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('oi');
img = oiGet(oi, 'photons');
rgb = imageSPD(img, oiGet(oi, 'wavelength'));
ieViewer(rgb);
return;

% --------------------------------------------------------------------
function menuOptics_Callback(hObject, eventdata, handles)
% (Optics) Top level Optics options
%
% Syntax:
%   menuOptics_Callback(hObject, eventdata, handles)
%
% Description:
%    Top level of Optics menu.
%
% Inputs:
%    hObject   - Handle. Handle to menuOptics (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --------------------------------------------------------------------
function menuOpticsHuman_Callback(hObject, eventdata, handles)
% (Optics | Human optics (MW)) The Human optics (mw) Optics option.
%
% Syntax:
%   menuOpticsHuman_Callback(hObject, eventdata, handles)
%
% Description:
%    Optics -> Human Optics (mw)
%
% Inputs:
%    hObject   - Handle. Handle to menuOpticsHuman (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = oiCreate('human');
ieAddObject(oi);
oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuHumanWVF_Callback(hObject, eventdata, handles)
% (Optics | Human (WVF)) Human (wvf) option under Optics
%
% Syntax:
%   menuHumanWVF_Callback(hObject, eventdata, handles)
%
% Description:
%    Human 'wvf' optics option.
%
% Inputs:
%    hObject   - Handle. Handle to menuHumanWVF (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = oiCreate('wvf human');
ieAddObject(oi);
oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function Diffraction_Callback(hObject, eventdata, handles)
% (Optics | Diffraction) Diffraction option under the optics menu
%
% Syntax:
%   Diffraction_Callback(hObject, eventdata, handles)
%
% Description:
%    The diffraction option under Optics.
%
% Inputs:
%    hObject   - Handle. Handle to Diffraction (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = oiCreate('diffraction');
ieAddObject(oi);
oiRefresh(hObject, eventdata, handles);
return;

% % --------------------------------------------------------------------
% function menuOpticsHalfInch_Callback(hObject, eventdata, handles)
%
% [val, oi] = vcGetSelectedObject('OPTICALIMAGE');
% oi = oiClearData(oi);
% optics = opticsCreate('standard (1/2-inch)');
% oi = oiSet(oi, 'optics', optics);
% vcReplaceObject(oi, val);
%
% oiRefresh(hObject, eventdata, handles);
%
% return;
%
% % --------------------------------------------------------------------
% function menuOpticsQuarterInch_Callback(hObject, eventdata, handles)
%
% [val, oi] = vcGetSelectedObject('OPTICALIMAGE');
% oi = oiClearData(oi);
% optics = opticsCreate('standard (1/4-inch)');
% oi = oiSet(oi, 'optics', optics);
% vcReplaceObject(oi, val);
% oiRefresh(hObject, eventdata, handles);
%
% return;
%
% % --------------------------------------------------------------------
% function menuOpticsThird_Callback(hObject, eventdata, handles)
% [val, oi] = vcGetSelectedObject('OPTICALIMAGE');
%
% optics = opticsCreate('standard (1/3-inch)');
% oi.optics = optics;
% oi.data = [];
% oi = sceneClearData(oi);
% vcReplaceObject(oi, val);
% oiRefresh(hObject, eventdata, handles);
% return;
%
% % --------------------------------------------------------------------
% function menuOpticsTwoThirds_Callback(hObject, eventdata, handles)
%
% [val, oi] = vcGetSelectedObject('OPTICALIMAGE');
% optics = opticsCreate('standard (2/3-inch)');
% oi.optics = optics;
% oi = sceneClearData(oi);
% vcReplaceObject(oi, val);
% oiRefresh(hObject, eventdata, handles);
% return;
%
% % --------------------------------------------------------------------
% function menuOpticsInch_Callback(hObject, eventdata, handles)
%
% [val, oi] = vcGetSelectedObject('OPTICALIMAGE');
%
% optics = opticsCreate('standard (1-inch)');
% oi.optics = optics;
% oi = sceneClearData(oi);
% vcReplaceObject(oi, val);
%
% oiRefresh(hObject, eventdata, handles);
% return;
% --------------------------------------------------------------------
% function menuHuman_Callback(hObject, eventdata, handles)
% return;

% % --------------------------------------------------------------------
% function menuMacular028_Callback(hObject, eventdata, handles)
% %
% [val, oi] = vcGetSelectedObject('OPTICALIMAGE');
% oi = humanMacularTransmittance(oi, 0.28);
% vcReplaceObject(oi, val);
%
% oiRefresh(hObject, eventdata, handles);
% return;
%
% % --------------------------------------------------------------------
% function menuMacular_Callback(hObject, eventdata, handles)
%
% [val, oi] = vcGetSelectedObject('OPTICALIMAGE');
%
% % Could use ieReadNumber here.
% dens = ieReadNumber('Enter macular density', 0.28, ' %.2f');
% % prompt = {'Enter macular density:'};
% % def = {'0.28'};
% % dlgTitle = 'Macular pigment density';
% % lineNo = 1; answer = inputdlg(prompt, dlgTitle, lineNo, def);
% % dens = str2num(answer{1});
% oi = humanMacularTransmittance([], dens);
% vcReplaceObject(oi, val);
%
% oiRefresh(hObject, eventdata, handles);
% return;


% --------------------------------------------------------------------
function menuOpticsImport_Callback(hObject, eventdata, handles)
% (Optics | Import Optics) Import Optics
%
% Syntax:
%   menuOpticsImport_Callback(hObject, eventdata, handles)
%
% Description:
%    Import Optics
%
% Inputs:
%    hObject   - Handle. Handle to menuOpticsImport (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
vcImportObject('OPTICS');
oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuOpticsRename_Callback(hObject, eventdata, handles)
% (Optics | Re-name Optics) Re-name the existing optics
%
% Syntax:
%   menuOpticsRename_Callback(hObject, eventdata, handles)
%
% Description:
%    Rename the existing optics
%
% Inputs:
%    hObject   - Handle. Handle to menuOpticsRename (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OI');
optics = oiGet(oi, 'optics');

if oiGet(oi, 'customCompute')
    name = ieReadString('Enter new ray trace optics name');
    if isempty(name), return; end
    optics = opticsSet(optics, 'rtname', name);
else
    name = ieReadString('Enter new diffraction limited optics name');
    if isempty(name), return; end
    optics = opticsSet(optics, 'name', name);
end
oi = oiSet(oi, 'optics', optics');

vcReplaceObject(oi);
oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuOpticsExports_Callback(hObject, eventdata, handles)
% (Optics | Export optics) Export the current optics data
%
% Syntax:
%   menuOpticsExports_Callback(hObject, eventdata, handles)
%
% Description:
%    Export the current Optics  to a file.
%
% Inputs:
%    hObject   - Handle. Handle to menuOpticsExport (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, optics] = vcGetSelectedObject('OPTICS');
vcExportObject(optics);
return;

% --------------------------------------------------------------------
function menuOpticsLoadSI_Callback(hObject, eventdata, handles)
% (Optics | Load SI data) Load the SI Data
%
% Syntax:
%   menuOpticsLoadSI_Callback(hObject, eventdata, handles)
%
% Description:
%    Load the SI Data.
%
% Inputs:
%    hObject   - Handle. Handle to menuOpticsLoadSI (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% The user selects a file containing the shift-invariant data.
[val, oi] = vcGetSelectedObject('OPTICALIMAGE');
optics = siSynthetic('custom', oi);
oi = oiSet(oi, 'optics', optics);
vcReplaceObject(oi, val);
return;

% --------------------------------------------------------------------
function PlotMenu_Callback(hObject, eventdata, handles)
% (Plot) Top level menu of Plot
%
% Syntax:
%   PlotMenu_Callback(hObject, eventdata, handles)
%
% Description:
%    This is the top level of the Plot Menu.
%
% Inputs:
%    hObject   - Handle. Handle to PlotMenu (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --------------------------------------------------------------------
function plotIrradiance_Callback(hObject, eventdata, handles)
% (Plot | Irradiance (photons)) Plot the irradiance in terms of photons.
%
% Syntax:
%   plotIrradiance_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot the irradiance in terms of photons.
%
% Inputs:
%    hObject   - Handle. Handle to plotIrradiance (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('oi');
oiPlot(oi, 'irradiance photons roi');
return;

% --------------------------------------------------------------------
function menuPlotIrradEnergy_Callback(hObject, eventdata, handles)
% (Plot | Irradiance (energy)) Plot the irradiance in terms of energy.
%
% Syntax:
%   menuPlotIrradEnergy_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot the irradiance in terms of energy.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotIrradEnergy (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('oi');
oiPlot(oi, 'irradiance energy roi');
return;

% --------------------------------------------------------------------
function menuPlotImageGrid_Callback(hObject, eventdata, handles)
% (Plot | Image | RGB | Grid) Plot the current RGB image in grid window.
%
% Syntax:
%   menuPlotImageGrid_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot the RGB image grid.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotImageGrid (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiPlot(vcGetObject('oi'), 'irradiance image grid');
return;

% --------------------------------------------------------------------
function menuPlotLens_Callback(hObject, eventdata, handles)
% (Plot | Lens transmittance)
%
% Syntax:
%   menuPlotLens_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot the lens transmittance.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotLens (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('oi');
oiPlot(oi, 'lens transmittance');
return;

% --------------------------------------------------------------------
function menuPlotDepthmap_Callback(hObject, eventdata, handles)
% (Plot | Depth Map) Plot the Depth Map
%
% Syntax:
%   menuPlotDepthmap_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot the Depth Map.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotDepthmap (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('oi');
if isempty(oiGet(oi, 'depth map'))
    handles = ieSessionGet('opticalimagehandle');
    ieInWindowMessage('No depth data.', handles, 3);
else
    oiPlot(oi, 'depth map');
end

return

% --------------------------------------------------------------------
function menuPlotDepthContour_Callback(hObject, eventdata, handles)
% (Plot | Depth Contour) Plot the Depth Contour
%
% Syntax:
%   menuPlotDepthContour_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot -> Depth Contour
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotVLContrast (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('oi');
if isempty(oiGet(oi, 'depth map'))
    handles = ieSessionGet('opticalimagehandle');
    ieInWindowMessage('No depth data.', handles, 3);
else
    oiPlot(oi, 'depth map contour');
end
return;

% --------------------------------------------------------------------
function menuPlotHLContrast_Callback(hObject, eventdata, handles)
% (Analyze | LxW Plots | Horizontal Line | Contrast) hLine Contrast
%
% Syntax:
%   menuPlotHLContrast_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze -> Line x wave plots -> Horizontal line (Contrast)
%
%    This may never be called. If it is, it is from the Analyze pulldown.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotHLContrast (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OPTICALIMAGE');
oiPlot(oi, 'hlinecontrast');
return;

% --------------------------------------------------------------------
function menuPlotVLContrast_Callback(hObject, eventdata, handles)
% (Analyze | LxW Plots | Vertical Line | Contrast) vLine Contrast
%
% Syntax:
%   menuPlotVLContrast_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze -> Line x wave plots -> Vertical line (Contrast)
%
%    This may never be called. If it is, it is from the Analyze pulldown.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotVLContrast (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OPTICALIMAGE');
oiPlot(oi, 'vlinecontrast');
return;

% --------------------------------------------------------------------
function menuPlotIllumLog_Callback(hObject, eventdata, handles)
% (Analyze | Illuminance Mesh (Log10)) Analyze the Log10 Illuminance Mesh
%
% Syntax:
%   menuPlotIllumLog_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze the log10 Illuminance Mesh option.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotIllumLog (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, oi] = vcGetSelectedObject('OI');
if ~checkfields(oi, 'data', 'illuminance')
    illuminance = oiCalculateIlluminance(oi);
    oi = oiSet(oi, 'illuminance', illuminance);
    vcReplaceObject(oi, val);
end
% Plots log10 or linear luminance,
% oiPlotIlluminance(oi, 'log');
oiPlot(oi, 'illuminance mesh log');
return;

% --------------------------------------------------------------------
function menuPlotIllumLin_Callback(hObject, eventdata, handles)
% (Analyze | Illuminance Mesh (Linear)) Analyze the Linear Illuminance Mesh
%
% Syntax:
%   menuPlotIllumLin_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze the Linear Illuminance Mesh option.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotIllumLin (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, oi] = vcGetSelectedObject('OPTICALIMAGE');

if ~checkfields(oi, 'data', 'illuminance')
    [oi.data.illuminance, oi.data.meanIll] = oiCalculateIlluminance(oi);
    vcReplaceObject(oi, val);
end
% Plots log10 or linear luminance,
oiPlot(oi, 'illuminance mesh linear');
return;

% --------------------------------------------------------------------
function menuPlotCIE_Callback(hObject, eventdata, handles)
% (Plot | ROI Summaries | Chromacity) Chromacity (xy) ROI Summary
%
% Syntax:
%   menuPlotCIE_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze -> Optics -> LS by Wavelength
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotCIE (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OI');
oiPlot(oi, 'chromaticity roi');
return

% --------------------------------------------------------------------
function menuPlotNewGraphWin_Callback(hObject, eventdata, handles)
% (Window) Open a new window
%
% Syntax:
%   menuPlotNewGraphWin_Callback(hObject, eventdata, handles)
%
% Description:
%    Open a new window using vcNewGraphWin.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotNewGraphWin (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
vcNewGraphWin;
return;

% --------------------------------------------------------------------
function menuPlOp_Callback(hObject, eventdata, handles)
% (Plot | Optics) Optics sub-heading of Plot
%
% Syntax:
%   menuPlOp_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot option. Plot -> Optics
%
% Inputs:
%    hObject   - Handle. Handle to menuPlOp (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --------------------------------------------------------------------
function menuTransmittance_Callback(hObject, eventdata, handles)
% (Analyze | Optics | Transmittance) Object transmittance
%
% Syntax:
%   menuTransmittance_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze -> Optics -> Transmittance
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotLSWave (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
opticsPlotTransmittance(vcGetObject('OPTICALIMAGE'));
return;

% --------------------------------------------------------------------
function menuAnPSFMovie_Callback(hObject, eventdata, handles)
% (Analyze | Optics | PSF Movie) PSF Movie Analysis
%
% Syntax:
%   menuAnPSFMovie_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze -> Optics -> PSF Movie
%
% Inputs:
%    hObject   - Handle. Handle to menuAnPSFMovie (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
psfMovie;
return;

% --------------------------------------------------------------------
function menuPlotPS550_Callback(hObject, eventdata, handles)
% (Analyze | Optics | PSF Mesh (550)) Analyze PSF Mesh at 550nm
%
% Syntax:
%   menuPlotPS550_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze -> Optics -> PSF Mesh 550nm
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotPS550 (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OPTICALIMAGE');
oiPlot(oi, 'psf 550');
return;

% --------------------------------------------------------------------
function menuPlotLSWave_Callback(hObject, eventdata, handles)
% (Analyze | Optics | LS by Wavelength) LS by Wavelength Analysis
%
% Syntax:
%   menuPlotLSWave_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze -> Optics -> LS by Wavelength
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotLSWave (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OPTICALIMAGE');
oiPlot(oi, 'ls Wavelength');
return;

% --------------------------------------------------------------------
function menuPlOTFWave_Callback(hObject, eventdata, handles)
% (Analyze | Optics | OTF 1d by wave) Analyze OTF 1DxWave
%
% Syntax:
%   menuPlOTFWave_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze the OTF (1D) by a wavelength.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlOTFWave (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OPTICALIMAGE');
oiPlot(oi, 'otf Wavelength');

return;

% --------------------------------------------------------------------
function menuOTFAnyWave_Callback(hObject, eventdata, handles)
% (Analyze | Optics | OTF) User selects wavelength & plots OTF
%
% Syntax:
%   menuOTFAnyWave_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze -> Optics -> OTF. User selects wavelength for
%    the OTF calculation.
%
% Inputs:
%    hObject   - Handle. Handle to menuOTFAnyWave (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OPTICALIMAGE');
oiPlot(oi, 'otf');
return;

% --------------------------------------------------------------------
function plotOTF_Callback(hObject, eventdata, handles)
% (Analyze | Optics | OTF (550)) The OTF 550nm option inside Optics
%
% Syntax:
%   menuPlotOTF_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze the OTF at 550nm in the Optics sub-menu.
%
% Inputs:
%    hObject   - Handle. Handle to plotOTF (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OPTICALIMAGE');
oiPlot(oi, 'otf 550');
return;

% --------------------------------------------------------------------
function menuPlotOffAxis_Callback(hObject, eventdata, handles)
% (Analyze | Optics | Off-Axis fall-off) Analyze the off-axis fall-off
%
% Syntax:
%   menuPlotOffAxis_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze the off-axis fall-off
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotOffAxis (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OPTICALIMAGE');
opticsPlotOffAxis(oi);  % If  no ray trace, cos4th.
return;

% --------------------------------------------------------------------
function menuPlCIE_Callback(hObject, eventdata, handles)
% ( | ) Short Description
%
% Syntax:
%   menuP1CIE_Callback(hObject, eventdata, handles)
%
% Description:
%    Longer Description.
%
% Inputs:
%    hObject   - Handle. Handle to menuP1CIE_Callback(see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --------------------------------------------------------------------
function menuAn_Callback(hObject, eventdata, handles)
% (Analyze) Menu Heading for Analyze
%
% Syntax:
%   menuAn_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Top level for Analyze
%
% Inputs:
%    hObject   - Handle. Handle to menuAnLineIllumVertical (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --------------------------------------------------------------------
function menuAnOpticsPSF_Callback(hObject, eventdata, handles)
% (Analyze | Optics | PSF) Analyze the PSF of the Optics
%
% Syntax:
%   menuAnOpticsPSF_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze the PSF of the Optics
%
% Inputs:
%    hObject   - Handle. Handle to menuAnOpticsPSF (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

oi  = vcGetObject('OPTICALIMAGE');
oiPlot(oi, 'psf');

return;

% --------------------------------------------------------------------
function menuAnalyzeLinePlots_Callback(hObject, eventdata, handles)
% (Analyze | Plots | Line) Analyze the line plots
%
% Syntax:
%   menuAnalyzeLinePlots_Callback(hObject, eventdata, handles)
%
% Description:
%    Analyze option. Analyze the Line plots
%
% Inputs:
%    hObject   - Handle. Handle to menuAnalyzeLinePlots (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --------------------------------------------------------------------
function menuAnLineIllum_Callback(hObject, eventdata, handles)
% (Plot | Illuminance) Plot line illuminance
%
% Syntax:
%   menuAnLineIllum_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Plot line illuminance.
%
% Inputs:
%    hObject   - Handle. Handle to menuAnLineIllum (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --------------------------------------------------------------------
function menuAnLineIllumHorizontal_Callback(hObject, eventdata, handles)
% (Plot | Illuminance | Horizontal) Plot a horizontal line illuminance
%
% Syntax:
%   menuAnLineIllumHorizontal_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Plot a Horizontal line for the illuminance.
%
% Inputs:
%    hObject   - Handle. Handle to menuAnLineIllumHorizontal (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiPlot(vcGetObject('OI'), 'illuminance hline');
return;

% --------------------------------------------------------------------
function menuAnLineIllumVertical_Callback(hObject, eventdata, handles)
% (Plot | Illuminance | Vertical) Plot a vertical line illuminance
%
% Syntax:
%   menuAnLineIllumVertical_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Plot a Veritcal line for the illuminance.
%
% Inputs:
%    hObject   - Handle. Handle to menuAnLineIllumVertical (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiPlot(vcGetObject('OI'), 'illuminance vline');
return;

% --------------------------------------------------------------------
function menuAnLineIllumHorFFT_Callback(hObject, eventdata, handles)
% (Plot | Illuminance | Horizontal | FFT) Plot illuminance hline FFT
%
% Syntax:
%   menuAnLineIllumHorFFT_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Plot a Horizontal Line FFT for the illuminance.
%
% Inputs:
%    hObject   - Handle. Handle to menuAnLineIllumHorFFT (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiPlot(vcGetObject('OI'), 'illuminance fft hline');
return;

% --------------------------------------------------------------------
function menuAnLineIllumVertFFT_Callback(hObject, eventdata, handles)
% (Plot | Illuminance | Vertical | FFT) Plot illuminance vertical line FFT
%
% Syntax:
%   menuAnLineIllumVertFFT_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Plot a Veritcal line FFT for the illuminance.
%
% Inputs:
%    hObject   - Handle. Handle to menuAnLineIllumVertFFT (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oiPlot(vcGetObject('OI'), 'illuminance fft vline');
return;

% --------------------------------------------------------------------
function menuAnOptSampling_Callback(hObject, eventdata, handles)
% (Analyze | Sampling Information) Returns Sampling Cutoff & Max Freq
%
% Syntax:
%   menuAnOptSampling_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. The image sampling rate supports a certain spatial
%    frequency. The diffraction limited optics supports a certain spatial
%    frequency. We only obtain a very accurate spatial representation when
%    the image sampling supports a representation as high as the
%    diffraction limited optics. Otherwise, the higher spatial frequencies
%    are not represented in the result.
%
%    In many cases, people will leave the lower sampling rate, which
%    provides speed but blurs the image, because they are interested in
%    other features of the simulation.
%
% Inputs:
%    hObject   - Handle. Handle to menuVLine (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('oi');
inCutoff = opticsGet(oiGet(oi, 'optics'), 'maxincutoff', 'mm');
maxFres = oiGet(oi, 'maxFreqRes', 'mm');
str = sprintf('DL cutoff: %.2f - Samp cutoff %.2f (cyc/mm)\n', ...
    inCutoff, maxFres);
ieInWindowMessage(str, handles);
return;

% --------------------------------------------------------------------
function menuROISummaries_Callback(hObject, eventdata, handles)
% (Analyze | ROI Summaries) ROI Summaries menu in Analyze
%
% Syntax:
%   menuROISummaries_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Analyzes -> ROI Summaries.
%
% Inputs:
%    hObject   - Handle. Handle to menuROISummaries (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --------------------------------------------------------------------
function menuPlotLuxHist_Callback(hObject, eventdata, handles)
% (Analyze | ROI Summaries | Illuminance) Plot the illuminance ROI
%
% Syntax:
%   menuPlotLuxHist_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Analyzes and Plot the ROI Summary & Illuminance.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotLuxHist (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OI');
oiPlot(oi, 'illuminance roi');
return;

% --------------------------------------------------------------------
function menuPlotRGB_Callback(hObject, eventdata, handles)
% (Plot | Image | RGB) Plots the current RGB image in a separate window.
%
% Syntax:
%   menuPlotRGB_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Plots the current RGB image in a separate window.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotRGB (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
imageMultiview('oi', vcGetSelectedObject('oi'));
return;

% --------------------------------------------------------------------
function menuPlotMultiRGB_Callback(hObject, eventdata, handles)
% (Plot | Image | RGB | Multiple) Plot selected RGB images from session OIs
%
% Syntax:
%   menuPlotMultiRGB_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot option. Plot the selected RGB images from all of the OIs in the
%    current session.
%
% Inputs:
%    hObject   - Handle. Handle to menuPlotMultiRGB (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
imageMultiview('oi');
return;

% --------------------------------------------------------------------
function menuHline_Callback(hObject, eventdata, handles)
% (Analyze | Line | Horizontal) Calls the Horizontal Line
%
% Syntax:
%   menuHLine_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Analyzes the Horizontal option of Line.
%
% Inputs:
%    hObject   - Handle. Handle to menuHLine (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OI');
oiPlot(oi, 'hline');
return;

% --------------------------------------------------------------------
function menuVLine_Callback(hObject, eventdata, handles)
% (Analyze | Line | Vertical) Calls the Vertical Line
%
% Syntax:
%   menuVLine_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Analyzes the Vertical option of Line.
%
% Inputs:
%    hObject   - Handle. Handle to menuVLine (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OI');
oiPlot(oi, 'vline');
return;

% --------------------------------------------------------------------
function menuFFTamp_Callback(hObject, eventdata, handles)
% (Analyze | FFT2dAmp) Analyze FFT2dAmp. Default image & mid wavelen
%
% Syntax:
%   menuFFTamp_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Clicking analyzes the FFT2dAmp. The default returns the
%    whole image and a middle wavelength.
%
% Inputs:
%    hObject   - Handle. Handle to menuFFTamp (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OI');
oiPlot(oi, 'irradiance fft');
return;

% --------------------------------------------------------------------
function menuStandForm_Callback(hObject, eventdata, handles)
% (Optics | Standard Format) Calls the Optics StandardFormat
%
% Syntax:
%   menuStandForm_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Clicking calls the optics standardFormat
%
% Inputs:
%    hObject   - Handle. Handle to menuStandForm (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return;

% --- Executes during object creation, after setting all properties.
function popCustom_CreateFcn(hObject, eventdata, handles)
% (Create | PopUp | Custom) Creates a custom popUp using defaults
%
% Syntax:
%   popCustom_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    PopUp option. Create a custom popUp using the defaults.
%
% Inputs:
%    hObject   - Handle. Handle to popCustom (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - empty - The handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

% --- Select type of optics model
function popOpticsModel_CreateFcn(hObject, eventdata, handles)
% (Create | PopUp | Optics Model Create) Creates Optics Model & Select Type
%
% Syntax:
%   popOpticsModel_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    PopUp Option. Creates the Optics model & selects type, following the
%    creation defaults.
%
% Inputs:
%    hObject   - Handle. Handle to popOpticsModel (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - empty - The handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end
return;

% --- Interpret the popup call back
function popOpticsModel_Callback(hObject, eventdata, handles)
% (PopUp | Optics Model) Calls the Optics Model & Interprets the popUp
%
% Syntax:
%   popOpticsModel_Callback(hObject, eventdata, handles)
%
% Description:
%    PopUp option. Interprets the optics model popUp callback
%
% Inputs:
%    hObject   - Handle. Handle to popOpticsModel (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
contents = get(handles.popOpticsModel, 'String');

% The method names are in the GUI of the window. As of June 06, 2010 the
% options were: Diffraction-limited, Shift-invariant, Ray trace, Skip OTF
method = contents{get(handles.popOpticsModel, 'Value')};

oi = vcGetObject('oi');
optics = oiGet(oi, 'optics');

switch lower(method)
    case 'diffraction-limited'
        optics = opticsSet(optics, 'model', 'diffractionLimited');
    case 'shift-invariant'
        optics = opticsSet(optics, 'model', 'shiftInvariant');
        if isempty(opticsGet(optics, 'otfdata'))
            % Warn the user
            ieInWindowMessage('Shift-invariant OTF data not loaded.', ...
                handles, 2);
            disp('Shift-invariant data not loaded')
        end
    case 'iset3d'
        optics = opticsSet(optics,'model','iset3d');

    otherwise
        error('Unknown optics method');
end

oi = oiSet(oi, 'optics', optics);
vcReplaceObject(oi);
oiRefresh(hObject, eventdata, handles);
return;

% --- Executes on selection change in popCustom.
function popCustom_Callback(hObject, eventdata, handles)
% (PopUp | Custom) Calls a custom render popUp.
%
% Syntax:
%   popCustom_Callback(hObject, eventdata, handles)
%
% Description:
%    PopUp option. Calls a custom render popUp. This is initialized with
%    vcimageRender, Add Custom, Delete Custom, -----. We will add and
%    delete routines from vcSESSION.CUSTOM.procMethod list.
%
% Inputs:
%    hObject   - Handle. Handle to popCustom (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
contents = get(handles.popCustom, 'String');
method = contents{get(handles.popCustom, 'Value')};

[val, oi] = vcGetSelectedObject('OPTICALIMAGE');

oi = oiSet(oi, 'oiMethod', method);

oi = oiSet(oi, 'consistency', 0);
vcReplaceObject(oi, val);
oiRefresh(hObject, eventdata, handles);
return;

% --- Executes on button press in btnOffAxis.
function btnOffAxis_Callback(hObject, eventdata, handles)
% (Button | Off Axis) The Off Axis button. Toggles cos4th (on or off).
%
% Syntax:
%   btnOffAxis_Callback(hObject, eventdata, handles)
%
% Description:
%    Button option. Clicking the radio button toggles the cos4th between on
%    and off.
%
% Inputs:
%    hObject   - Handle. Handle to btnOffAxis (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, oi] = vcGetSelectedObject('OPTICALIMAGE');
optics = oiGet(oi, 'optics');
if get(hObject, 'Value')
    optics = opticsSet(optics, 'offaxismethod', 'cos4th');
    ieInWindowMessage([], handles, []);
else
    optics = opticsSet(optics, 'offaxismethod', 'skip');
    ieInWindowMessage([], handles, []);
end
oi = oiSet(oi, 'optics', optics);
vcReplaceObject(oi, val);
oiRefresh(hObject, eventdata, handles);
return;

% --- Executes during object creation, after setting all properties.
function popDiffuser_CreateFcn(hObject, eventdata, handles)
% (Create | PopUp | Diffuser) Creates the diffuser popUp using defaults.
%
% Syntax:
%   popDiffuser_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    PopUp option. Create the diffuser's popUp
%
% Inputs:
%    hObject   - Handle. Handle to popDiffuser (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - empty - handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: popupmenu controls usually have a white background on Windows.
%      See ISPC and COMPUTER.
%
if ispc && ...
        isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end
return;

% --- Executes on selection change in popDiffuser.
function popDiffuser_Callback(hObject, eventdata, handles)
% (PopUp | Diffuser) Selects the Diffuser method using a popUp.
%
% Syntax:
%   popDiffuser_Callback(hObject, eventdata, handles)
%
% Description:
%    PopUp option. PopUp allows user to select the diffuser method. Options
%    include: skip, blur, and birefringent.
%
% Inputs:
%    hObject   - Handle. Handle to menuHelpAppNotes (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, oi] = vcGetSelectedObject('OPTICALIMAGE');
contents = get(handles.popDiffuser, 'String');
dMethod = contents{get(handles.popDiffuser, 'Value')};

oi = oiSet(oi, 'diffuserMethod', dMethod);
vcReplaceObject(oi, val);
oiRefresh(hObject, eventdata, handles);
return;

% --- Executes on button press in btnDiffuser.
function btnDiffuser_Callback(hObject, eventdata, handles)
% (Button | Diffuser) Toggles the diffuser simulation
%
% Syntax:
%   btnDiffuser_Callback(hObject, eventdata, handles)
%
% Description:
%    Button option. Clicking toggles the diffuser simulation (on, off).
%
%    *May be obsolete?? Possibly replaced by the popup, popDiffuser.*
%
% Inputs:
%    hObject   - Handle. Handle to menuHelpAppNotes (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, oi] = vcGetSelectedObject('OPTICALIMAGE');
if get(hObject, 'Value')
    oi = oiSet(oi, 'diffuserMethod', 'blur');
    ieInWindowMessage([], handles, []);
else
    oi = oiSet(oi, 'diffuserMethod', 'skip');
    ieInWindowMessage([], handles, []);
end

vcReplaceObject(oi, val);
oiRefresh(hObject, eventdata, handles);
return;

% --- Executes during object creation, after setting all properties.
function editDiffuserBlur_CreateFcn(hObject, eventdata, handles)
% (Create | Diffuser Blur) Creates the Diffuser using defaults
%
% Syntax:
%   editDiffuserBlur_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Menu option. Clicking calls the web page help file.
%
% Inputs:
%    hObject   - Handle. Handle to menuHelpAppNotes (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

function editDiffuserBlur_Callback(hObject, eventdata, handles)
% (Edit | Diffuser Blur Call) Set the FWHM (in um) of the Diffuser
%
% Syntax:
%   editDiffuserBlur_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit option. Edit the Diffuser's FWHM in micrometers.
%
% Inputs:
%    hObject   - Handle. Handle to editDiffuserBlur (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, oi] = vcGetSelectedObject('OPTICALIMAGE');

% returns contents of editDiffuserBlur as a double
blur = str2double(get(hObject, 'String'));

oi = oiSet(oi, 'diffuserBlur', blur * 10 ^ -6);  % Stored in meters
vcReplaceObject(oi, val);
oiRefresh(hObject, eventdata, handles);
return;

% --- Executes on button press in btnOTF.
function btnOTF_Callback(hObject, eventdata, handles)
% (Button | Diff. Ltd. OTF) Button for the diffraction limited OTF
%
% Syntax:
%   btnOTF_Callback(hObject, eventdata, handles)
%
% Description:
%    Button. Pressing calls diffraction limited OTF.
%
% Inputs:
%    hObject   - Handle. Handle to btnOTF (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[val, oi] = vcGetSelectedObject('OPTICALIMAGE');
optics = oiGet(oi, 'optics');

if get(hObject, 'Value')
    optics = opticsSet(optics, 'otfmethod', 'dlMTF');
    ieInWindowMessage([], handles, []);
else
    optics = opticsSet(optics, 'otfmethod', 'skip');
    ieInWindowMessage([], handles, []);
end
oi = oiSet(oi, 'optics', optics);
vcReplaceObject(oi, val);
oiRefresh(hObject, eventdata, handles);
return;

% --------------------------------------------------------------------
function menuHelp_Callback(hObject, eventdata, handles)
% (Help | Top) Top menu for help
%
% Syntax:
%   menuHelp_Callback(hObject, eventdata, handles)
%
% Description:
%    Top Menu option. Clicking opens the help menu
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
return;

% --------------------------------------------------------------------
function menuHelpAppNotes_Callback(hObject, eventdata, handles)
% (Help | Documentation (Web)) Calls the online help text
%
% Syntax:
%   menuHelpAppNotes_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Clicking calls the web page help file.
%
% Inputs:
%    hObject   - Handle. Handle to menuHelpAppNotes (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
web('https://github.com/isetbio/isetbio/wiki', '-browser');
return;

% --------------------------------------------------------------------
function menuHelpOpticsOnline_Callback(hObject, eventdata, handles)
% (Help | Optics Online) Calls the online help text for Optics
%
% Syntax:
%   menuHelpOpticsOnline_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Clicking calls the web page help file for Optics.
%
% Inputs:
%    hObject   - Handle. Handle to menuHelpOpticsOnline (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
web('https://github.com/isetbio/isetbio/wiki', '-browser');
return;

% --------------------------------------------------------------------
function menuHelpOIOnline_Callback(hObject, eventdata, handles)
% (Help | Optical Image Functions) Calls the online help text for OI
%
% Syntax:
%   menuHelpOIOnline_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Clicking calls the web page help file for OI.
%
% Inputs:
%    hObject   - Handle. Handle to menuHelpOIOnline (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
web('https://github.com/isetbio/isetbio/wiki', '-browser');
return;

% --------------------------------------------------------------------
function menuHelpISETOnline_Callback(hObject, eventdata, handles)
% (Help | ISET Functions) Calls the online help text for ISET
%
% Syntax:
%   menuHelpISETOnline_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option. Clicking calls the web page help file for ISET.
%
% Inputs:
%    hObject   - Handle. Handle to menuHelpISETOnline (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
web('https://github.com/isetbio/isetbio/wiki', '-browser');
return;


% --- Executes on button press in btnNext.
function btnNext_Callback(hObject, eventdata, handles)
% Button to move to next image (->)
%
% Syntax:
%   btnPrev_Callback(hObject, eventdata, handles)
%
% Description:
%    The call to the button to navigate to the next image.
%
% Inputs:
%    hObject   - Handle. Handle to btnNext (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
thisOI = ieSessionGet('selected', 'oi');
nS  = ieSessionGet('nobjects', 'oi');
thisOI = min(thisOI + 1, nS);
thisOI = max(thisOI, 1);
vcSetSelectedObject('oi', thisOI);
oiRefresh(hObject, eventdata, handles);
return;

% --- Executes on button press in btnPrev.
function btnPrev_Callback(hObject, eventdata, handles)
% Button to move to previous image (<-)
%
% Syntax:
%   btnPrev_Callback(hObject, eventdata, handles)
%
% Description:
%    The call to the button to navigate to the previous image.
%
% Inputs:
%    hObject   - Handle. Handle to btnPrev (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
thisOI  = ieSessionGet('selected', 'oi');
nS = ieSessionGet('nobjects', 'oi');
thisOI = min(thisOI - 1, nS);
thisOI = max(thisOI, 1);
vcSetSelectedObject('oi', thisOI);
oiRefresh(hObject, eventdata, handles);
return;