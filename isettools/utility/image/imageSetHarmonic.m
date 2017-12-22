function varargout = imageSetHarmonic(varargin)
% UI to read in the parameters of a harmonic pattern.
%
% Syntax:
%   varargout = imageSetHarmonic(varargin)
%
% Description:
%    This functions is the M-file for imageSetHarmonic.fig
%
%	 imageSetHarmonic, by itself, creates a new IMAGEHARMONIC or raises
%    the existing singleton*.
%
%    H = imageSetHarmonic returns the handle to a new IMAGEHARMONIC or the
%    handle to the existing singleton*.
%
%    Parameters:
%       parms.freq, parms.contras, parms.ph, parms.ang, parms.row,
%       parms.col, parms.GaborFlag
%
%	 imageSetHarmonic('CALLBACK', hObject, eventData, handles, ...) calls
%	 the local function named CALLBACK in IMAGESETHARMONIC.M with the given
%    input arguments.
%
%    imageSetHarmonic('Property', 'Value', ...) creates a new
%    IMAGESETHARMONIC or raises the existing singleton*. Starting from
%    the left, property value pairs are applied to the GUI before
%    imageSetHarmonic_OpeningFunction gets called. An unrecognized
%    property name or invalid value makes property application stop. All
%    inputs are passed to imageSetHarmonic_OpeningFcn via varargin.
%
%	 *See GUI Options on GUIDE's Tools menu. Choose "GUI allows only one
%    instance to run (singleton)".
%
% Inputs:
%    varargin - Variable list of input parameters.
%
% Outputs:
%    varargout - Variable list of output parameters.
%
% Notes:
%    * [Note: XXX - I have not found a good way to return the parameters of
%      the harmonic function. These are frequency, contrast, phase, angle,
%      Gaussian window, and row, col. To get them back, I create global
%      parms, read it, and then destroy it. This is really dumb. There must
%      be a better way.]
%
% See Also:
%    guide, guiData, guiHandles
%

% History:
%    xx/xx/03         Copyright ImagEval Consultants, LLC, 2003.
%    08/23/03  GUIDE  v2.5 09:34:24 Edit the above text to modify the
%                     response to help imageSetHarmonic
%    12/11/17  jnm    Formatting

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
                   'gui_Singleton', gui_Singleton, ...
                   'gui_OpeningFcn', @imageSetHarmonic_OpeningFcn, ...
                   'gui_OutputFcn', @imageSetHarmonic_OutputFcn, ...
                   'gui_LayoutFcn', [] , ...
                   'gui_Callback', []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
return;

% --- Executes just before imageSetHarmonic is made visible.
function imageSetHarmonic_OpeningFcn(hObject, ~, handles, varargin)
% Opening function for imageSetHarmonic
%
% Syntax:
%   imageSetHarmonic_OpeningFcn(hObject, [], handles, [varargin])
%
% Description:
%    Function executes just before imageSetHarmonic becomes visible. The
%    opening function for imageSetHarmonic to determine the default command
%    line output.
%
% Inputs:
%    hObject  - The Harmonic object
%    ~        - Empty space?
%    handles  - The handle associated with the object
%    varargin - Optional variable input array for function
%
% Outputs:
%    None.
%

% Choose default command line output for imageSetHarmonic
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imageSetHarmonic wait for user response (see UIRESUME)
% uiwait(handles.figure1);
return;

% --- Outputs from this function are returned to the command line.
function varargout = imageSetHarmonic_OutputFcn(hObject, eventdata, ...
    handles)
% Create outputs for command line
%
% Syntax:
%   varargout = imageSetHarmonic_OutputFcn(hObject, eventdata, handles)
%
% Description:
%    Capture imageSetHarmonic function output for the command line.
%
% Inputs:
%    hObject   - Harmonic Object
%    eventdata - Function event data
%    handles   - The handle associated with the object
%
% Outputs:
%    varargout - Output arguments of variable length. Specified below as a
%                two-argument list. The first argument is the handle output
%                property. The second argument is the output from the
%                btnDone_Callback function.
%
varargout{1} = handles.output;
varargout{2} = btnDone_Callback(hObject, eventdata, handles, 0);
return;

% --- Executes during object creation, after setting all properties.
function editFreq_CreateFcn(hObject, ~, ~)

if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

% --- Executes during object creation, after setting all properties.
function editContrast_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

% --- Executes during object creation, after setting all properties.
function editRow_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

% --- Executes during object creation, after setting all properties.
function editCol_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

% --- Executes during object creation, after setting all properties.
function editPh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

% --- Executes during object creation, after setting all properties.
function editAng_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject, 'BackgroundColor', 'white');
else
    set(hObject, 'BackgroundColor', ...
        get(0, 'defaultUicontrolBackgroundColor'));
end
return;

% --- Executes on button press in btnDone.
function params = btnDone_Callback(hObject, eventdata, handles, closeMe)

if notDefined('closeMe'), closeMe = 1; end

params.freq = str2double(get(handles.editFreq, 'String'));
params.contrast = str2double(get(handles.editContrast, 'String'));
params.ph = str2double(get(handles.editPh, 'String'));
params.ang = str2double(get(handles.editAng, 'String'));
params.row = str2double(get(handles.editRow, 'String'));
params.col = str2double(get(handles.editCol, 'String'));
params.GaborFlag = get(handles.btnGabor, 'Value');

if closeMe, close(gcbf); end

return;

function editFreq_Callback(hObject, eventdata, handles), return;

function editContrast_Callback(hObject, eventdata, handles), return;

function editPh_Callback(hObject, eventdata, handles), return;

function editAng_Callback(hObject, eventdata, handles), return;

function editRow_Callback(hObject, eventdata, handles), return;

function editCol_Callback(hObject, eventdata, handles), return;

% --- Executes on button press in btnGabor.
function btnGabor_Callback(hObject, eventdata, handles), return;