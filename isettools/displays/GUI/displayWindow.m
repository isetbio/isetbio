function varargout = displayWindow(varargin)
% displayWindow main window
%
% Syntax:
%   [varargout] = displayWindow([varargin])
%
% Description:
%    This is the main GUI window for interfacing with the Clear Type or
%    Display Simulator design functions.  From this window you can
%    visualize the sub-pixels, load simple images, and perform various
%    analytical calculations using the display.  The display radiance data
%    can also be converted into an ISET Scene format and thus transferred
%    into the ISET analysis tools.
%
%    This function brings up the window to edit display properties
%
%       DISPLAYWINDOW, by itself, creates a new DISPLAYWINDOW or raises the
%       existing singleton*.
%
%       H = DISPLAYWINDOW returns the handle to a new DISPLAYWINDOW or the
%       handle to the existing singleton*.
%
%       DISPLAYWINDOW('Property', 'Value', ...) creates a new DISPLAYWINDOW
%       using the given property value pairs. Unrecognized properties are
%       passed via varargin to displayWindow_OpeningFcn.  This calling
%       syntax produces a warning when there is an existing singleton*.
%
%       DISPLAYWINDOW('CALLBACK') & DISPLAYWINDOW('CALLBACK', hObject, ...)
%       call the local function named CALLBACK in DISPLAYWINDOW with the
%       given input arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% Inputs:
%    None required.
%
% Outputs:
%    None required.
%
% Optional key/value pairs:
%    **NEEDS TO BE FILLED OUT**
%
% See Also:
%    GUIDE, GUIDATA, GUIHANDLES
%

% History:
%    xx/xx/10  BW   (c) Stanford, PDCSOFT, Wandell, 2010
%    xx/xx/14  HJ   HJ, PDCSOFT, 2014
%    05/18/18  jnm  Formatting

%#ok<*DEFNU> % suppress unused warnings

% Edit the above text to modify the response to help displayWindow

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
                   'gui_Singleton', gui_Singleton, ...
                   'gui_OpeningFcn', @displayWindow_OpeningFcn, ...
                   'gui_OutputFcn', @displayWindow_OutputFcn, ...
                   'gui_LayoutFcn', [], ...
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

% --- Executes just before displayWindow is made visible.
function displayWindow_OpeningFcn(hObject, eventdata, handles, varargin)
% Set defaults for displayWindow
%
% Syntax:
%   displayWindow_OpeningFcn(hObject, eventdata, handles, [varargin])
%
% Description:
%    Set up defaults for displayWindow just before it is made visible.
%
% Inputs:
%    hObject   - handle to displayWindow (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    **NEEDS TO BE ADDED**
%

% Update handles structure
% Not yet properly implemented for the display.
%
% handles.output = hObject;
% guidata(hObject, handles);
% vcSetFigureHandles('DISPLAY',hObject,eventdata,handles);
%

%  Check the preferences for ISET and adjust the font size.
ieFontInit(hObject);

if ~isempty(varargin)
    thisDisplay = varargin{1};
    if strcmp(thisDisplay.type,'display'), ieAddObject(thisDisplay);
    else, warning('Unexpected variable input.\n');
    end
end

handles.output = hObject;

% Refresh handles structure
guidata(hObject, handles);

% add to vcSESSION
global vcSESSION
if ~checkfields(vcSESSION, 'GUI', 'vcDisplayWindow')
    vcSESSION.GUI.vcDisplayWindow.hObject = hObject;
    vcSESSION.GUI.vcDisplayWindow.eventdata = eventdata;
    vcSESSION.GUI.vcDisplayWindow.handles = handles;
end

% Refresh image window
if ~isfield(vcSESSION, 'imgData') || isempty(vcSESSION.imgData)
    I = imread(fullfile(isetbioDataPath, 'images', 'rgb', 'macbeth.tif'));
    I = im2double(I);
    vcSESSION.imgData = I;
end

% Refresh other components
displayRefresh(hObject, eventdata, handles);

return;

% --- Outputs from this function are returned to the command line.
function varargout = displayWindow_OutputFcn(~, ~, handles)
% Get default command line output from handles structure
%
% Syntax:
%   [varargout] = displayWindow_OutputFcn(~, ~, handles)
%
% Description:
%    Get the default command line output from the handles structure.
%
% Inputs:
%    hObject   - handle to displayWindow (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Struct. Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    **NEEDS TO BE ADDED**
%
% Optional key/value pairs:
%    None.
%
varargout{1} = handles.output;

return;

% --------------------------------------------------------------------
function menuFileLoadImage_Callback(hObject, eventdata, handles)
% (File | Load Image) Load the image file
%
% Syntax:
%   menuFileLoadImage_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu under file to load an image.
%
% Inputs:
%    hObject   - handle to menuFileLoadImage (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
[fname, p] = uigetfile('*.*', 'Choose Image File');
if fname == 0, return; end
I = im2double(imread(fullfile(p, fname)));

ind = vcGetSelectedObject('display');
if isempty(ind), disp('No display selected'); return; end
d = vcGetObject('display', ind);
d = displaySet(d, 'main image', I);
d_list = vcGetObjects('display');
d_list{ind} = d;
vcSetObjects('display', d_list);

% Refresh other components
displayRefresh(hObject, eventdata, handles);

return;

% --------------------------------------------------------------------
function menuEditNew_Callback(~, ~, handles)
% (Edit | New Model) Create a new display model
%
% Syntax:
%   menuEditNew_Callback(~, ~, handles)
%
% Description:
%    Menu option under Edit to create a new display model.
%
% Inputs:
%    hObject   - handle to menuEditNew (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
d = displayCreate('LCD-Apple');
vcAddAndSelectObject('display', d);
displayRefresh([], [], handles);
return;

% --------------------------------------------------------------------
function menuDeleteCurrent_Callback(~, ~, handles)
% (Edit | Delete) Delete a display
%
% Syntax:
%   menuDeleteCurrent_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option under Edit to delete the current display.
%
% Inputs:
%    hObject   - handle to menuDeleteCurrent (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
vcDeleteSelectedObject('display');
displayRefresh([], [], handles);
return;

% --------------------------------------------------------------------
function subpixeImage_Callback(~, ~, ~)
% (Plot | Sub pixel image) Plot Subpixel PSF
%
% Syntax:
%   subpixeImage_Callback(~, ~, ~)
%
% Description:
%    Menu option under Plot to plot the subpixel PSF
%
% Inputs:
%    hObject   - handle to subpixeImage (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ind = vcGetSelectedObject('display');
if isempty(ind), disp('No display selected'); return; end
d = vcGetObject('display', ind);
psfs = displayGet(d, 'dixel image');
if isempty(psfs)
    disp('no psfs in display model');
    return;
else
    vcNewGraphWin;
    if size(psfs, 3) == 3
        imshow(psfs / max(psfs(:)));
    else
        gam = 1;
        wave = displayGet(d, 'wave');
        photons = vcReadImage(psfs, 'rgb', d);
        imageSPD(photons, wave, gam, [], [], 1);
    end
end

return;

% --------------------------------------------------------------------
function menuAnalyzeOutputImage_Callback(~, ~, handles)
% (Analyze | Output image (by units)) The output image by selected units
%
% Syntax:
%   menuAnalyzeOutputImage_Callback(~, ~, handles)
%
% Description:
%    Menu option under Analyze to calculate the output image by units. Will
%    analyze by whichever units are selected.
%
% Inputs:
%    hObject   - handle to menuAnalyzeOutputImage (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ind = vcGetSelectedObject('display');
d = vcGetObject('display', ind);
if isempty(ind), disp('No display selected'); return; end

h = get(handles.axes1, 'Children');
I = get(h, 'CData');

% select region
set(handles.txtMessage, 'String', 'Select Region in Image');
rect = floor(getrect);
set(handles.txtMessage, 'String', 'Original Image');

% round width and height to multiple of pixelsperdixel
ppd = displayGet(d, 'pixels per dixel');
rect(3:4) = ppd .* ceil(rect(3:4) ./ ppd);
I = I(rect(2): rect(2) + rect(4) - 1, rect(1) : rect(1) + rect(3) - 1, :);

% Create a scene:
% For the rgb image, we can use displayCompute the generate the output
% image However, the output image will have same number of primaries as the
% display (could be more than 3) and then imshow will not work. Thus, we
% compute the scene radiance and get the srgb of the scene and show it to
% the window.
scene = sceneFromFile(I, 'rgb', [], d, [], 1, [], [20 20]);
vcNewGraphWin;
imshow(sceneGet(scene, 'rgb image'));

return;

% --------------------------------------------------------------------
function menuAnalyzeSceneSubpixel_Callback(~, ~, handles)
% (Analyze | Scene) Analyze scene by subpixel
%
% Syntax:
%   menuAnalyzeSceneSubpixel_Callback(~, ~, handles)
%
% Description:
%    Menu option under Analyze for the subpixel analysis of the scene.
%
% Inputs:
%    hObject   - handle to menuAnalyzeSceneSubpixel (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Opens up a Scene Window
ind = vcGetSelectedObject('display');
d = vcGetObject('display', ind);

% Get current image
I = displayGet(d, 'main image');
if isempty(I)
    % Should go away
    global vcSESSION;
    if isfield(vcSESSION, 'imgData')
        I = displayGet(d, 'main image');
    else
        warning('No image set');
        return;
    end
end

% select region
set(handles.txtMessage, 'String', 'Select Region in Image');
rect = floor(getrect);
set(handles.txtMessage, 'String', 'Original Image');
% round width and height to multiple of pixelsperdixel
ppd = displayGet(d, 'pixels per dixel');
rect(3:4) = ppd .* ceil(rect(3:4) ./ ppd);
I = I(rect(2): rect(2) + rect(4) - 1, rect(1) : rect(1) + rect(3) - 1, :);

% Get scene info
overSample = displayGet(d, 'over sample');
answer = inputdlg({'Scene name', 'Up sample (col)', 'Up sample (row)'}, ...
    'Scene Info', 1, {sprintf('scene-%s', displayGet(d, 'name')), ...
    num2str(overSample(1)), num2str(overSample(2))});
if isempty(answer)
    return;
else
    name = answer{1};
    overSample = [str2double(answer{2}) str2double(answer{3})];
end

% Generate scene
if isempty(I), disp('No image set'); return; end
scene = sceneFromFile(I, 'rgb', [], d, [], 1, [], overSample);
scene = sceneSet(scene, 'name', name);

% Compute scene size
[r, c, ~] = size(I);
vDist = sceneGet(scene, 'distance');
fov = atand(max(r, c) * displayGet(d, 'metersperdot') / vDist);
scene = sceneSet(scene, 'fov', fov);

% Show scene window
vcAddAndSelectObject('scene', scene);
sceneWindow;

% --------------------------------------------------------------------
function menuAnalyzeScene_Callback(~, ~, ~)
% (Analyze | Scene Pixel) Call to analyze a scene.
%
% Syntax:
%   menuAnalyzeScene_Callback(~, ~, ~)
%
% Description:
%    Menu option under Analyze that deals with scenes.
%
% Inputs:
%    hObject   - handle to menuAnalyzeScene (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Opens up a Scene Window
ind = vcGetSelectedObject('display');
if isempty(ind), disp('no display selected'); return; end
d = vcGetObject('display', ind);

% Get current image
global vcSESSION;
if isfield(vcSESSION, 'imgData')
    I = vcSESSION.imgData;
else
    warning('No image set');
    return;
end

% Get scene info
name = inputdlg({'Scene name'}, 'Scene Info', 1, ...
    {sprintf('scene-%s', displayGet(d, 'name'))});
if isempty(name), return; end

% Generate scene
if isempty(I), warning('No image set'); return; end
scene = sceneFromFile(I, 'rgb', [], d, [], 0);
scene = sceneSet(scene, 'name', name{1});

% Compute scene size
[r, c, ~] = size(I);
vDist = sceneGet(scene, 'distance');
fov = atand(max(r, c) * displayGet(d, 'metersperdot') / vDist);
scene = sceneSet(scene, 'fov', fov);

% Show scene window
vcAddAndSelectObject('scene', scene);
sceneWindow;

% --------------------------------------------------------------------
function menuClose_Callback(~, ~, handles)
% (File | Close) Close function call
%
% Syntax:
%   menuClose_Callback(~, ~, handles)
%
% Description:
%    Menu option to close the window
%
% Inputs:
%    hObject   - handle to menuClose (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
displayClose(handles.figure1);

% --------------------------------------------------------------------
function menuSaveDisplayModel_Callback(~, ~, ~)
% (File | Save Display) Save the current display
%
% Syntax:
%   menuSaveDisplayModel_Callback(~, ~, ~)
%
% Description:
%    Menu option under File to save the current display model.
%
% Inputs:
%    hObject   - handle to menuSaveDisplayModel (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ind = vcGetSelectedObject('display');
if isempty(ind), disp('no display selected'); return; end
d = vcGetObject('display', ind); %#ok
fname = uiputfile('*.mat');
if fname == 0, return; end
save(fname, 'd');

% --------------------------------------------------------------------
function displayRefresh(~, ~, handles)
% (Refresh) Refresh current display
%
% Syntax:
%   displayRefresh(~, ~, handles)
%
% Description:
%    Refresh the current display
%
% Inputs:
%    hObject   - handle to displayRefresh (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
displaySetEditsAndButtons(handles);

% --------------------------------------------------------------------
function figure1_CloseRequestFcn(~, ~, handles)
% (Close) Function request to close figure1
%
% Syntax:
%   figure1_CloseRequestFcn(~, ~, handles)
%
% Description:
%    Close request for figure1
%
% Inputs:
%    hObject   - handle to figure1 (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
displayClose(handles.figure1);

% --------------------------------------------------------------------
function menuCRT_Callback(hObject, eventdata, handles)
% (Display | CRT) Display CRT Menu option
%
% Syntax:
%   menuCRT_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option under Display to change the display to a CRT type.
%
% Inputs:
%    hObject   - handle to menuCRT (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
d = displayCreate('CRT-Dell');
vcAddAndSelectObject('display', d);
displayRefresh(hObject, eventdata, handles);

% --------------------------------------------------------------------
function menuShowInNewWindow_Callback(~, ~, handles)
% (Edit | Show in Zoomed View) New window
%
% Syntax:
%   menuShowInNewWindow_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option under Edit to open in a new window.
%
% Inputs:
%    hObject   - handle to menuShowInNewWindow (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
pos = get(handles.figure1, 'Position');
figure('Position', pos);
% Get current image
h = get(handles.axes1, 'children');
I = get(h, 'CData');
% Show image
imshow(I);

% --------------------------------------------------------------------
function menuLoadDisplayModel_Callback(hObject, eventdata, handles)
% (File | Load Display) Load a saved display
%
% Syntax:
%   menuLoadDisplayModel_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu option under File to load a saved display model.
%
% Inputs:
%    hObject   - handle to menuLoadDisplayModel (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
fname = uigetfile('*.mat', 'Load Display From MAT file');
if fname == 0, return; end
tmp = load(fname);
if isfield(tmp, 'd'), d = tmp.d; else, error('Not display file'); end
vcAddAndSelectObject('display', d);
displayRefresh(hObject, eventdata, handles);

% --------------------------------------------------------------------
function menuPlotDisplaySPD_Callback(~, ~, ~)
% (Plot | SPD Primaries) Plot the SPD Primaries
%
% Syntax:
%   menuPlotDisplaySPD_Callback(~, ~, ~)
%
% Description:
%    Menu option under Plot to display the SPD primaries.
%
% Inputs:
%    hObject   - handle to menuPlotDisplaySPD (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ind = vcGetSelectedObject('display');
if isempty(ind), disp('no display selected'); return; end
d = vcGetObject('display', ind);
displayPlot(d, 'spd');

% --------------------------------------------------------------------
function menuPlotGamut_Callback(~, ~, ~)
% (Plot | Gamut) Plot Gamut
%
% Syntax:
%   menuPlotGamut_Callback(~, ~, ~)
%
% Description:
%    Plot gamut display.
%
% Inputs:
%    hObject   - handle to menuPlotGamut (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ind = vcGetSelectedObject('display');
if isempty(ind), disp('no display selected'); return; end
d = vcGetObject('display', ind);
displayPlot(d, 'gamut');

% --------------------------------------------------------------------
function menuPlotGamma_Callback(~, ~, ~)
% (Plot | Gamma) Gamma plot
%
% Syntax:
%   menuPlotGamma_Callback(~, ~, ~)
%
% Description:
%    Menu under Plot to plot the gamma.
%
% Inputs:
%    hObject   - handle to menuPlotGamma (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ind = vcGetSelectedObject('display');
if isempty(ind), disp('no display selected'); return; end
d = vcGetObject('display', ind);
displayPlot(d, 'gamma');

% --------------------------------------------------------------------
function menuLCDHorizontalStripesRGB_Callback(~, ~, ~)
% NYI- (Display | LCD | Horizontal Stripes RGB) RGB Horizontal Stripes call
%
% Syntax:
%   menuLCDHorizontalStripesRGB_Callback(~, ~, ~)
%
% Description:
%    Display RGB Horizontal Stripes on an LCD. Not yet implemented.
%
% Inputs:
%    hObject   - handle to menuLCDHorizontalStripesRGB (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% displayGD = ctGetObject('display');
% dpi = 72;
% dSpacing = 0.001; % sample spacing in mm
% vDisplayLCD = vDisplayCreate('lcd', dpi, dSpacing, 'h', 'rgb');
%
% % Add? Replace? for displayGD
% displayGD = ctDisplaySet(displayGD, 'vDisplay', vDisplayLCD);
% ctSetObject('display', displayGD);
% ctdpRefreshGUIWindow(hObject);
warning('NYI');

% --------------------------------------------------------------------
function menuLCDVerticalStripesRGB_Callback(hObject, eventdata, handles)
% (Display | LCD | Vertical Stripes RGB) Display LCD's RGB Vertical stripes
%
% Syntax:
%   menuLCDVerticalStripesRGB_Callback(hObject, eventdata, handles)
%
% Description:
%    Display RGB vertical stripes on an LCD.
%
% Inputs:
%    hObject   - handle to menuLCDVerticalStripesRGB (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
d = displayCreate('LCD-Dell');
vcAddAndSelectObject('display', d);
displayRefresh(hObject, eventdata, handles);

% --------------------------------------------------------------------
function menuLCDVerticalStripesBGR_Callback(~, ~, ~)
% NYI- (Display | LCD | Vertical BGR) Display LCD vertical BGR stripes.
%
% Syntax:
%   menuLCDVerticalStripesBGR_Callback(~, ~, ~)
%
% Description:
%    Display vertical bgr stripes on an LCD
%
% Inputs:
%    hObject   - handle to menuLCDVerticalStripesBGR (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% displayGD = ctGetObject('display');
%
% dpi = 72;
% dSpacing = 0.001; % sample spacing in mm
% vDisplayLCD = vDisplayCreate('lcd', dpi, dSpacing, 'v', 'bgr');
%
% Add? Replace?
% displayGD = ctDisplaySet(displayGD, 'vDisplay', vDisplayLCD);
% ctSetObject('display', displayGD);
% ctdpRefreshGUIWindow(hObject);
warning('NYI');

% --------------------------------------------------------------------
function menuLCDHorizontalStripesBGR_Callback(~, ~, ~)
% NYI- (Display | LCD | Horizontal BGR) Horizontal BGR stripes on LCD.
%
% Syntax:
%   menuLCDHorizontalStripesBGR_Callback(~, ~, ~)
%
% Description:
%    Not yet implemented. Display Horizontal BGR stripes on an LCD.
%
% Inputs:
%    hObject   - handle to menuLCDHorizontalStripesBGR (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% displayGD = ctGetObject('display');
%
% dpi = 72;
% dSpacing = 0.001; % sample spacing in mm
% vDisplayLCD = vDisplayCreate('lcd', dpi, dSpacing, 'h', 'bgr');
%
% Add? Replace?
% displayGD = ctDisplaySet(displayGD, 'vDisplay', vDisplayLCD);
% ctSetObject('display', displayGD);
% ctdpRefreshGUIWindow(hObject);
warning('NYI');

% --------------------------------------------------------------------
function menuEditChangeFontSize_Callback(~, ~, handles)
% (Edit | Font Size) Change the font size
%
% Syntax:
%   menuEditChangeFontSize_Callback(~, ~, handles)
%
% Description:
%    Under the Edit menu, select the Font Size option, to change the font
%    size in the display.
%
% Inputs:
%    hObject   - handle to menuEditChangeFontSize (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% ctFontChangeSize(handles.figure1);
ieFontSizeSet(handles.figure1);

% answer = inputdlg('New Font Size (7~25)');
% if isempty(answer), return; end
% answer = str2double(answer);
% assert(answer > 6 && answer < 26, 'Front size out of range');
% set(handles.text61, 'FontSize', answer);
% set(handles.uipanelSummary, 'FontSize', answer);
% set(handles.txtSummary, 'FontSize', answer);
% set(handles.txtMaxLum, 'FontSize', answer);
% set(handles.txtPosVar, 'FontSize', answer);
% set(handles.txtAmp, 'FontSize', answer);
% set(handles.uipanel2, 'FontSize', answer);

% --- Executes on selection change in popupSelectDisplay.
function popupSelectDisplay_Callback(hObject, eventdata, handles)
% (Popup | Select Display) Select a display model
%
% Syntax:
%   popupSelectDisplay_Callback(hObject, eventdata, handles)
%
% Description:
%    Popup menu to select a display model.
%
% Inputs:
%    hObject   - handle to popupSelectDisplay (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Called when the 'Selected Display' popup is chosen
val = get(handles.popupSelectDisplay, 'value');
vcSetSelectedObject('display', val);
displayRefresh(hObject, eventdata, handles);

% --------------------------------------------------------------------
function menuEditRenameDisplay_Callback(hObject, eventdata, handles)
% (Edit | Rename) Rename the display model
%
% Syntax:
%   menuEditRenameDisplay_Callback(hObject, eventdata, handles)
%
% Description:
%    Menu under Edit to rename the display model.
%
% Inputs:
%    hObject   - handle to menuEditRenameDisplay (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ind = vcGetSelectedObject('display');
if isempty(ind), disp('no display selected'); return; end
answer = inputdlg('New Display Name');
if isempty(answer), return; end
d = vcGetObject('display', ind);
d = displaySet(d, 'name', answer{1});
vcDeleteSelectedObject('display');
vcAddAndSelectObject('display', d);
displayRefresh(hObject, eventdata, handles);

%-----------------------------------------------------
function editMaxLum_Callback(hObject, eventdata, handles)
% (Edit | Text box) Edit text box for maximum luminance.
%
% Syntax:
%   editMaxLum_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit the text box for Maximum Luminance.
%
% Inputs:
%    hObject   - handle to editMaxLum (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Edit box for max luminance
ind = vcGetSelectedObject('display');
if isempty(ind), disp('no display selected'); return; end
d = vcGetObject('display', ind);
xyz = displayGet(d, 'white xyz');
newLum = str2double(get(handles.editMaxLum, 'String'));
spd = displayGet(d, 'spd') * newLum / xyz(2);
d = displaySet(d, 'spd', spd);
vcDeleteSelectedObject('display');
vcAddAndSelectObject('display', d);
displayRefresh(hObject, eventdata, handles);

%-----------------------------------------------------
function editMaxLum_CreateFcn(hObject, ~, ~)
% (Create | Textbox) Create the textbox for Max Lum according to defaults.
%
% Syntax:
%   editMaxLum_CreateFcn(hObject, ~, ~)
%
% Description:
%    Create the textbox to control the Maximum Luminance.
%
% Inputs:
%    hObject   - handle to editMaxLum (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
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

function editVar_Callback(~, ~, ~)
% NYI- (Edit | Var) Call to modify var
%
% Syntax:
%   editVar_Callback(~, ~, ~)
%
% Description:
%    Not yet implemented. Call to edit the Var.
%
% Inputs:
%    hObject   - handle to editVar (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
disp('Not yet implemented')

function editVar_CreateFcn(hObject, ~, ~)
% (Create | Textbox) Create textbox for Var
%
% Syntax:
%   editVar_CreateFcn(hObject, ~, ~)
%
% Description:
%    Create function for Var text box
%
% Inputs:
%    hObject - handle to editVar (see GCBO)
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
%    * Hint: edit controls usually have a white background on Windows. See
%      ISPC and COMPUTER.
%

if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

function editPPI_Callback(hObject, eventdata, handles)
% (Edit | Textbox) Edit value for Display PPI.
%
% Syntax:
%   editPPI_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit value for the Display PPI textbox.
%
% Inputs:
%    hObject   - handle to editPPI (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ind = vcGetSelectedObject('display');
if isempty(ind), disp('no display selected'); return; end
d = vcGetObject('display', ind);
newPPI = str2double(get(handles.editPPI, 'String'));
d = displaySet(d, 'dpi', newPPI);
vcDeleteSelectedObject('display');
vcAddAndSelectObject('display', d);
displayRefresh(hObject, eventdata, handles);

function editPPI_CreateFcn(hObject, ~, ~)
% (Create | Textbox) Create the textbox to control Display PPI
%
% Syntax:
%   editPPI_CreateFcn(hObject, ~, ~)
%
% Description:
%    Create function for the Display PPI textbox.
%
% Inputs:
%    hObject - handle to editPPI (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty - The handles are not created until all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Hint: edit controls usually have a white background on Windows. See
%      ISPC and COMPUTER.
%

if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

% --------------------------------------------------------------------
function menuOLED_Callback(hObject, eventdata, handles)
% (Display | OLED) Select the OLED type of display.
%
% Syntax:
%   menuOLED_Callback(hObject, eventdata, handles)
%
% Description:
%    Select an OLED type for a display model.
%
% Inputs:
%    hObject   - handle to menuOLED (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
d = displayCreate('OLED-Sony');
vcAddAndSelectObject('display', d);
displayRefresh(hObject, eventdata, handles);

% --------------------------------------------------------------------
function menuProductHelp_Callback(~, ~, ~)
% NYI- (Help | Product Help) Product Help menu item under Help
%
% Syntax:
%   menuProductHelp_Callback(~, ~, ~)
%
% Description:
%    Not yet implemented. Menu under Help for Product Help.
%
% Inputs:
%    hObject   - handle to menuProductHelp (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
warning('NYI');

% --------------------------------------------------------------------
function menuAbout_Callback(~, ~, ~)
% NYI- (Help | About) About menu call
%
% Syntax:
%   menuAbout_Callback(~, ~, ~)
%
% Description:
%    Not yet implemented. Menu under Help, the 'About' selection.
%
% Inputs:
%    hObject   - handle to menuAbout (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
warning('NYI')

% --------------------------------------------------------------------
function menuCopyCurrent_Callback(hObject, eventdata, handles)
% (Edit | Copy) Copy the current display model.
%
% Syntax:
%   menuCopyCurrent_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit menu item to copy the current display model.
%
% Inputs:
%    hObject   - handle to menuCopyCurrent (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ind = vcGetSelectedObject('display');
if isempty(ind), disp('no display selected'); return; end
d = vcGetObject('display', ind);
dname = displayGet(d, 'name');
d = displaySet(d, 'name', ['copy - ' dname]);
vcAddAndSelectObject('display', d);
displayRefresh(hObject, eventdata, handles);

% --------------------------------------------------------------------
function popupSelectDisplay_CreateFcn(hObject, ~, ~)
% (Create | popUp) Create function for the select display popup.
%
% Syntax:
%   popupSelectDispay_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Create the popUp dropdown for Select Display according to defaults.
%
% Inputs:
%    hObject   - handle to popupSelectDisplay (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - empty - the handles are not created until after all of the
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
if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

% --------------------------------------------------------------------
function SubpixelScale_Callback(hObject, eventdata, handles)
% (Popup | Subpixel) Call to the subpixel scale plot
%
% Syntax:
%   SubpixelScale_Callback(hObject, eventdata, handles)
%
% Description:
%    Call to the subpixel plot in the bottom right corner.
%
% Inputs:
%    hObject   - handle to SubpixelScale (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ind = vcGetSelectedObject('display');
if isempty(ind), disp('no display selected'); return; end
answer = inputdlg('New samples per pixel (10 ~ 30)');
if isempty(answer), return; end
answer = str2double(answer);
assert(answer > 9 && answer < 31, 'samples per pixel out of range');
d = vcGetObject('display', ind);
psfs = displayGet(d, 'psfs');
d = displaySet(d, 'psfs', imresize(psfs, [answer answer]));
vcDeleteSelectedObject('display');
vcAddAndSelectObject('display', d);
displayRefresh(hObject, eventdata, handles);

% --------------------------------------------------------------------
function menuPlotGamut3D_Callback(~, ~, ~)
% (Analyze | Gamut3D) Analyze menu selection for Gamut LAB
%
% Syntax:
%   editPPI_Callback(hObject, eventdata, handles)
%
% Description:
%    Callback function for the Gamut LAB analyze option.
%
% Inputs:
%    hObject   - handle to menuPlotGamut3D (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ind = vcGetSelectedObject('display');
if isempty(ind), disp('no display selected'); end
d = vcGetObject('display', ind);
displayPlot(d, 'gamut3d');

% --- Executes on button press in Create Scene
function pushbutton9_Callback(hObject, eventdata, handles)
% (Button | 9) Call to push on button 9
%
% Syntax:
%   pushbutton9_Callback(hObject, eventdata, handles)
%
% Description:
%    Callback function for button press on 9.
%
% Inputs:
%    hObject   - handle to pushbotton9 (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
menuAnalyzeSceneSubpixel_Callback(hObject, eventdata, handles);

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, ~, ~)
% (KeyPress | Figure1) Key press analysis for figure1
%
% Syntax:
%   figure1_KeyPressFcn(hObject, eventdata, handles)
%
% Description:
%    Key press analysis function for figure1.
%
% Inputs:
%    hObject   - handle to figure1 (see GCBO)
%    eventdata - Struct. Structure contained the fields below (see FIGURE)
%         Key: Name of the key pressed, in lower case.
%         Character: The character interpretation of the pressed key(s).
%         Modifier: Pressed modifier key name(s) (i.e. control, shift)
%    handles   - structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(~, ~, ~)
% (WindowKeyPress | Figure1) Window key presses for figure1.
%
% Syntax:
%   figure1_WindowKeyPress(hObject, eventdata, handles)
%
% Description:
%    Function to analyze window key presses for figure1.
%
% Inputs:
%    hObject   - handle to figure1 (see GCBO)
%    eventdata - Struct. Structure contained the fields below (see FIGURE)
%         Key: Name of the key pressed, in lower case.
%         Character: The character interpretation of the pressed key(s).
%         Modifier: Pressed modifier key name(s) (i.e. control, shift)
%    handles   - structure with handles and user data (see GUIDATA)
%
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
