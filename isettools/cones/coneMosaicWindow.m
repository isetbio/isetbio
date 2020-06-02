function varargout = coneMosaicWindow(varargin)
% Cone image coneMosaicWindow interface
%
% Syntax:
%	varargout = coneMosaicWindow(varargin)
%
% Description:
%	 This is the CONEMOSAICWINDOW M-file for coneMosaicWindow.fig
%
%    Graphical user interface to manage the various Image Sensor Array
%    (ISA) properties.
%
%    CONEMOSAICWINDOW, by itself, creates a new CONEMOSAICWINDOW or raises
%    the existing singleton*.
%
%    H = CONEMOSAICWINDOW returns the handle to a new CONEMOSAICWINDOW or
%    the handle to the existing singleton*.
%
%    CONEMOSAICWINDOW('CALLBACK', hObject, eventData, handles, ...) calls
%    the local function named CALLBACK in CONEMOSAICWINDOW.M with the given
%    input arguments.
%
%    CONEMOSAICWINDOW('Property', 'Value', ...) creates a new
%    CONEMOSAICWINDOW or raises the existing singleton*. Starting from the
%    left, property value pairs are applied to the GUI before
%    sensorImageWindow_OpeningFunction gets called. An unrecognized
%    property name or invalid value makes property application stop. All
%    inputs are passed to coneMosaicWindow_OpeningFcn via varargin.
%
%  * See GUI Options on GUIDE's Tools menu. Choose "GUI allows only one
%    instance to run (singleton)".
%
% Inputs:
%    varargin  - (Optional) Input argument(s). Default is empty. Can
%                contain such information as existing cone mosaics.
%
% Outputs:
%    varargout - ???
%
% Optional key/value pairs:
%    ???
%

% History:
%    xx/xx/05         Copyright ImagEval Consultants, LLC, 2005.
%    06/08/17  GUIDE  Last Modified by GUIDE v2.5 08-Jun-2017 21:21:47
%    02/05/18  jnm    Formatting
%    04/07/17  dhb    Get rid of examples in subfunctions, they can't work.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton', gui_Singleton, ...
    'gui_OpeningFcn', @coneMosaicWindow_OpeningFcn, ...
    'gui_OutputFcn', @coneMosaicWindow_OutputFcn, ...
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
end

% --- Executes just before coneMosaicWindow is made visible.
function coneMosaicWindow_OpeningFcn(hObject, eventdata, handles, varargin)
% (Window | Open) Initialize and open the Cone Mosaic Window
%
% Syntax:
%   coneMosaicWindow_OpeningFcn(hObject, eventdata, handles, [varargin])
%
% Description:
%    This function will initialize and open the cone mosaic window using
%    the provided information
%
% Inputs:
%    hObject   - The handle to menuPlotCone (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%    varargin  - (Optional) Additional information for initializing the
%                window, such as an existing cone mosaic.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

%#ok<*DEFNU>
%#ok<*INUSD>
%#ok<*ST2NM>

% check inputs
if isempty(varargin) || ~isa(varargin{1}, 'coneMosaic')
    error('cone mosaic object required');
end
plotType = 'meanabsorptions';
if length(varargin) > 1, plotType = varargin{2}; end

% Choose default command line output for coneMosaicWindow
handles.output = hObject;
handles.cMosaic = varargin{1};
handles.mov = [];  % absorption movie
handles.curMov = [];  % photocurrent movie

% Update handles structure
guidata(hObject, handles);
handles.cMosaic.hdl = hObject;

% Adjust the database and bring this figure to the front
vcSetFigureHandles('conemosaic', hObject, eventdata, handles);
figure(hObject);

% Set the popup default image selection to mean absorptions for when the
% window opens.
str = get(handles.popupImageType, 'String');
for ii=1:numel(str), str{ii} = ieParamFormat(str{ii}); end

% Refresh and move on
coneMosaicGUIRefresh(hObject, eventdata, handles, plotType);

% Set the font size based on the ISETBIO preferences
ieFontInit(hObject);

% Very important for good rendering speed
set(hObject, 'Renderer', 'OpenGL')

end

% --- Outputs from this function are returned to the command line.
function varargout = coneMosaicWindow_OutputFcn(~, ~, handles)
% (Cone | Output) Outputs returned to the command line
%
% Syntax:
%   vargout = coneMosaicWindow_OutputFcn(hObject, eventdata, handles)
%
% Description:
%    Return the outputs from the window to the Command Line
%
% Inputs:
%    hObject   - The handle to coneMosaicWindow(see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    varargout - Output of varying length containing return text from the
%                included handles.
%
% Optional key/value pairs:
%    None.
%
varargout{1} = handles.output;
end

function btnComputeImage_Callback(hObject, eventdata, handles)
% (Button | Absorptions) Compute image from optics data - button press
%
% Syntax:
%   btnComputeImage_Callback(hObject, eventdata, handles)
%
% Description:
%    Function called at press of the button at the bottom of the window,
%    which will compute the image from the provided optics data.
%
% Inputs:
%    hObject   - The handle to menuPlotCone (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OI');
if isempty(oi) || isempty(oiGet(oi, 'photons'))
    warning('No optical image. Use ieAddObject(oi) to store.');
    return;
end

handles.cMosaic.compute(oi);
handles.cMosaic.name = oiGet(oi, 'name');
set(handles.popupImageType, 'Value', 2); % mean absorptions

coneMosaicGUIRefresh(hObject, eventdata, handles);

end

function menuAnComputeFromOI_Callback(hObject, eventdata, handles)
% (Compute | Absorptions) Compute absorptions from OI in database
%
% Syntax:
%   menuAnComputeFromOI_Callback(hObject, eventdata, handles)
%
% Description:
%    Press the button in the window, and compute the cone absorptions from
%    an OI in the database.
%
% Inputs:
%    hObject   - The handle to menuAnComputeFromOI (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
btnComputeImage_Callback(hObject, eventdata, handles);
end

function editRows_Callback(hObject, eventdata, handles)
% (Edit | Textbox) Change the value of the Rows text box
%
% Syntax:
%   editRows_Callback(hObject, eventdata, handles)
%
% Description:
%    Change the value of the Rows text box.
%
% Inputs:
%    hObject   - The handle to editRows (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.rows = str2double(get(hObject, 'String'));
menuEditClearData_Callback(hObject, eventdata, handles)
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

function editCols_Callback(hObject, eventdata, handles)
% (Edit | Textbox) Change the value of the Columns text box
%
% Syntax:
%   editCols_Callback(hObject, eventdata, handles)
%
% Description:
%    Change the value of the Columns text box.
%
% Inputs:
%    hObject   - The handle to editCols (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.cols = str2double(get(hObject, 'String'));
menuEditClearData_Callback(hObject, eventdata, handles)
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

function editExpTime_Callback(hObject, eventdata, handles)
% (Edit | Textbox) Change the value of the integration time text box
%
% Syntax:
%   editExpTime_Callback(hObject, eventdata, handles)
%
% Description:
%    Change the value of the Integration Time text box.
%
% Inputs:
%    hObject   - The handle to editExpTime (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.integrationTime = 1e-3 * ...
    str2double(get(hObject, 'String'));
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% GUI object create functions
function editRows_CreateFcn(hObject, eventdata, handles)
% (Create | Textbox) Initialize the Rows text box
%
% Syntax:
%   menuPlotCone_Callback(hObject, eventdata, handles)
%
% Description:
%    Initialize the Rows text box, unless specified, set the background
%    color to default.
%
% Inputs:
%    hObject   - The handle to editRows (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
set(hObject, 'BackgroundColor', get(0, 'defaultUicontrolBackgroundColor'));
end

function editCols_CreateFcn(hObject, eventdata, handles)
% (Create| Textbox) Initialize the Columns text box
%
% Syntax:
%   editCols_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Initialize the Columns text box, unless specified, set the background
%    color to default.
%
% Inputs:
%    hObject   - The handle to editCols (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
set(hObject, 'BackgroundColor', get(0, 'defaultUicontrolBackgroundColor'));
end

function editExpTime_CreateFcn(hObject, eventdata, handles)
% (Create| Textbox) Initialize the Integration Time text box
%
% Syntax:
%   editExpTime_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Initialize the Integration Time text box, and unless specified, set
%    the background color to default.
%
% Inputs:
%    hObject   - The handle to editExpTime (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles are not created until all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
set(hObject, 'BackgroundColor', get(0, 'defaultUicontrolBackgroundColor'));
end

% menu call back functions
function menuFile_Callback(hObject, eventdata, handles)
% (Menu | File) Top level File object
%
% Syntax:
%   menuFile_Callback(hObject, eventdata, handles)
%
% Description:
%    Top level of File menu object
%
% Inputs:
%    hObject   - The handle to menuFile (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

function menuFileClose_Callback(~, ~, handles)
% (Menu | File) Close the window
%
% Syntax:
%   menuFileClose_Callback(hObject, eventdata, handles)
%
% Description:
%    Close the window (delete the handle to the window)
%
% Inputs:
%    hObject   - The handle to menuFileClose (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
delete(handles.coneMosaicWindow);
end

function menuEdit_Callback(hObject, eventdata, handles)
% (Menu | Edit) Top level Edit menu
%
% Syntax:
%   menuEdit_Callback(hObject, eventdata, handles)
%
% Description:
%    Top level edit menu
%
% Inputs:
%    hObject   - The handle to menuEdit (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

function menuEditName_Callback(hObject, eventdata, handles)
% (Menu | Edit | Rename) Rename an existing object
%
% Syntax:
%   menuEditName_Callback(hObject, eventdata, handles)
%
% Description:
%    Rename an existing object
%
% Inputs:
%    hObject   - The handle to menuPlotCone (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
str = ieReadString('New name', handles.cMosaic.name);
if ~isempty(str), handles.cMosaic.name = str; end
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

function menuEditClearData_Callback(hObject, eventdata, handles)
% (Menu | Edit | Clear Data) Clear the data
%
% Syntax:
%   menuEditClearData_Callback(hObject, eventdata, handles)
%
% Description:
%    Clear the data
%
% Inputs:
%    hObject   - The handle to menuEditClearData (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.clearData();  % Clear absorptions and current

% mosaicImage is the colorful cone mosaic image.
% We used to save the absorption and current movie.  But no more.  It
% doesn't seem to save any time.
uData = get(handles.axes2,'UserData');
if isfield(uData,'mosaicImage')
    uData.mosaicImage = [];
    set(handles.axes2,'UserData',uData);
end

guidata(hObject, handles); % Put back the modified handles

coneMosaicGUIRefresh(hObject, eventdata, handles); % Update the window
end

% --- Executes during object creation, after setting all properties.
function editGam_CreateFcn(hObject, eventdata, handles)
% (Edit | Create) Text book Gamm Initialization
%
% Syntax:
%   editGam_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Text book Gamm creation
%
% Inputs:
%    hObject   - The handle to editGam (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
set(hObject, 'BackgroundColor', get(0, 'defaultUicontrolBackgroundColor'));
end

function editGam_Callback(hObject, eventdata, handles)
% (Edit | ) Text book Gamm
%
% Syntax:
%   menuPlotCone_Callback(hObject, eventdata, handles)
%
% Description:
%    "Text book Gamm" - edit
%
% Inputs:
%    hObject   - The handle to editGam (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

set(handles.editGam, 'value', str2double(get(handles.editGam, 'string')));
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

function coneMosaicGUIRefresh(~, ~, handles, plotType)
% (File | Refresh) Update the Cone mosaic window interface - main pulldown
%
% Syntax:
%   coneMosaicGUIRefresh(handles)
%
% Description:
%    Update the cone mosaic window's GU interface, via the main pulldown.
%
% Inputs:
%    hObject   - The handle to menuPlotCone (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/16  HJ, BW  ISETBIO TEAM, 2016
%    02/08/18  jnm     Formatting

% get coneMosaic object
cm = handles.cMosaic;

% Place name in text string box
set(handles.mosaicName, 'string', sprintf('%s', cm.name));

% set row and cols
set(handles.editRows, 'string', num2str(cm.rows));
set(handles.editCols, 'string', num2str(cm.cols));

% set integration time
set(handles.editExpTime, 'string', ...
    sprintf('%.1f', cm.integrationTime * 1e3));

% set KLMS ratio
str = sprintf('[%.1f, %.1f, %.1f, %.1f]', cm.spatialDensity(1), ...
    cm.spatialDensity(2), cm.spatialDensity(3), cm.spatialDensity(4));
set(handles.editKLMS, 'string', str);

% set description strings
str = cm.description('skipMacular', true, 'skipPigment', true);
set(handles.txtMosaic, 'string', str);

% set outersegment description
str = cm.descriptionOS;
set(handles.txtOS, 'string', str);

% set photopigment properties
set(handles.editConeWidth, 'string', num2str(cm.pigment.width * 1e6));
set(handles.editConeHeight, 'string', num2str(cm.pigment.height * 1e6));

str = sprintf('[%.1f, %.1f, %.1f]', cm.pigment.opticalDensity(1), ...
    cm.pigment.opticalDensity(2), cm.pigment.opticalDensity(3));
set(handles.editConeOpticalDensity, 'string', str);

str = sprintf('[%.2f, %.2f, %.2f]', cm.pigment.peakEfficiency(1), ...
    cm.pigment.peakEfficiency(2), cm.pigment.peakEfficiency(3));
set(handles.editConePeakEfficiency, 'string', str);

% set macular density
set(handles.editMacularDensity, 'string', num2str(cm.macular.density));

% Set eccentricity in the window based on the center. This is specified in
% meters, and we convert it to deg for the window.
ecc = sqrt(sum(cm.center .^ 2));   % Meters

% Why don't we have a builtin variable for this?  Or a way to compute it?
% I think we do ... help! (BW)
deg2m = 3.3333e-04;
ecc = ecc / deg2m;
set(handles.editEccentricity, 'string', num2str(ecc, 2));

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

%% Set the strings in the popup menu for what to plot
str = {'Cone mosaic'};

% Add absorption options
if ~isempty(cm.absorptions)
    str = [str {'Mean absorptions', 'Absorption movie'}];
end

% Add photocurrent options
if ~isempty(cm.current)
    str = [str {'Mean photocurrent', 'Photocurrent movie'}];
end

if ~exist('plotType','var')
    % The user did not specify the plot type.  So we get it from the
    % window.
    index = get(handles.popupImageType, 'Value');
    if index > length(str), index = 1; end  % Cone mosaic
    plotType = str{index};
else
    plotType = ieParamFormat(plotType);
    tmp = cell(size(str));
    for ii=1:numel(str), tmp{ii} = ieParamFormat(str{ii}); end
    % The user did specify plotType.  Find it and use it
    [~,index] = ismember(plotType,tmp);
end

% Set the index and then string in the window popup menu
set(handles.popupImageType, 'Value', index);
set(handles.popupImageType, 'String', str);

%% Here are the different plotting options for the main axis (handles.axes2)

% Gamma is managed within plot for the mean images. For the video, we are
% handling it this way so we can, in the future, set additional parameters
% for the movies that get passed to ieMovie.
g = str2num(get(handles.editGam, 'string'));

switch ieParamFormat(plotType)
    case 'conemosaic'
        % cone mosaic image
        % TODO:  For large mosaics, the computation is slow. We should
        % compute it once and store it.
        nCones  = size(cm.coneLocs, 1);
        maxCones = 1e4;
        if  nCones > maxCones
            ieInWindowMessage('Initializing large mosaic image', handles);
        end
        cm.plot('cone mosaic', 'hf', handles.axes2);
        ieInWindowMessage('', handles)

    case 'meanabsorptions'
        % mean cone absorptions
        cm.plot('mean absorptions', 'hf', handles.axes2);
        
    case 'absorptionmovie'
        ieInWindowMessage('Showing absorption movie', handles)
        cm.plot('movie absorptions', 'hf', handles.axes2, ...
            'gamma', g);
        % Saved the movie data in the axis as userdata.data
        ieInWindowMessage('', handles)
        
    case 'meanphotocurrent'
        cm.plot('mean current', 'hf', handles.axes2);

    case 'photocurrentmovie'
        ieInWindowMessage('Showing photocurrent movie', handles)
        cm.plot('movie current', 'hf', handles.axes2, 'gamma', g);
        ieInWindowMessage('', handles)
        
    otherwise
        error('Unknown plot type');
end

enable.hLine = 'on';
enable.vLine = 'on';
enable.hLineLMS = 'on';
enable.vLineLMS = 'on';
enable.timeSeries = 'on';
contextMenuInit(handles);
contextMenuEnable(handles, enable);

end

function c = contextMenuInit(handles)
% Set up the context (right click) menu.
%
% Syntax:
%   c = contextMenuInit(handles)
%
% Description:
%    Set up the context (right click) menu.
%
% Inputs:
%    handles - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    c       - The created context menu
%
% Optional key/value pairs:
%    None.
%

c = uicontextmenu;
if ~isempty(handles.axes2.Children)
    for ichild = 1:size(handles.axes2.Children, 1)
        handles.axes2.Children(ichild).UIContextMenu = c;
    end
    uimenu(c, 'Label', 'hLine response', 'Callback', @contextMenuPlot);
    uimenu(c, 'Label', 'vLine response', 'Callback', @contextMenuPlot);
    uimenu(c, 'Label', 'hLine LMS', 'Callback', @contextMenuPlot);
    uimenu(c, 'Label', 'vLine LMS', 'Callback', @contextMenuPlot);
    uimenu(c, 'Label', 'time series', 'Callback', @contextMenuPlot);
end

end

function contextMenuEnable(handles, enable)
% Set and/or enable plot options in menu
%
% Syntax:
%   menuPlotCone_Callback(hObject, eventdata, handles)
%
% Description:
%    Set the enable plot options in the menu. Enable is a structure
%    containing 'on' and 'off' values.
%
% Inputs:
%    handles - Structure with handles and user data (see GUIDATA)
%    enable  - A structure containing 'on' and 'off' values in reference to
%              the menu options.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% enable is a structure of 'on' and 'off' values
set(handles.menuPlotHLine, 'Enable', enable.hLine);
set(handles.menuPlotVLine, 'Enable', enable.vLine);
set(handles.menuPlotHLineLMS, 'Enable', enable.hLineLMS);
set(handles.menuPlotVLineLMS, 'Enable', enable.vLineLMS);
set(handles.menuPlotTimeSeries, 'Enable', enable.timeSeries);

end

function contextMenuPlot(source, callbackdata)
% Callback function for five of the context menu plots
%
% Syntax:
%   contextMenuPlot(source, callbackdata)
%
% Description:
%    Callback function for the following five context menu plots:
%       hline, vline, hLineLMS, vLineLMS, timeSeries
%
%    Data may cover absorptions or current for any of the five
%    aforementioned plots.
%
% Inputs:
%    hObject   - The handle to menuPlotCone (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    source       - The structure with handles and user data (see GUIDATA),
%                   also known as handles in other functions.
%    callbackdata - Data for callback
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% The guidata return of the variable 'source' contains the gui objects
handles = guidata(source);

% determine which data to use (absorption or current)
contents = get(handles.popupImageType, 'String');
index = get(handles.popupImageType, 'Value');
if index > length(contents), index = 1; end
plotType = contents{index};

% Identify the data type, absorptions or current
switch plotType
    case {'Mean absorptions', 'Cone mosaic', 'Absorption movie'}
        dataType = 'absorptions';
    case {'Mean photocurrent', 'Photocurrent movie'}      
        dataType = 'current';
end

% Figure out which plot was requested and build the command
switch ieParamFormat(source.Label)
    case 'hlineresponse'
        cmd = ['hline', dataType];
    case 'vlineresponse'
        cmd = ['vline', dataType];
    case 'hlinelms'
        cmd = ['hline', dataType, 'lms'];
    case 'vlinelms'
        cmd = ['vline', dataType, 'lms'];
    case 'timeseries'
        cmd = ['time series', dataType];
    otherwise
        error('Unknown plot type %s\n', source.label);
end

% Call the plot command, setting the main window axis for the first place
% to start.
handles.cMosaic.plot(cmd, 'hf', handles.axes2);

end

function menuPlot_Callback(hObject, eventdata, handles)
% (Menu | Plot) Top level Plot menu
%
% Syntax:
%   menuPlot_Callback(hObject, eventdata, handles)
%
% Description:
%    Top level of Plot
%
% Inputs:
%    hObject   - The handle to menuPlot (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

function menuEditFontSize_Callback(~, ~, handles)
% (Menu | Edit | Font Size) Modify the font size
%
% Syntax:
%   menuEditFontSize_Callback(hObject, eventdata, handles)
%
% Description:
%    Modify the font size
%
% Inputs:
%    hObject   - The handle to menuEditFontSize (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
ieFontSizeSet(handles.coneMosaicWindow);
end

function menuHelp_Callback(hObject, eventdata, handles)
% (Menu | Help) Help Menu
%
% Syntax:
%   menuHelp_Callback(hObject, eventdata, handles)
%
% Description:
%    Top level of Help menu
%
% Inputs:
%    hObject   - The handle to menuHelp (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

function menuAppNotes_Callback(hObject, eventdata, handles)
% (Menu | Help | Web Documentation) Call web documentation
%
% Syntax:
%   menuAppNotes_Callback(hObject, eventdata, handles)
%
% Description:
%    Open the web documentation
%
% Inputs:
%    hObject   - The handle to menuPlotCone (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
web('https://github.com/isetbio/isetbio/wiki', '-browser');
end

function editKLMS_Callback(hObject, eventdata, handles)
% (Edit | KLMS) Modify the KLMS
%
% Syntax:
%   editKLMS_Callback(hObject, eventdata, handles)
%
% Description:
%    Modify the KLMS
%
% Inputs:
%    hObject   - The handle to menuPlotCone (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
density = str2num(get(handles.editKLMS, 'String'));
assert(numel(density) == 4, 'invalid input');

density = density / sum(density);
handles.cMosaic.spatialDensity = density;
menuEditClearData_Callback(hObject, eventdata, handles);
end

function editKLMS_CreateFcn(hObject, eventdata, handles)
% (Edit | KLMS) Initialize the KLMS
%
% Syntax:
%   editKLMS_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Initialize the KLMS
%
% Inputs:
%    hObject   - The handle to editKLMS (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles are not created until after all of the
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
if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

end

function txtMosaic_CreateFcn(hObject, eventdata, handles)
% (Mosaic | Text) Create the txtMosaic
%
% Syntax:
%   txtMosaic_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Create the txtMosaic object
%
% Inputs:
%    hObject   - The handle to txtMosaic (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

function editConeWidth_Callback(hObject, eventdata, handles)
% (Edit | Cone Width) Modify the width of an existing cone
%
% Syntax:
%   editConeWidth_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit the cone width
%
% Inputs:
%    hObject   - The handle to editConeWidth (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
newWidth = 1e-6 * str2double(get(hObject, 'String'));

if handles.cMosaic.pigment.width ~= newWidth
    handles.cMosaic.pigment.width = newWidth;
    menuEditClearData_Callback(hObject, eventdata, handles);
end

end

% --- Executes during object creation, after setting all properties.
function editConeWidth_CreateFcn(hObject, eventdata, handles)
% (Edit | Cone Width) Initialize the cone width and color default bg
%
% Syntax:
%   editConeWidth_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Initialize the cone width, and if not specified, color the background
%    to white.
%
% Inputs:
%    hObject   - The handle to editConeWidth (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles are not created until after all
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
if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end
end

function editConeHeight_Callback(hObject, eventdata, handles)
% (Edit | Cone Height) Modify the cone height
%
% Syntax:
%   editConeHeight_Callback(hObject, eventdata, handles)
%
% Description:
%    Modify the cone height.
%
% Inputs:
%    hObject   - The handle to editConeHeight (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
newHeight = 1e-6 * str2double(get(hObject, 'String'));
if handles.cMosaic.pigment.height ~= newHeight
    handles.cMosaic.pigment.height = newHeight;
    menuEditClearData_Callback(hObject, eventdata, handles);
end
end

% --- Executes during object creation, after setting all properties.
function editConeHeight_CreateFcn(hObject, eventdata, handles)
% (Edit | Cone Height) Create the cone height after properties set.
%
% Syntax:
%   editConeHeight_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    After the properties have been set, create the cone height, and set
%    default background color to white (unless specified otherwise).
%
% Inputs:
%    hObject   - The handle to editConeHeight (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles are not created until after all of the
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
end

function editConeOpticalDensity_Callback(hObject, eventdata, handles)
% (Edit | Cone Optical Density) Edit the cone optical density
%
% Syntax:
%   menuPlotConeOpticalDensity_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot the cone optical density.
%
% Inputs:
%    hObject   - The handle to coneOpticalDensity (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
val = str2num(get(hObject, 'String'));
assert(numel(val) == 3, 'invalid input for optical density');

if any(handles.cMosaic.pigment.opticalDensity(:) ~= val(:))
    handles.cMosaic.pigment.opticalDensity = val;
    menuEditClearData_Callback(hObject, eventdata, handles);
end
end

% --- Executes during object creation, after setting all properties.
function editConeOpticalDensity_CreateFcn(hObject, eventdata, handles)
% (Edit | Cone Optical Density) Create the Cone Optical Density
%
% Syntax:
%   editConeOpticalDensity_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    After all of the properties have been set, create the object
%    containing the cone optical density
%
% Inputs:
%    hObject   - The handle to editConeOpticalDensity (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles are not created until after all of the
%                CreateFcns are called.
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
if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end
end

function editMacularDensity_Callback(hObject, eventdata, handles)
% (Edit | Macular Density) Edit the existing macular density
%
% Syntax:
%   editMacularDensity_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit the existing macular density
%
% Inputs:
%    hObject   - The handle to editMacularDensity (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
val = str2double(get(hObject, 'String'));
if handles.cMosaic.macular.density ~= val
    handles.cMosaic.macular.density = val;
    menuEditClearData_Callback(hObject, eventdata, handles);
end
end

function editMacularDensity_CreateFcn(hObject, eventdata, handles)
% (Edit | Macular Density) Set default Macular Density background to white
%
% Syntax:
%   editMacularDensity_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    If following defaults, set the Macular Density background color to
%    white after all of the createFcns have completed.
%
% Inputs:
%    hObject   - The handle to editMacularDensity (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles are not created until after all of the
%                CreateFcns are called.
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
if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end
end

function editConePeakEfficiency_Callback(hObject, eventdata, handles)
% (Edit | Cone Peak Efficiency) Edit the Cone peak efficiency
%
% Syntax:
%   editConePeakEfficiency_Callback(hObject, eventdata, handles)
%
% Description:
%    Edit the cone peak efficiency
%
% Inputs:
%    hObject   - The handle to editConePeakEfficiency (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
val = str2num(get(hObject, 'String'));
assert(numel(val) == 3, 'invalid input for peak efficiency');
if any(handles.cMosaic.pigment.peakEfficiency(:) ~= val(:))
    handles.cMosaic.pigment.peakEfficiency = val;
    menuEditClearData_Callback(hObject, eventdata, handles);
end
end

function editConePeakEfficiency_CreateFcn(hObject, eventdata, handles)
% (Edit | Cone Peak Efficiency) Set default background to white
%
% Syntax:
%   editConePeakEfficiency_Callback(hObject, eventdata, handles)
%
% Description:
%    Initialize the cone peak efficiency - if default, set the background
%    color to white.
%
% Inputs:
%    hObject   - The handle to editConePeakEfficiency (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, handles are not created until after all of the
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
if ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end
end

function menuPlotMacular_Callback(hObject, eventdata, handles)
% (Edit | Plot | Macular) Plot Macular top level
%
% Syntax:
%   menuPlotMacular_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot the Macular.
%
% Inputs:
%    hObject   - The handle to menuPlotMacular (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

function menuPlotCone_Callback(hObject, eventdata, handles)
% (Menu | Plot | Cone) Plot the cone
%
% Syntax:
%   menuPlotCone_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot the Cone.
%
% Inputs:
%    hObject   - The handle to menuPlotCone (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

function menuPlotMosaic_Callback(hObject, eventdata, handles)
% (Menu | Mosaic) Mosaic Top Level menu
%
% Syntax:
%   menuPlotMosaic_Callback(hObject, eventdata, handles)
%
% Description:
%    Plot the mosaic
%
% Inputs:
%    hObject   - The handle to menuPlotMosaic (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

function menuPlotMosaicConeMosaic_Callback(~, ~, handles)
% (Menu | Plot | Mosaic ) Plot the Cone Mosaic from the Mosaic menu
%
% Syntax:
%   menuPlotMosaicConeMosaic_Callback(~, ~, handles)
%
% Description:
%    Plot the cone mosaic from the Mosaic menu
%
% Inputs:
%    hObject   - The handle to menuPlotMosaicConeMosaic (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.plot('cone mosaic', ...
    'showCorrespondingRectangularMosaicInstead', false);
end

function menuPlotMosaicMeanAbsorptions_Callback(~, ~, handles)
% (Menu | Plot | Mosaic) Plot the mean absorptions from Mosaic
%
% Syntax:
%   menuPlotMosaicMeanAbsorptions_Callback(~, ~, handles)
%
% Description:
%    Plot the mean absorptions from the Mosaic menu under Plot
%
% Inputs:
%    hObject   - The handle to menuPlotMosaicMeanAbsorptions (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.plot('mean absorptions');
end

function menuPlotMosaicMeanCurrent_Callback(~, ~, handles)
% (Menu | Plot | Mosaic) Plot Mean Current from Mosaic
%
% Syntax:
%   menuPlotMosaicMeanCurrent_Callback(~, ~, handles)
%
% Description:
%    Plot the mean current from the mosaic menu
%
% Inputs:
%    hObject   - The handle to menuPlotMosaicMeanCurrent (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.plot('mean current');
end

function menuPlotConeAbsorptance_Callback(~, ~, handles)
% (Menu | Plot | Cone Absorptance) Plot cone absorptance
%
% Syntax:
%   menuPlotConeAbsorptance_Callback(~, ~, handles)
%
% Description:
%    Plot the cone absorptance.
%
% Inputs:
%    hObject   - The handle to menuPlotConeAbsorptance (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.plot('cone fundamentals');
end

function menuPlotMacularTransmittance_Callback(~, ~, handles)
% (Menu | Plot | Macular) Plot the macular transmittance
%
% Syntax:
%   menuPlotMacularTransmittance_Callback(~, ~, handles)
%
% Description:
%    Plot the macular transmittance
%
% Inputs:
%    hObject   - The handle to menuPlotMacularTransmittance (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.plot('macular transmittance');
end

function menuPlotMacularAbsorptance_Callback(~, ~, handles)
% (Menu | Plot | Macular) Plot the macular absorptance
%
% Syntax:
%   menuPlotMacularAbsorptance_Callback(~, ~, handles)
%
% Description:
%    Plot the macular absorptance
%
% Inputs:
%    hObject   - The handle to menuPlotMacularAbsorptance (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.plot('macular absorptance');
end

function menuPlotMacularAbsorbance_Callback(~, ~, handles)
% (Menu | Plot | Macular) Plot the macular absorbance
%
% Syntax:
%   menuPlotEMPath_Callback(~, ~, handles)
%
% Description:
%    Plot the macular absorbance
%
% Inputs:
%    hObject   - The handle to menuPlotMacularAbsorbance (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.plot('macular absorbance');

end

function menuFileRefresh_Callback(hObject, eventdata, handles)
% (Menu | File | Refresh) Refresh button
%
% Syntax:
%   menuFileRefresh_Callback(~, ~, handles)
%
% Description:
%    Refresh button in file menu
%
% Inputs:
%    hObject   - The handle to menuFileRefresh (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

function menuPlotEMPath_Callback(~, ~, handles)
% (Menu | Plot | EMPath) Plot Eye Movement Path
%
% Syntax:
%   menuPlotEMPath_Callback(~, ~, handles)
%
% Description:
%    Plot the eye movement path.
%
% Inputs:
%    hObject   - The handle to menuPlotEMPath (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.plot('eye movement path');
end

% --------------------------------------------------------------------
function menuCones_Callback(hObject, eventdata, handles)
% (Menu | Cones) Main Pull down for cones menu.
%
% Syntax:
%   menuCones_Callback(hObject, eventdata, handles)
%
% Description:
%    The main pulldown function for computing cones.
%
% Inputs:
%    hObject   - The handle to menuCones (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
end

function menuConesGenerateEM_Callback(hObject, eventdata, handles)
% (Menu | Cones | Eye Movement) Generates eye movements
%
% Syntax:
%   menuConesGenerateEM_Callback(hObject, eventdata, handles)
%
% Description:
%    Generate the cone's eye movements
%
% Inputs:
%    hObject   - The handle to menuConesGenerateEM (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
str = ieReadString('Number of frames', '500');
if ~isempty(str)
    % Changed by BW on Nov 12, 2019.  The path was not randomized when
    % called from the pulldown.  This rSeed calls shuffle.
    % Must check with others about this!
    handles.cMosaic.emGenSequence(str2double(str),'rSeed',[]);
    menuEditClearData_Callback(hObject, eventdata, handles);
    set(handles.popupImageType, 'Value', 2); % mean absorptions
    coneMosaicGUIRefresh(hObject, eventdata, handles);
end

end

function menuConesAbsorptions_Callback(hObject, eventdata, handles)
% (Menu | Cones | Absorptions) Computes absorptions
%
% Syntax:
%   menuConesAbsorptions_Callback(hObject, eventdata, handles)
%
% Description:
%    Load the current oi to compute the absorptions. If no oi is selected,
%    the function will complain.
%
% Inputs:
%    hObject   - The handle to menuConesAbsorptions (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
oi = vcGetObject('OI');
if isempty(oi) || isempty(oiGet(oi, 'photons'))
    warning('No optical image. Use ieAddObject(oi) to store.');
    return;
end

fprintf('Calculating with optical image %s\n', oiGet(oi, 'name'));
handles.cMosaic.compute(oi);
handles.cMosaic.name = oiGet(oi, 'name');
set(handles.popupImageType, 'Value', 2); % mean absorptions
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

function menuConesPhotocurrent_Callback(hObject, eventdata, handles)
% (Menu | Cones | Photocurrent) Compute photocurrent
%
% Syntax:
%   menuConePhotocurrent_Callback(hObject, eventdata, handles)
%
% Description:
%    Compute the photocurrent
%
% Inputs:
%    hObject   - The handle to menuConePhotocurrent (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
handles.cMosaic.computeCurrent;
set(handles.popupImageType, 'Value', 4); % mean current
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

function menuConePhotocurrentNoise_Callback(hObject, eventdata, handles)
% (Menu | Cones | Photocurrent) Toggle photocurrent noise
%
% Syntax:
%   menuConePhotocurrentNoise_Callback(hObject, eventdata, handles)
%
% Description:
%    Toggle the photocurrent noise setting. The toggle command produces the
%    following changes in the noise setting:
%       Random -> None
%       Frozen -> Random
%       None   -> Frozen
%
% Inputs:
%    hObject   - The handle to menuConePhotocurrentNoise (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Also executes computeCurrent
%

% set(handles.btnPlayPause, 'Value', 0);  % Turn off any movie.

% Flip from whatever state to the other
switch handles.cMosaic.os.noiseFlag
    case 'random'
        handles.cMosaic.os.noiseFlag = 'none';
    case 'frozen'
        handles.cMosaic.os.noiseFlag = 'random';
    case 'none'
        handles.cMosaic.os.noiseFlag = 'frozen';
end

% We used to use this. But now, I think it should always be off.
handles.menuConePhotocurrentNoise.Checked = 'off';
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

%---------------
function popupImageType_Callback(hObject, eventdata, handles)
% (Image | Refresh Image) Refresh image object
%
% Syntax:
%   popupImageType_Callback(hObject, eventdata, handles)
%
% Description:
%    Refresh the image object
%
% Inputs:
%    hObject   - The handle to popupImageType (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Stop any movies.
% set(handles.btnPlayPause, 'Value', 0);  % Turn off the movie.

% Refresh.
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

function popupImageType_CreateFcn(hObject, eventdata, handles)
% (Image | Background Color) Set background color to white if default
%
% Syntax:
%   popupImageType_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    If background set to default, make the image background white.
%
% Inputs:
%    hObject   - The handle to popupImageType (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles not created until after all CreateFcns
%                have been called.
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
end

function sliderMovieProgress_Callback(~, ~, handles)
% (Movie | Set Progress Bar) Set slider to point in timeline of movie
%
% Syntax:
%   menuPlotVLine_Callback(hObject, eventdata, handles)
%
% Description:
%    Set the slider to a specific time point on the bar.
%
% Inputs:
%    hObject   - The handle to sliderMovieProgress (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
index = get(handles.popupImageType, 'Value');
if index == 3, mov = handles.cMosaic.absorptions;
elseif index == 5, mov = handles.cMosaic.current;
end

gam = get(handles.editGam, 'value');

mov = ieScale(mov, 0, 1) .^ gam;
mind = min(mov(:));
maxd = max(mov(:));

cnt = round(get(handles.sliderMovieProgress, 'Value'));
assert(cnt <= size(mov, ndims(mov)), 'slider choice out of range');
axes(handles.axes2);
imagesc(mov(:, :, cnt));
axis image;
set(gca, 'xticklabel', '', 'yticklabel', '');
caxis([mind maxd]);
drawnow;
set(handles.txtMovieFrame, 'string', cnt);

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
% (Create | Slider Background) Set default slider background to black
%
% Syntax:
%   sliderMovieProgress_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Set the slider background to black after all attributes and parameters
%    have been set.
%
% Inputs:
%    hObject   - The handle to sliderMovieProgress (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, handles are not created until after all of the
%                CreateFcns have been called.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
if isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', [.9 .9 .9]);
end
end

function editEccentricity_Callback(hObject, eventdata, handles)
% (Edit | ???) An editing event call
%
% Syntax:
%   menuPlotVLine_Callback(hObject, eventdata, handles)
%
% Description:
%    An empty editing event call, with "hints" below?
%
% Inputs:
%    hObject   - The handle to editEccentricity (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data (see GUIDATA)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Hints:
%   get(hObject, 'String')             returns contents of editEccentricity
%                                      as text
%   str2double(get(hObject, 'String')) returns contents of editEccentricity
%                                      as a double
end

% --- Executes during object creation, after setting all properties.
function editEccentricity_CreateFcn(hObject, eventdata, handles)
% (Create | Background) If using defaults, set background color to white.
%
% Syntax:
%   editEccentricity_CreateFcn(hObject, eventdata, handles)
%
% Description:
%    Function executes during object creation, after setting all of the
%    associated properties.
%
%    Function sets the background color to white if default values are to
%    be used.
%
% Inputs:
%    hObject   - The handle to editEccentricity (see GCBO)
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Empty, the handles are not created until after all
%                CreateFcns are called.
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

end

function menuPlotHLine_Callback(hObject, eventdata, handles)
% (Menu | Plot | H Line Response) Plot H Cone responses from a point
%
% Syntax:
%   menuPlotHLine_Callback(hObject, eventdata, handles)
%
% Description:
%    From a single point in the cone mosaic window, calculate the
%    horizontal cone responses per frame in a new window.
%
% Inputs:
%    hObject   - The handle to the time series
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
contextMenuPlot(hObject, []);
end

function menuPlotVLine_Callback(hObject, eventdata, handles)
% (Menu | Plot | V Line Response) Plot V Cone responses from a point
%
% Syntax:
%   menuPlotVLine_Callback(hObject, eventdata, handles)
%
% Description:
%    From a single point in the cone mosaic window, calculate the vertical
%    cone responses per frame in a new window.
%
% Inputs:
%    hObject   - The handle to the time series
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
contextMenuPlot(hObject, []);
end

function menuPlotHLineLMS_Callback(hObject, eventdata, handles)
% (Menu | Plot | H Line LMS Absorptions) Plot H LMS cone absorptions
%
% Syntax:
%   menuPlotHLineLMS_Callback(hObject, eventdata, handles)
%
% Description:
%    From a single point in the cone mosaic window, calculate the L, M, and
%    S absorptions per frame by cone horizontal positions in a new window
%    containing a vertical subplot for each cone type.
%
% Inputs:
%    hObject   - The handle to the time series
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
contextMenuPlot(hObject, []);
end

function menuPlotVLineLMS_Callback(hObject, eventdata, handles)
% (Menu | Plot | V Line LMS Absorptions) Plot V LMS cone absorptions
%
% Syntax:
%   menuPlotVLineLMS_Callback(hObject, eventdata, handles)
%
% Description:
%    From a single point in the cone mosaic window, calculate the L, M, and
%    S absorptions per frame by cone vertical positions in a new window
%    containing a vertical subplot for each cone type.
%
% Inputs:
%    hObject   - The handle to the time series
%    eventdata - reserved - to be defined in a future version of MATLAB
%    handles   - Structure with handles and user data
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
contextMenuPlot(hObject, []);
end

function menuPlotTimeSeries_Callback(hObject, ~, handles)
% (Menu | Plot | Time Series) A red circle deposited on cone mosaic window.
%
% Syntax:
%	menuPlotTimeSeries_Callback(hObject, ~, handles)
%
% Description:
%    Plot a single point in the window and create a graph of absorptions
%    per frame over time.
%
% Inputs:
%    hObject - The handle to the time series
%    handles - Structure with handles and user data
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% set(handles.btnPlayPause, 'Value', 0);  % Pause the movie
contextMenuPlot(hObject, []);
end
