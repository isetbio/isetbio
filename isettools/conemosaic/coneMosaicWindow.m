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

% Last Modified by GUIDE v2.5 28-Jun-2016 16:17:51

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

% --- Executes on button press in btnComputeImage.
function btnComputeImage_Callback(hObject, eventdata, handles)
% Button press computes the image from the optics data
[~, OI] = vcGetSelectedObject('OI');
if isempty(oiGet(OI, 'photons'))
    ieInWindowMessage('No optical image photon data.',handles);
    return
else
    ieInWindowMessage([],handles);
end
[val,ISA] = vcGetSelectedObject('ISA');

% The custom compute button is checked inside of sensorCompute.
ISA = sensorCompute(ISA,OI);

% For the moment, data and ISA are consistent.
ISA = sensorSet(ISA,'consistency',1);
vcReplaceObject(ISA,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuAnComputeFromOI_Callback(hObject, eventdata, handles)
btnComputeImage_Callback(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuAnComputeFromScene_Callback(hObject, eventdata, handles)

% Recompute the sensor data starting all the way back with the current
% scene and optical image.
[~, scene] = vcGetSelectedObject('scene');
if isempty(scene)
    warndlg('Creating default scene');
    scene = sceneCreate; val = 1;
    vcReplaceAndSelectObject(scene,val);
end

[val,oi] = vcGetSelectedObject('oi');
if isempty(oi),
    warndlg('Creating default OI');
    oi = oiCreate; val = 1;
end
oi = oiCompute(scene,oi);
vcReplaceAndSelectObject(oi,val);

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

% --------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)
end

function menuFileLoad_Callback(hObject, eventdata, handles, varargin)
end

function menuFileSave_Callback(hObject, eventdata, handles, varargin)
end

% --------------------------------------------------------------------
function menuFileLoadVoltsMat_Callback(hObject, eventdata, handles)
% Load Volts (MAT file)
% Load in new voltage data from a MAT file
%
fullName = vcSelectDataFile('stayput','r','mat');
if isempty(fullName), return; end

tmp = load(fullName,'volts');
isa = vcGetObject('ISA');
isa = sensorSet(isa,'volts',tmp.volts);

vcReplaceAndSelectObject(isa);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuFileSaveVoltsMat_Callback(hObject, eventdata, handles)
% Save Volts (MAT file)
isa = vcGetObject('ISA');
volts = sensorGet(isa,'volts');
if isempty(volts), errordlg('No voltage data.'); end

fullName = vcSelectDataFile('stayput','w','mat');
if ~isempty(fullName), save(fullName,'volts'); end

end

% --------------------------------------------------------------------
function menuSaveImage_Callback(hObject, eventdata, handles)
% Save Display (RGB Image)
%
% An option to save other data types will be needed some day, including DV
% in particular.

[val,isa] = vcGetSelectedObject('ISA');
gam = str2double(get(handles.editGam,'String'));
scaleMax = get(handles.btnDisplayScale,'Value');

sensorSaveImage(isa,[],'volts',gam,scaleMax);

end

% --------------------------------------------------------------------
function menuFileSpecsheet_Callback(hObject, eventdata, handles)
% File | Spec sheet
% Write out an Excel (or text?) spec sheet describing the sensor
%

disp('Spec sheet not yet implemented')

% The idea is to create a set of sensor spec values and then use
% xlswrite to printout an Excel spreadsheet summarizing them.

end

% --------------------------------------------------------------------
function menuFileClose_Callback(hObject, eventdata, handles)
sensorClose
end

% --------------------------------------------------------------------
function menuEdit_Callback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
function menuEditName_Callback(hObject, eventdata, handles)
[val,sensor] = vcGetSelectedObject('ISA');

newName = ieReadString('New sensor name','new-isa');
if isempty(newName), return;
else    sensor = sensorSet(sensor, 'name', newName);
end

vcReplaceObject(sensor,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuCopySensor_Callback(hObject, eventdata, handles)

sensor = vcGetObject('ISA');

newName = ieReadString('New ISA name','new-isa');
if isempty(newName),  return
else    sensor = sensorSet(sensor,'name',newName);
end

vcAddAndSelectObject('ISA',sensor);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuEditCreate_Callback(hObject, eventdata, handles)

createNewSensor(hObject, eventdata, handles);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuEditDelete_Callback(hObject, eventdata, handles)
sensorDelete(hObject,eventdata,handles);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuEditResWave_Callback(hObject, eventdata, handles)
% Edit | Resample Wavelength
isa = vcGetObject('isa');
isa = sensorResampleWave(isa);
vcReplaceObject(isa);
end

% --------------------------------------------------------------------
function sensorEditClearData_Callback(hObject, eventdata, handles)
% Edit | clear data
[val,ISA] = vcGetSelectedObject('ISA');
ISA = sensorClearData(ISA);
vcReplaceObject(ISA,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuEditClearMessage_Callback(hObject, eventdata, handles)
ieInWindowMessage('',ieSessionGet('sensorWindowHandles'),[]);
end

% --------------------------------------------------------------------
function menuEditZoom_Callback(hObject, eventdata, handles)
% Edit | Zoom
% Toggle the zoom state.  The zoom state is stored in the status of the
% checkbox of the zoom menu item.
if isequal(get(hObject,'checked'),'off'),
    zoom('on');
    set(hObject,'checked','on');
else % Must be on
    zoom('off');
    set(hObject,'checked','off');
end

end

% --------------------------------------------------------------------
function menuEditViewer_Callback(hObject, eventdata, handles)
% Edit | Viewer
sensor = vcGetObject('ISA');
img = sensorData2Image(sensor,'volts');
ieViewer(img);
end

% --- Executes on button press in btnAutoExp.
function btnAutoExp_Callback(hObject, eventdata, handles)
% Auto exposure button
[val,ISA] = vcGetSelectedObject('ISA');
ISA.AE = get(hObject,'Value');
vcReplaceObject(ISA,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuAnalyzePOVignetting_Callback(hObject, eventdata, handles)
% Analyze | Pixel Optics | Relative Illumination
sensor = vcGetObject('sensor');
sensorPlot(sensor,'etendue');
end
% --------------------------------------------------------------------
function menuSensorHumanCones_Callback(hObject, eventdata, handles)
% hObject    handle to menuSensorHumanCones (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function humanCones631_Callback(hObject, eventdata, handles)
% Sensor | Human cones | Cones-631
params.sz = [128,192];
params.rgbDensities = [0.0 .6 .3 .1]; % Empty, L,M,S
params.coneAperture = [2 2]*1e-6;     % In meters
sensor = sensorCreate('human', params);
vcAddAndSelectObject(sensor);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function humanConesKLMS1631_Callback(hObject, eventdata, handles)
% Sensor | Human cones | Cones-631
params.sz = [128,192];
params.rgbDensities = [0.1 .6 .3 .1]; % Empty, L,M,S
params.coneAperture = [2 2]*1e-6;     % In meters
sensor = sensorCreate('human', params);
vcAddAndSelectObject(sensor);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuSensorPixelVignetting_Callback(hObject, eventdata, handles)
% Sensor | Pixel OE Method
%
% Set check the Pixel OE computation.

ISA = vcGetObject('ISA');

% Set the vignetting
pvFlag = sensorGet(ISA,'pixelVignetting');
str = sprintf('Help: 0=skip,1=bare,2=centered,3=optimal');
pvFlag = ieReadNumber(str,pvFlag,'%.0f');
if isempty(pvFlag), return; end

ISA = sensorSet(ISA,'pixelvignetting',pvFlag);
ISA = sensorSet(ISA,'etendue',[]);

vcReplaceObject(ISA);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuSensorMicrolens_Callback(hObject, eventdata, handles)
% Sensor | Design Microlens
if ~exist('microLensWindow.m','file')
    errordlg('You do not have the optics toolkit.  Please contact ImagEval.');
else
    microLensWindow;
end

end

% --------------------------------------------------------------------
function menuSensorCDS_Callback(hObject, eventdata, handles)
%  Set check on CDS menu item.  We should probably display this in the
%  information box at the right, also.

ISA = vcGetObject('ISA');
state = get(hObject,'Check');

switch state
    case 'on'
        set(hObject,'Check','off');
        ISA = sensorSet(ISA,'cds',0);
    case 'off'
        set(hObject,'Check','on');
        ISA = sensorSet(ISA,'cds',1);
end
vcReplaceObject(ISA);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end


% --------------------------------------------------------------------
function menuSensorColFPN_Callback(hObject, eventdata, handles)

ISA = vcGetObject('ISA');
state = get(hObject,'Check');

switch state
    case 'on'
        set(hObject,'Check','off');
        ISA = sensorSet(ISA,'columnFPN',[]);
        ISA = sensorSet(ISA,'columnDSNU',[]);
        ISA = sensorSet(ISA,'columnPRNU',[]);
        % Should we clear the sensor data here?
    case 'off'
        % Store the column FPN value.
        set(hObject,'Check','on');
        
        % Should we pull this out as a separate routine?
        prompt={'Enter column DSNU (sd in millivolts)', ...
            'Enter column PRNU (sd. around unity gain)'};
        def={'1','0.01'};
        dlgTitle= sprintf('ISET read number');
        answer = inputdlg(prompt,dlgTitle,2,def);
        
        if   isempty(answer),  val = []; return;
        else
            colOffsetFPN = eval(answer{1})/1000;     % Read in mV, stored in volts
            colGainFPN   = eval(answer{2});          % Store as sd around unity gain
        end
        
        % Create and store the column noise parameters and an instance of
        % the noise itself
        nCol = sensorGet(ISA,'cols');
        colDSNU = randn(1,nCol)*colOffsetFPN;       % Offset noise stored in volts
        colPRNU = randn(1,nCol)*colGainFPN + 1;             % Column gain noise
        
        % Set the parameters.  Could combine the two reads into one.
        ISA = sensorSet(ISA,'columnFPN',[colOffsetFPN,colGainFPN]);
        ISA = sensorSet(ISA,'columnDSNU',colDSNU);
        ISA = sensorSet(ISA,'columnPRNU',colPRNU);
        
end
vcReplaceObject(ISA);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end


% --------------------------------------------------------------------
function menuSensorSetComputeGrid_Callback(hObject, eventdata, handles)
% Sensor | Set Compute Grid
% Sets pixel samples in signalCurrent

ISA = vcGetObject('ISA');
currentPixelSamples = sensorGet(ISA,'nSamplesPerPixel');
nPixelSamples = ieReadNumber('Enter odd integer specifying number of samples/pixel (default=1)',currentPixelSamples,'%.0f');

if isempty(nPixelSamples), return
elseif nPixelSamples < 1, nPixelSamples = 1;
elseif ~mod(nPixelSamples,2), nPixelSamples = nPixelSamples+1;
end

ISA = sensorSet(ISA,'nSamplesPerPixel',nPixelSamples);
vcReplaceObject(ISA);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --- Executes during object creation, after setting all properties.
function popupSelect_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
end

% --- Executes on selection change in popupSelect.
function popupSelect_Callback(hObject, eventdata, handles)
% Main popup coneMosaicWindow at top.  Create a New ISA or select a different one.
contents = get(hObject,'String');

switch  (contents{get(hObject,'Value')})
    
    case 'New'
        createNewSensor(hObject, eventdata, handles);
        
    otherwise,
        % The first two entries is always New.  The selections,
        % therefore, begin with the entry number in the list - 1.
        val = get(hObject,'Value') - 1;
        vcSetSelectedObject('ISA',val);
        coneMosaicGUIRefresh(hObject, eventdata, handles);
end

end


% --- Executes on selecting new from the popup.
function createNewSensor(hObject, eventdata, handles)
% New ISA
% Defaults to current values, except new color, as per sensorCreate.

[val, newISA] = vcGetSelectedObject('ISA');
newVal = vcNewObjectValue('ISA');

sensorArrayNames = get(handles.popISA,'String');
thisArrayName = sensorArrayNames{get(handles.popISA,'Value')};
switch lower(thisArrayName)
    case {'bayer-grbg'}
        newISA = sensorCreate('bayer-grbg');
    case {'bayer-rggb'}
        newISA = sensorCreate('bayer-rggb');
    case {'bayer-bggr'}
        newISA = sensorCreate('bayer-bggr');
    case 'bayer-ycmy'
        newISA = sensorCreate('bayer-ycmy');
    case {'Four Color'}
        newISA = sensorCreate('fourcolor');
    case 'monochrome'
        newISA = sensorCreate('monochrome');
    otherwise
        warning('Creating default sensor.')
        newISA = sensorCreate('default');
end

vcReplaceAndSelectObject(newISA,newVal);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --- Executes on button press in btnDelete.
function sensorDelete(hObject, eventdata, handles)

vcDeleteSelectedObject('ISA');
[val,isa] = vcGetSelectedObject('ISA');
if isempty(val)
    isa = sensorCreate('default');
    vcReplaceAndSelectObject(isa,1);
end

coneMosaicGUIRefresh(hObject, eventdata, handles);

end

function editDeleteSome_Callback(hObject, eventdata, handles)
% Edit delete some sensors
%
vcDeleteSomeObjects('sensor');
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
% Handle value is read during refresh.

% Don't change the red consistency button
sensor = vcGetObject('sensor');
sensor = sensorSet(sensor,'consistency',-1);
vcReplaceObject(sensor);
coneMosaicGUIRefresh(hObject,eventdata,handles);
end

% --- Executes on button press in btnDisplayScale.
function btnDisplayScale_Callback(hObject, eventdata, handles)
% Handle value is read during refresh.

% Don't change the red consistency button
sensor = vcGetObject('sensor');
sensor = sensorSet(sensor,'consistency',-1);
vcReplaceObject(sensor);

coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function popupColor_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
end

% --- Executes during object creation, after setting all properties.
function popFormat_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
end

% --------------------------------------------------------------------
function menuFileRefresh_Callback(hObject, eventdata, handles)
coneMosaicGUIRefresh(hObject, eventdata, handles);
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

end

% --- Executes during object creation, after setting all properties.
function popQuantization_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
end

% --- Executes on selection change in popQuantization.
function popQuantization_Callback(hObject, eventdata, handles)

contents = get(hObject,'String');
qMethod = contents{get(hObject,'Value')};
[val,isa] = vcGetSelectedObject('ISA');

isa = sensorSet(isa,'quantization',qMethod);

isa = sensorClearData(isa);
vcReplaceObject(isa,val);

coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuSensor_Callback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
function menuSensorCIF_Callback(hObject, eventdata, handles)

end
% --------------------------------------------------------------------
function menuSensorQQCIFSixteenthInch_Callback(hObject, eventdata, handles)
%

[val,isa] = vcGetSelectedObject('ISA');
isa = sensorRescale(isa,sensorFormats('qqcif'),sensorFormats('sixteenthinch'));
vcReplaceObject(isa,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuSensorQQCIFQuartInch_Callback(hObject, eventdata, handles)
%

[val,isa] = vcGetSelectedObject('ISA');
isa = sensorRescale(isa,sensorFormats('qqcif'),sensorFormats('quarterinch'));
vcReplaceObject(isa,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuSensorQCIFQuartInch_Callback(hObject, eventdata, handles)
%
% QCIF format, half-inch sensor

[val,isa] = vcGetSelectedObject('ISA');
isa = sensorRescale(isa,sensorFormats('qcif'),sensorFormats('quarterinch'));
vcReplaceObject(isa,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuSensorCIFHalfInch_Callback(hObject, eventdata, handles)
%
[val,isa] = vcGetSelectedObject('ISA');
isa = sensorRescale(isa,sensorFormats('cif'),sensorFormats('halfinch'));
vcReplaceObject(isa,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuSensorVGA_Callback(hObject, eventdata, handles)
% Sensor->VGA
end

% --------------------------------------------------------------------
function menuSensorQQVGAQuartInch_Callback(hObject, eventdata, handles)
%
[val,isa] = vcGetSelectedObject('ISA');
isa = sensorRescale(isa,sensorFormats('qqvga'),sensorFormats('quarterinch'));
vcReplaceObject(isa,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuSensorQVGAQuartInch_Callback(hObject, eventdata, handles)
%
[val,isa] = vcGetSelectedObject('ISA');
isa = sensorRescale(isa,sensorFormats('qvga'),sensorFormats('quarterinch'));
vcReplaceObject(isa,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuSensorQVGAHalfInch_Callback(hObject, eventdata, handles)
[val,isa] = vcGetSelectedObject('ISA');
isa = sensorRescale(isa,sensorFormats('qvga'),sensorFormats('halfinch'));
vcReplaceObject(isa,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuSensorVGAHalfInch_Callback(hObject, eventdata, handles)
%
[val,isa] = vcGetSelectedObject('ISA');
isa = sensorRescale(isa,sensorFormats('vga'),sensorFormats('halfinch'));
vcReplaceObject(isa,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuSensorQQVGASixteenthInch_Callback(hObject, eventdata, handles)
%
[val,isa] = vcGetSelectedObject('ISA');
isa = sensorRescale(isa,sensorFormats('qqvga'),sensorFormats('sixteenthinch'));
vcReplaceObject(isa,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuDesignCFA_Callback(hObject, eventdata, handles)
%      I first use the 'Design CFA' option in the Sensor Menu to design the CFA I desire. I make
%      sure I position the individual color filters appropriately in the 2x2 square as 'RGGB'.
%
%      Then after I load all the respective transmissivity values for each of the filters, I save
%      the CFA to a file (say we call this file micronRGGB.mat).
%
%      Now during a simulation run, when I load this CFA using the "Load CFA" option in the Sensor
%      menu, the "Standard CFA" cell in the upper right hand side of the ISET-Sensor conemosaicwindow
%      gets automatically updated as "Other" instead of 'bayer-rggb' as I would expect.
%
sensorDesignCFA;
end


% --------------------------------------------------------------------
function menuSensorExport_Callback(hObject, eventdata, handles)
% File | Save | Save sensor (.mat)
[val,isa] = vcGetSelectedObject('ISA');
fullName = vcExportObject(isa);
end

% --------------------------------------------------------------------
function menuSensorImport_Callback(hObject, eventdata, handles)
% File | Load  | Sensor (.mat)
%
newVal = vcImportObject('ISA');
if isempty(newVal), return; end
vcSetSelectedObject('ISA',newVal);
coneMosaicGUIRefresh(hObject, eventdata, handles);
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

% ---------------Plot Menu------------------------
function menuPlot_Callback(hObject, eventdata, handles)
% Menu Plot
end
% --------------------------------------------------------------------
function menuPlotSpectra_Callback(hObject, eventdata, handles)
% Plot-> SpectralInformation
end

% --------------------------------------------------------------------
function menuPlotPixelSR_Callback(hObject, eventdata, handles)
sensorPlot([],'pixel spectral sr');
end

% --------------------------------------------------------------------
function menuPlotColorFilters_Callback(hObject, eventdata, handles)
sensorPlot([],'color filters');
end

% --------------------------------------------------------------------
function plotSpecCFApattern_Callback(hObject, eventdata, handles)
% Plot | Spectral Information | CFA Pattern
sensorPlot([],'CFA');
end
% --------------------------------------------------------------------
function menuPlotIR_Callback(hObject, eventdata, handles)
sensorPlot([],'irfilter');
end

% --------------------------------------------------------------------
function menuPlotPDSpectralQE_Callback(hObject, eventdata, handles)
sensorPlot([],'pixel spectral QE');
end

% --------------------------------------------------------------------
function menuPlotSpecResp_Callback(hObject, eventdata, handles)
sensorPlot([],'sensor spectral qe');
end

% --------------------------------------------------------------------
function menusensorPlotImageTSize_Callback(hObject, eventdata, handles)
% Plot | SensorImage (True Size)
gam      = str2double(get(handles.editGam,'String'));
scaleMax = get(handles.btnDisplayScale,'Value');
sensor   = vcGetObject('sensor');

% Get voltages or digital values
bits     = sensorGet(sensor,'bits');
if isempty(bits)
    img      = sensorData2Image(sensor,'voltage',gam,scaleMax);
else
    img      = sensorData2Image(sensor,'dv',gam,scaleMax);
end

if ndims(img) == 2;
    % imtool needs monochrome images scaled between 0 and 1
    w = vcNewGraphWin; img = img/max(img(:));
    imshow(img); truesize(w);
    set(w,'Name',sensorGet(sensor,'name'));
else
    ieViewer(img);
end


end

% --------------------------------------------------------------------
function plotMccOverlay_Callback(hObject, eventdata, handles)
% Plot | MCC overlay off
% Delete the MCC boxes showing the selection

sensor = vcGetObject('sensor');
macbethDrawRects(sensor,'off');
vcReplaceObject(sensor);
coneMosaicGUIRefresh(hObject, eventdata, handles);

end

% --------------------------------------------------------------------
function menuPlotHumanCone_Callback(hObject, eventdata, handles)
% Plot | Human Cone

sensor = vcGetObject('sensor');
if sensorCheckHuman(sensor), sensorConePlot(sensor)
else ieInWindowMessage('Not a human cone sensor',handles,3);
end

end

% --------------------------------------------------------------------
function menuPlotNewGraphWindow_Callback(hObject, eventdata, handles)
vcNewGraphWin;
end

% --------------------------------------------------------------------
function menuAnalyze_Callback(hObject, eventdata, handles)
% Analyze menu
end

% --------------------------------------------------------------------
function menuAnLine_Callback(hObject, eventdata, handles)
% Analyze->Line menu
end

% --------------------------------------------------------------------
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
function menuAnPixelSNR_Callback(hObject, eventdata, handles)
% Graph the pixel SNR as a function of voltage swing
sensorPlot(vcGetObject('sensor'),'pixel snr');
% plotPixelSNR;
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

% ----------------------Pixel Optics-------------------
function menuAnPO_Callback(hObject, eventdata, handles)
% Analysis->Pixel Optics
end


% --------------------------------------------------------------------
function menuAnPixOptLoadUL_Callback(hObject, eventdata, handles)
% Analyze | Pixel Optics | Load uL

fullName = vcSelectDataFile('stayPut','r');
if isempty(fullName), disp('User canceled'); return; end
tmp = load(fullName);
if isfield(tmp,'ml')
    ISA = vcGetObject('ISA');
    ISA = sensorSet(ISA,'microLens',tmp.ml);
    vcReplaceObject(ISA);
else
    error('No microlens structure (named ml) in the file.');
end

end

% --------------------------------------------------------------------
function menuAnPixOptSaveUl_Callback(hObject, eventdata, handles)
% Analyze | Pixel Optics | Save uL

ISA = vcGetObject('ISA');
ml = sensorGet(ISA,'microLens');
fullName = vcSelectDataFile('stayPut','w');
if isempty(fullName), disp('User canceled'); return; end
save(fullName,'ml');

end

% --------------------------------------------------------------------
function menuAnColor_Callback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
function menuAnColorRB_Callback(hObject, eventdata, handles)
% Analyze | Color | RB Analysis
sensorPlotColor(vcGetObject('ISA'),'rb');
end

% --------------------------------------------------------------------
function menuAnColorRG_Callback(hObject, eventdata, handles)
% Analyze | Color | RG Analysis
sensorPlotColor(vcGetObject('ISA'),'rg');
end

% --------------------------------------------------------------------
function menuAnColCCM_Callback(hObject, eventdata, handles)
% Analyze | Color | Color Conversion Matrix

sensor = vcGetObject('sensor');
[L,corners] = sensorCCM(sensor); %#ok<ASGLU>

fprintf('    ==  MCC to XYZ_D65 matrix  ==\n');
disp(L)

% Store the selection of the corners
sensor = sensorSet(sensor,'mcc corner points',corners);
vcReplaceObject(sensor);

end

% --------------------------------------------------------------------
function menuEdgeOp_Callback(hObject, eventdata, handles)
% Start a new process to initiate the edge operator coneMosaicWindow.
% This coneMosaicWindow works with monochrome sensor images and investigate how
% various operators perform with different types of sensors.
edgeOperatorWindow;
end

% ------------------PlotImage------------------------------
function menuIm_Callback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
function menuImDSNU_Callback(hObject, eventdata, handles)
sensorPlot([],'dsnu');
end

% --------------------------------------------------------------------
function menuImPRNU_Callback(hObject, eventdata, handles)
sensorPlot([],'prnu');
end

% --------------------------------------------------------------------
function menuImShotNoise_Callback(hObject, eventdata, handles)
sensorPlot([],'shotnoise');
end

% --------------------------------------------------------------------
function menuSensorLuxSec_Callback(hObject, eventdata, handles)
sensorSNRluxsec;
end

% --------------------------------------------------------------------
function menuPixelLuxSec_Callback(hObject, eventdata, handles)
pixelSNRluxsec;
end

% --------------------------------------------------------------------
function menuPDSize_Callback(hObject, eventdata, handles)
pixelGeometryWindow;
end

% --------------------------------------------------------------------
function menuPixelLayers_Callback(hObject, eventdata, handles)
pixelOEWindow;
end

% --------------------------------------------------------------------
function menuLoadSF_Callback(hObject, eventdata, handles)
% Load spectral functions heading.
end

% --------------------------------------------------------------------
function menuLoadCFA_Callback(hObject, eventdata, handles)
%
[val,ISA] = vcGetSelectedObject('ISA');
ISA = sensorReadFilter('cfa',ISA);
ISA = sensorClearData(ISA);
vcReplaceObject(ISA,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuColorFilters_Callback(hObject, eventdata, handles)
% Sensor | Load Color Filters
[val,ISA] = vcGetSelectedObject('ISA');
ISA = sensorReadFilter('colorfilters',ISA);
ISA = sensorClearData(ISA);
vcReplaceObject(ISA,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuInfraRed_Callback(hObject, eventdata, handles)
[val,ISA] = vcGetSelectedObject('ISA');
ISA = sensorReadFilter('infrared',ISA);
ISA = sensorClearData(ISA);
vcReplaceObject(ISA,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuPDQE_Callback(hObject, eventdata, handles)
[val,ISA] = vcGetSelectedObject('ISA');
ISA = sensorReadFilter('pdspectralqe',ISA);
ISA = sensorClearData(ISA);
vcReplaceObject(ISA,val);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuEditFontSize_Callback(hObject, eventdata, handles)
ieFontChangeSize(handles.coneMosaicWindow);
end

% --------------------------------------------------------------------
function menuAnalyzeMicroLens_Callback(hObject, eventdata, handles)
% Analyze | Pixel Optics | Microlens coneMosaicWindow
if isempty(which('microLensWindow'))
    warndlg('No micro lens analysis software on your path.  Contact ImagEval for a license.');
    return;
else
    % Should test for license here
end

microLensWindow;
end

% --------------------------------------------------------------------
function menuAnPOShowUL_Callback(hObject, eventdata, handles)
ISA = vcGetObject('ISA');
mlPrint(sensorGet(ISA,'microLens'));
end

% --------------------------------------------------------------------
function menuHelp_Callback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
function menuAppNotes_Callback(hObject, eventdata, handles)
% Help | Documentation (web)
web('http://imageval.com/documentation/','-browser');
end

% --------------------------------------------------------------------
function menuHelpSensorOnline_Callback(hObject, eventdata, handles)
% Help | Sensor (online)
web('http://www.imageval.com/public/ISET-Functions/ISET/sensor/index.html','-browser');
end

% --------------------------------------------------------------------
function menuHelpPixelOnline_Callback(hObject, eventdata, handles)
% Help | Pixel (online)
web('http://www.imageval.com/public/ISET-Functions/ISET/sensor/pixel/index.html','-browser');
end

% --------------------------------------------------------------------
function menuHelpISETOnlineManual_Callback(hObject, eventdata, handles)
% Help | ISET (online)
web('http://www.imageval.com/public/ISET-Functions/','-browser');
end

% --- Executes during object creation, after setting all properties.
function popupExpMode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end
% --- Executes on selection change in popupExpMode.
function popupExpMode_Callback(hObject, eventdata, handles)
% Exposure popup control
%
% Hints: contents = get(hObject,'String') returns popupExpMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupExpMode

sensor = vcGetObject('sensor');

% Determine which case we popped into
contents = get(hObject,'String');
switch contents{get(hObject,'Value')}
    case 'Single'
        set(handles.btnAutoExp,'visible','on');
        set(handles.editExpTime,'visible','on');
        set(handles.btnShowCFAExpDurations,'visible','off');
        set(handles.editNExposures,'visible','off');
        set(handles.editExpFactor,'visible','off');
        set(handles.sliderSelectBracketedExposure,'visible','off');
        set(handles.txtBracketExposure,'visible','off');
        
        eTime  = sensorGet(sensor,'geometricMeanExposureTime');
        sensor = sensorSet(sensor,'expTime',eTime);
        
    case 'Bracketing'
        set(handles.editExpTime,'visible','on');
        set(handles.sliderSelectBracketedExposure,'visible','on');
        set(handles.editNExposures,'visible','on');
        set(handles.editExpFactor,'visible','on');
        set(handles.btnShowCFAExpDurations,'visible','off');
        
        set(handles.btnAutoExp,'visible','off');
        set(handles.txtBracketExposure,'visible','on');
        
        % Manage the bracket times
        sensor = sensorAdjustBracketTimes(handles,sensor);
        
    case 'CFA Exposure'
        set(handles.btnAutoExp,'visible','on');
        set(handles.btnShowCFAExpDurations,'visible','on');
        set(handles.editExpTime,'visible','off');
        set(handles.editNExposures,'visible','off');
        set(handles.editExpFactor,'visible','off');
        set(handles.sliderSelectBracketedExposure,'visible','off');
        set(handles.txtBracketExposure,'visible','off');
        
        sensor = sensorAdjustCFATimes(handles,sensor);
        
    otherwise
        error('Unknown exposure condition %s\n',contents{get(hObject,'Value')});
end

sensor = sensorClearData(sensor);
vcReplaceObject(sensor);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function editNExposures_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function editNExposures_Callback(hObject, eventdata, handles)
% Set the number of exposures in the bracketing mode
sensor = vcGetObject('sensor');
sensor = sensorAdjustBracketTimes(handles,sensor);
sensor = sensorClearData(sensor);
vcReplaceObject(sensor);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function editExpFactor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function editExpFactor_Callback(hObject, eventdata, handles)
% Set the scale factor between exposures in bracketing mode
%
% Hints: get(hObject,'String') returns contents of editExpFactor as text
%        str2double(get(hObject,'String')) returns contents of editExpFactor as a double
sensor = sensorAdjustBracketTimes(handles);
sensor = sensorClearData(sensor);
vcReplaceObject(sensor);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --- Executes on button press in btnShowCFAExpDurations.
function btnShowCFAExpDurations_Callback(hObject, eventdata, handles)
% Bring up a gui that shows the CFA exposures

sensor = vcGetObject('sensor');

% Make sure we are in the CFA format
sensor = sensorAdjustCFATimes(handles,sensor);

% Put the exposureData matrix into the base workspace in ms
fmt    = '%.1f';
prompt = 'Time (ms)';
defMatrix = sensorGet(sensor,'expTime')*1e3;
saturation = 0.3;
filterRGB = sensorFilterRGB(sensor,saturation);
ieReadSmallMatrix(size(defMatrix),defMatrix,fmt,prompt,[],'msExposureData',filterRGB);

%Read base space data in seconds
secExposureData = evalin('base','msExposureData')*1e-3;

% Put the sensor back
sensor = sensorSet(sensor,'expTime',secExposureData);
vcReplaceObject(sensor);

end

% --- Executes during object creation, after setting all properties.
function sliderSelectBracketedExposure_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
    
end
end

% --- Executes on slider movement.
function sliderSelectBracketedExposure_Callback(hObject, eventdata, handles)
% Chooses which of the bracketed exposures should be display
sensor = vcGetObject('sensor');

exposurePlane = get(handles.sliderSelectBracketedExposure,'value');
sensor = sensorSet(sensor,'exposurePlane',exposurePlane);
sensor = sensorSet(sensor,'consistency',-1);  % Don't change the red consistency button
vcReplaceObject(sensor);
coneMosaicGUIRefresh(hObject, eventdata, handles);
end

% --- Executes any time we need to update the bracketed exposures.
function sensor = sensorAdjustBracketTimes(handles,sensor)
% Examine the bracket settings edit boxes and adjust the exposure times
%

% This can be called by editing the exp time box or the nExposures or the
% scale factor on exposures.
%
% If it is called from the exp time box, then we get the sensor passed in.
% Otherwise we use the default sensor.
if notDefined('sensor'), sensor = vcGetObject('sensor'); end

% Managing the GUI should probably be done inside of the refresh command,
% not here.

% Require odd exposure number.  Update the edit box with the right number
nExposures = str2double(get(handles.editNExposures,'String'));
if ~isodd(nExposures), nExposures = nExposures + 1; end
set(handles.editNExposures,'String',num2str(nExposures));

% Create the exposure list from GUI data
sFactor         = str2double(get(handles.editExpFactor,'String'));
centralExposure = sensorGet(sensor,'Geometric Mean Exposure Time');
nBelow          = floor(nExposures/2);
shortExposure   = centralExposure/(sFactor^nBelow);
expTimes = zeros(1,nExposures);
for ii=1:nExposures
    expTimes(ii) = shortExposure*sFactor^(ii-1);
end

% Update the sensor
sensor        = sensorSet(sensor,'Exp Time',expTimes);
exposurePlane = floor(nExposures/2) + 1;
sensor = sensorSet(sensor,'Exposure Plane',exposurePlane);

% Set the slider for bracketed exposures now done in sensorEditsAndButtons
% set(handles.sliderSelectBracketedExposure,'max',nExposures);
% set(handles.sliderSelectBracketedExposure,'value',exposurePlane);
% if nExposures > 1, ss = 1/(nExposures-1);
% else               ss = 0; end
% set(handles.sliderSelectBracketedExposure,'sliderStep',[ss ss]);

end

% --- Executes any time we need to update the CFA exposures
function sensor = sensorAdjustCFATimes(handles,sensor)
% Adjust the sensor exposure time slot to match the CFA size, putting the
% sensor into the CFA exposure mode.
%

if notDefined('sensor'), sensor = vcGetObject('sensor'); end

mSize = size(sensorGet(sensor,'pattern'));
eTimes = sensorGet(sensor,'expTimes');

% If the eTimes has the wrong size, use the geometric mean of the current
% eTime values and assign it to a matrix of the right size.
if ~isequal(size(eTimes),mSize)
    eTimes = ones(mSize)*sensorGet(sensor,'geometricMeanExposuretime');
    sensor = sensorSet(sensor,'expTime',eTimes);
end

%#ok<*DEFNU>
%#ok<*INUSD>
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
