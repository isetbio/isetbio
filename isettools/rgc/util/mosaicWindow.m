function varargout = mosaicWindow(varargin)
% MOSAICWINDOW MATLAB code for mosaicWindow.fig
%      MOSAICWINDOW, by itself, creates a new MOSAICWINDOW or raises the existing
%      singleton*.
%
%      H = MOSAICWINDOW returns the handle to a new MOSAICWINDOW or the handle to
%      the existing singleton*.
%
%      MOSAICWINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOSAICWINDOW.M with the given input arguments.
%
%      MOSAICWINDOW('Property','Value',...) creates a new MOSAICWINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mosaicWindow_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mosaicWindow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mosaicWindow

% Last Modified by GUIDE v2.5 28-Oct-2016 13:46:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mosaicWindow_OpeningFcn, ...
                   'gui_OutputFcn',  @mosaicWindow_OutputFcn, ...
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

% --- Executes just before mosaicWindow is made visible.
function mosaicWindow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mosaicWindow (see VARARGIN)

% check inputs and get the rgcMosaic object
if isempty(varargin) || ~isa(varargin{1}, 'rgcMosaic')
    error('rgc mosaic object required');
end
rgcM = varargin{1};
rgcM.figureHandle = hObject;   % Store this figure handle

% Choose default command line output for mosaicWindow
handles.output = hObject;
handles.rgcMosaic = varargin{1};
handles.spikes = [];  % spike movie

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mosaicWindow wait for user response (see UIRESUME)
% uiwait(handles.mosaicWindow);

% Refresh/Initialize window information
mosaicWindowRefresh(handles);

% Very important for good rendering speed
set(hObject, 'Renderer', 'OpenGL')

handles.linearMov = [];
handles.psthMov = [];
end

% --- Outputs from this function are returned to the command line.
function varargout = mosaicWindow_OutputFcn(hObject, eventdata, handles) 
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

handles.linearMov = [];
handles.psthMov = [];
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

% These are all the strings in the popup
contents = cellstr(get(hObject,'String'));

% This is the selected string
str = contents{get(hObject,'Value')};

% Perform the action given the selection
switch str
    case 'Receptive field mosaic'
        % Should be ellipsoids
        handles.rgcMosaic.plot('mosaic')
        axis image;
        
    case 'Spike movie'
        % Spike Movie pulldown
        psthTest = handles.rgcMosaic.get('spikes');
        clear vParams; vParams = [];
        vParams.FrameRate = 30; vParams.show = true; %vParams.step = 2; 
        frameSkip = round(1./handles.rgcMosaic.get('dt'));
        ieMovie(psthTest(:,:,1:frameSkip:end),vParams);
        axis image
        
    case 'Linear movie'
        % Linear movie, needs units if possible, colorbar if possible
        responseLinear = handles.rgcMosaic.get('responseLinear');
        clear vParams; vParams = [];
        vParams.FrameRate = 30; vParams.show = true; %vParams.step = 2; 
        
        % Should we add controls like the cone mosaic window, or just play
        % it once?
        ieMovie(responseLinear,vParams);

        %         set(handles.btnPlayPause, 'Visible', 'on');
        %         set(handles.btnPlayPause, 'Value', 1);  % Auto start the movie
        %         set(handles.sliderMovieProgress, 'Visible', 'on');
        %         if ~exist('handles.linearMov','var')
        % %         if isempty(handles.linearMov)
        % %             ieInWindowMessage('Building movie',handles,2);
        %             % generate movie
        % %             handles.mov = cm.plot('movie absorptions', 'hf','none',...
        % %                 'show',true, ...
        % %                 'gamma', str2double(get(handles.editGam, 'String')));
        %             handles.linearMov = ieMovie(responseLinear,vParams);
        %             guidata(hObject, handles);
        %         end
        %
        %         % play movie if more than one frame
        %         btnPlayPause_Callback(hObject, eventdata, handles);
        
    case 'Spike mean (image)'
        disp(str)
        spikes = handles.rgcMosaic.get('spikes');
        imagesc(mean(spikes,3)); 
        colorbar; axis image; set(gca,'xticklabels','','yticklabels','');
        drawnow;
        
    case 'Linear plot'
        disp(str)        
        responseLinear = handles.rgcMosaic.get('responseLinear');
        plot(RGB2XWFormat(responseLinear)');
        
    case 'PSTH plot'
        responsePsth = handles.rgcMosaic.get('psth');
        plot(RGB2XWFormat(responsePsth)');
        
    case 'PSTH movie'
        % PSTH movie shows all the cells as a PSTH
        responsePsth = handles.rgcMosaic.get('psth');
        clear vParams; vParams = [];
        vParams.FrameRate = 30; vParams.show = true; %vParams.step = 2;
        frameSkip = round(1./handles.rgcMosaic.get('dt'));
        
        % Same as above ... loop?  Control pause?
        ieMovie(responsePsth(:,:,1:frameSkip:end),vParams);
        
        %         set(handles.btnPlayPause, 'Visible', 'on');
        %         set(handles.btnPlayPause, 'Value', 1);  % Auto start the movie
        %         set(handles.sliderMovieProgress, 'Visible', 'on');
        %         if ~exist('handles.psthMov','var')
        % %         if isempty(handles.psthMov)
        % %             ieInWindowMessage('Building movie',handles,2);
        %             % generate movie
        % %             handles.mov = cm.plot('movie absorptions', 'hf','none',...
        % %                 'show',true, ...
        % %                 'gamma', str2double(get(handles.editGam, 'String')));
        %             handles.psthMov = ieMovie(responsePsth,vParams);
        %             guidata(hObject, handles);
        %         end
        %
        %         % play movie if more than one frame
        %         btnPlayPause_Callback(hObject, eventdata, handles);
        
    case 'PSTH mean (image)'
        % This is odd.  Maybe it should exist.
        psthTest = handles.rgcMosaic.get('psth');
        imagesc(mean(psthTest,3)); drawnow
        
    otherwise
        error('Unknown string %s\n',str);
end

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
function menuPlotPSTH_Callback(hObject, eventdata, handles)
% Plot | PSTH
%
% hObject    handle to menuPlotPSTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Plot | PSTH')
end

% --------------------------------------------------------------------
function menFileSave_Callback(hObject, eventdata, handles)
% File | Save
%
% hObject    handle to menFileSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Save')
end

% --------------------------------------------------------------------
function menuFileRefresh_Callback(hObject, eventdata, handles)
% hObject    handle to menuFileRefresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mosaicWindowRefresh(handles)
end

% --------------------------------------------------------------------
function menuFileClose_Callback(hObject, eventdata, handles)
% File | Close
% hObject    handle to menuFileClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.rgcM.figureHandle = [];
delete(handles.mosaicWindow);
end

%% Internal functions

function mosaicWindowRefresh(handles)
% Update all the text fields and such with the data in the mosaic
disp('Refresh')

rgcM  = handles.rgcMosaic;
fig   = figure(rgcM.figureHandle);
gdata = guidata(fig);
% Show the appropriate response axis plot
axes(gdata.axisResponse);
cla(gdata.axisResponse,'reset');

% Switch depending on the state of the pull down
spikes = rgcM.get('response spikes');
img = mean(spikes(:,:,:,1),3);
imagesc(img); colormap(gray); colorbar;
xlabel('Distance (um)');

% RF shape overlay on the response window
% % Update the geometry axis plot
% axes(gdata.axisGeometry);
% cla(gdata.axisGeometry,'reset');

% Make a button for rfOverlay
rfOverlay = false;
if rfOverlay
    rgcM.plot('mosaic');
end


% Text description
str = rgcM.describe;
set(gdata.rgcProperties,'string',str);

end


% --- Executes on button press in btnPlayPause.
function btnPlayPause_Callback(hObject, eventdata, handles)
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