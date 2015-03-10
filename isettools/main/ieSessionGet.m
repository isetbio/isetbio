function val = ieSessionGet(param)
% Get fields in vcSESSION, including figure handles and custom routines
%
%  val = ieSessionGet(param);
%
%  The vcSESSION parameter is a global variable at present.  It contains
%  information about the windows, custom processing routines, and related
%  ISET processing information.
%
%  This get routine retrieves that information.  The information is stored
%  in a global variable for now.  In the future, this information will be
%  obtained using findobj().
%
%  The tag 'handle' refers to the guihandles.  The tag 'figure' refers to
%  the figure number.  The guidhandles can be retrieved by using 
%  h = guihanles(f);
%
%  A list of the parameters is:
%
%      {'version'}
%      {'name','session name'}
%      {'dir','session dir'}
%      {'help','init help'}
%
%  Matlab pref variables
%      {'delta font size'} -- This value determines whether we change the
%        font size in every window by this increment, calling
%        ieFontChangeSize when the window is opened. 
%      {'waitbar'}  - Show compute waitbars or not.  
%
%  Figure handles
%      {'graphwin structure'}
%      {'graphwin figure'}
%      {'main figure'}
%      {'scene figure'}
%      {'oi figure'}
%      {'sensor figure'}
%      {'vcimage figure'}
%      {'metrics figure'}
%
% Guidata
%      {'main  guidata'}
%      {'scene  guidata'}
%      {'oi  guidata'}
%      {'sensor  guidata'}
%      {'ip guidata'}
%      {'metrics  guidata'}
%      {'graphwin guidata'}
%
%  ** Custom algorithms list is deprecated
%      {'custom','customall','customstructure'}
%         val = vcSESSION.CUSTOM;
%         % These are cell arrays of user-defined routines that implement
%         % these various types of operations.
%      {'customdesmoaiclist'}
%      {'customcolorbalancelist'}
%      {'customcolorconversionlist'}
%      {'processingmethods'}
%         % These routines are a complete processing chain that replace the
%         % entire call to vcimageCompute
%      {'oicomputelist'}
%         % These routines replace the standard oiCompute call.  They
%         % customize the signal processing chain from the optical image to
%         % the OI data.
%      {'sensorcomputelist'}
%         % These routines replace the standard sensorCompute call.  They
%         % customize the signal processing chain from the optical image to
%         % the ISA data.
%      {'edgealgorithmlist'}
%         % These routines replace the standard sensorCompute call.  They
%         % customize the signal processing chain from the optical image to
%         % the ISA data.
%
% Example:
%   h = ieSessionGet('scene window handle')
%   f = ieSessionGet('scene figure')
%   guihandles(f)
%
%   ieSessionGet('version')
%   d = ieSessionGet('fontsize'); ieFontChangeSize(sceneWindow,d);
%
%   hobj = ieSessionGet('opticalimagefigure');
%
% Copyright ImagEval Consultants, LLC, 2005.

global vcSESSION

if notDefined('param'), error('You must specify a parameter.'); end
val = [];

% Eliminate spaces and make lower case
param = ieParamFormat(param);

switch param
    case {'version'}
        val = vcSESSION.VERSION;
    case {'name','sessionname'}
        val = vcSESSION.NAME;
    case {'dir','sessiondir'}
        val = vcSESSION.DIR;
    case {'help','inithelp'}
        % Default for help is true, if the initHelp has not been set.
        if checkfields(vcSESSION,'initHelp'), val = vcSESSION.initHelp; 
        else vcSESSION.initHelp = 1; val = 1; 
        end
        
    % Matlab setpref/getpref 
    case {'detlafontsize', 'fontsize', 'fontincrement', ...
            'increasefontsize', 'fontdelta', 'deltafont'}
        % This value determines whether we change the font size in every
        % window by this increment, calling ieFontChangeSize when the
        % window is opened. if checkfields(vcSESSION,'FONTSIZE'), val =
        % vcSESSION.FONTSIZE;  end
        isetPref = getpref('ISET');
        if ~isempty(isetPref)
            if checkfields(isetPref,'fontDelta'), val = isetPref.fontDelta; 
            end
        else 
            val = 0; 
        end
        if isempty(val), val = 0; end
        
    case {'waitbar'}
        % Used to decide whether we show the waitbars.
        if checkfields(vcSESSION,'GUI','waitbar')
            val = vcSESSION.GUI.waitbar;
        else
            % The getpref is slow.  So, we attach it to the session 
            % at start up.  Otherwise, loops that test for it take too
            % long.
            iePref = getpref('ISET');
            if ~checkfields(iePref,'waitbar')
                setpref('ISET','waitbar',0);
                val = 0;
            else val = iePref.waitbar;
            end
            vcSESSION.GUI.waitbar = val;
        end
        
    case {'gpu', 'gpucompute', 'gpucomputing'}
        % Whether or not to use gpu compute.  Always false now, but in the
        % future we may do more with this.
        val = false;
    case {'imagesizethreshold'}
        % Used by the display code.  This sets a value for when we loop
        % over wavelength instead of doing a large matrix multiplication.
        % HJ - more comments later.
        if isfield(vcSESSION, 'imagesizethreshold')
            val = vcSESSION.imagesizethreshold;
        else
            val = 1e6;
        end
    case {'graphwinstructure'}
        val = vcSESSION.GRAPHWIN;
    case {'graphwinfigure'}
        if checkfields(vcSESSION,'GRAPHWIN','hObject') 
            val = vcSESSION.GRAPHWIN.hObject; 
        end  
    case {'graphguidata','graphwinhandles','graphwinhandle'}
        if checkfields(vcSESSION,'GRAPHWIN','handle') 
            val = vcSESSION.GRAPHWIN.handle; 
        end  
        
        % Handles to the various windows
    case {'mainguidata','mainwindowhandle','mainhandle','mainhandles'}
        v = ieSessionGet('mainfigure');
        if ~isempty(v), val = guihandles(v); end
    case {'sceneguidata','scenewindowhandle','sceneimagehandle','sceneimagehandles','scenewindowhandles'}
        v = ieSessionGet('sceneimagefigure');
        if ~isempty(v), val = guihandles(v); end
    case {'oiguidata','oiwindowhandle','oihandle','opticalimagehandle','oihandles','opticalimagehandles','oiwindowhandles'}
        v = ieSessionGet('opticalimagefigure');
        if ~isempty(v), val = guihandles(v); end
    case {'sensorguidata','sensorwindowhandle','sensorimagehandle','sensorhandle','isahandle','sensorhandles','isahandles','sensorwindowhandles'}
        v = ieSessionGet('sensorfigure');
        if ~isempty(v), val = guihandles(v); end
    case {'ipguidata','vciguidata','vciwindowhandle','vcimagehandle','vcimagehandles','processorwindowhandles','processorhandles','processorhandle','processorimagehandle'}
        v = ieSessionGet('vcimagefigure');
        if ~isempty(v), val = guihandles(v); end
    case {'metricguidata','metricshandle','metricshandles','metricswindowhandles','metricswindowhandle'}
        v = ieSessionGet('vcimagefigure');
        if ~isempty(v), val = guihandles(v); end
        
        % Figure numbers of the various windows.  I am not sure these are
        % properly updated, but I think so.
    case {'mainfigure','mainfigures','mainwindow'}
        if checkfields(vcSESSION,'GUI','vcMainWindow')
            val = vcSESSION.GUI.vcMainWindow.hObject;
        end
    case {'scenefigure','sceneimagefigure','sceneimagefigures','scenewindow'}
        if checkfields(vcSESSION,'GUI','vcSceneWindow')
            val = vcSESSION.GUI.vcSceneWindow.hObject;
        end
    case {'oifigure','opticalimagefigure','oifigures','opticalimagefigures','oiwindow'}
        if checkfields(vcSESSION,'GUI','vcOptImgWindow')
            val = vcSESSION.GUI.vcOptImgWindow.hObject;
        end
    case {'sensorfigure','isafigure','sensorfigures','isafigures','sensorwindow','isawindow'}
        if checkfields(vcSESSION,'GUI','vcSensImgWindow')
            val = vcSESSION.GUI.vcSensImgWindow.hObject;
        end
        
    case {'vcimagefigure','vcimagefigures','vcimagewindow'}
        if checkfields(vcSESSION,'GUI','vcImageWindow')
            val = vcSESSION.GUI.vcImageWindow.hObject;
        end
        
    case {'metricsfigure','metricswindow','metricsfigures'}
        if checkfields(vcSESSION,'GUI','metricsWindow')
            val = vcSESSION.GUI.metricsWindow.hObject;
        end
        
        % I think all this custom routine management from the GUI has not
        % been used and is not likely to be helpful.  It should go away.
        % We should allow people to use scripts to set custom routines, but
        % they should not manage them from the GUI. - BW 2010.
    case {'custom','customall','customstructure'}
        val = vcSESSION.CUSTOM;
        % These are cell arrays of user-defined routines that implement
        % these various types of operations.
    case {'demosaiclist','customdesmoaiclist','customdemosaic'}
        if checkfields(vcSESSION,'CUSTOM','demosaic')
            val = vcSESSION.CUSTOM.demosaic;
        end
    case {'balancelist','colorbalancelist','customcolorbalancelist','customcolorbalance'}
        if checkfields(vcSESSION,'CUSTOM','colorBalance')
            val = vcSESSION.CUSTOM.colorBalance;
        end
    case {'conversionlist','colorconversionlist','customcolorconversionlist','customcolorconversion'}
        if checkfields(vcSESSION,'CUSTOM','colorConversion')
            val = vcSESSION.CUSTOM.colorConversion;
        end
    case {'renderlist','processinglist','processingmethods'}
        % These routines are a complete processing chain that replace the
        % entire call to vcimageCompute
        if checkfields(vcSESSION,'CUSTOM','render')
            val = vcSESSION.CUSTOM.render;
        end
    case {'oicomputelist'}
        % These routines replace the standard oiCompute call.  They
        % customize the signal processing chain from the optical image to
        % the OI data.
        if checkfields(vcSESSION,'CUSTOM','oicompute')
            val = vcSESSION.CUSTOM.oicompute;
        end
    case {'sensorcomputelist'}
        % These routines replace the standard sensorCompute call.  They
        % customize the signal processing chain from the optical image to
        % the ISA data.
        if checkfields(vcSESSION,'CUSTOM','sensorcompute')
            val = vcSESSION.CUSTOM.sensorcompute;
        end
    case {'edgealgorithmlist'}
        % These routines replace the standard sensorCompute call.  They
        % customize the signal processing chain from the optical image to
        % the ISA data.
        if checkfields(vcSESSION,'CUSTOM','edgeAlgorithm')
            val = vcSESSION.CUSTOM.edgeAlgorithm;
        end
        
    otherwise
        error('Unknown parameter')
end
end