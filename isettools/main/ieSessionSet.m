function ieSessionSet(param, val, varargin)
% Set vcSESSION parameters.
%
% Syntax:
%   ieSessionSet(param, val, varargin);
%
% Description:
%    While vcSESSION parameters are generally set by this routine. There
%    remain some places, however, where vcSESSION is touched directly. This
%    will change over time. At some point in time, vcSESSION is likely to
%    be hidden in the main window, and not be global.
%
%    The function contains examples of usage. To access, type 'edit
%    ieSessionSet.m' into the command window.
%
% Inputs:
%    param - String. The parameter name you wish to assign the value to.
%            Some possible options are listed below, with the type for each
%            one's value (val) included:
%       General
%           {'version'}: Numeric. The vcSession version.
%           {'sessionname'}: String. The vcSession name.
%           {'sessiondir'}: String. The vcSession working directory.
%           {'inithelp'}: Boolean. Whether or not to initialize help.
%       Matlab setpref variables
%           {'detlafontsize'}: Numeric. This value determines whether we
%                              change the font size in every window by this
%                              increment, calling ieFontChangeSize when the
%                              window is opened.
%           {'waitbar'}: String. Show waitbars (or not) during certain
%                        long computations. Options are 'on', and 'off'.
%           {'gpu computing'}: Boolean. Whether/not to use GPU computing.
%       Figure handles
%           {'main window'}: Handle. Store handles for main window
%           {'scene window'}: Handle. Store handles for scene window
%           {'oi window'}: Handle. Store handles of optical image window
%           {'sensor window'}: Handle. Store handles for sensor window
%           {'cone mosaic window'}: Handle. Store handle for the cone
%                                   mosaic window.
%           {'graphwin val'}: Numeric. The number for graphics window.
%           {'graphwin handle'}: Handle. Not currently used, in future will
%                                be as named.
%           {'graphwin figure'}: Object. hObject for graphics window. Why
%                                is this not handle?
%    val   - VARIES. The value to assign to the specified parameter. Find
%            the specific parameter above to see the required value's type.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    Needs to be added.
%
% Notes:
%    * PROGRAMMING TODO:
%      It seems that not all of the session sets are handled through this
%      call yet. Must find more. Look at the routine vcSetFigureHandles for
%      some clues.
%      Here's a clue:  vcReplaceObject. Maybe vcAddandSelectObject ...
%      Rather than vcGetObject, we should be using ieSessionGet('scene') or
%      ieSessionGet('scene', 3);  Sigh. For historical reasons, the
%      vcSESSION was not properly protected. Hence, there are still way too
%      many vcSESSION. calls in the sub-routines.
%
% See Also:
%   ieSessionGet
%

% History:
%    XX/XX/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/28/19  JNM  Documentation pass

% Examples:
%{
    % ETTBSkip - skipping broken examples
    %
    % This needs code to define newALg before it could
    % possibly work.
    ieSessionSet('addrender', newAlg);
    ieSessionSet('main window', hObject, eventdata, handles);
%}

global vcSESSION

if notDefined('param'), error('You must specify a parameter.'); end
if ~exist('val', 'var'), error('You must specify a value.'); end

param = ieParamFormat(param);
switch param
    case {'version'}
        vcSESSION.VERSION = val;
    case {'name', 'sessionname'}
        vcSESSION.NAME = val;
    case {'dir', 'sessiondir'}
        vcSESSION.DIR = val;
    case {'help', 'inithelp'}
        % Default for help is true, if the initHelp has not been set.
        if checkfields(vcSESSION, 'initHelp'), vcSESSION.initHelp = val;
        else, vcSESSION.initHelp = 1;
        end
    % Matlab setpref values
    case {'detlafontsize', 'fontincrement', 'increasefontsize', ...
            'fontdelta', 'deltafont'}
        % This value determines whether we change the font size in every
        % window by this increment, calling ieFontChangeSize when the
        % window is opened.
        setpref('ISET', 'fontDelta', val);
    case {'waitbar'}
        % 0 means off, 1 means on
        if ischar(val)
            switch val
                case 'on', val = 1;
                case 'off', val = 0;
            end
        end
        setpref('ISET', 'waitbar', val);
        % Because getpref is slow, we also attach it to the session. Then
        % looping and checking doesn't cost us much time.
        vcSESSION.GUI.waitbar = val;
    case {'gpu', 'gpucompute', 'gpucomputing'}
        vcSESSION.GPUCOMPUTE = val;
    case {'imagesizethreshold'}
        vcSESSION.imagesizethreshold = val;
    % Set window information at startup
    case {'mainwindow'}
        if length(varargin) < 2
            error('main window requires hObject, eventdata, handles');
        end
        vcSESSION.GUI.vcMainWindow.hObject = val;
        vcSESSION.GUI.vcMainWindow.eventdata = varargin{1};
        vcSESSION.GUI.vcMainWindow.handles = varargin{2};
    case {'scenewindow'}
        if length(varargin) < 2
            error('scene window requires hObject, eventdata, handles');
        end
        vcSESSION.GUI.vcSceneWindow.hObject = val;
        vcSESSION.GUI.vcSceneWindow.eventdata = varargin{1};
        vcSESSION.GUI.vcSceneWindow.handles = varargin{2};
    case {'oiwindow'}
        if length(varargin) < 2
            error(strcat('optical image window requires hObject, ', ...
                'eventdata, handles'));
        end
        vcSESSION.GUI.vcOptImgWindow.hObject = val;
        vcSESSION.GUI.vcOptImgWindow.eventdata = varargin{1};
        vcSESSION.GUI.vcOptImgWindow.handles = varargin{2};
    case {'sensorwindow'}
        if length(varargin) < 2
            error('sensor window requires hObject, eventdata, handles');
        end
        vcSESSION.GUI.vcSensImgWindow.hObject = val;
        vcSESSION.GUI.vcSensImgWindow.eventdata = varargin{1};
        vcSESSION.GUI.vcSensImgWindow.handles = varargin{2};
    case {'conemosaicwindow'}
        if length(varargin) < 2
            error(strcat('cone mosaic window requires hObject, ', ...
                'eventdata, handles'));
        end
        vcSESSION.GUI.vcConeImgWindow.hObject = val;
        vcSESSION.GUI.vcConeImgWindow.eventdata = varargin{1};
        vcSESSION.GUI.vcConeImgWindow.handles = varargin{2};
    % This graphics window stuff is a mess
    case {'graphwinstructure', 'graphwinval'}
        vcSESSION.GRAPHWIN = val;
    case {'graphwinhandle'}
        % At present we don't add any objects with handles. So this is
        % empty. But we might some day.
        vcSESSION.GRAPHWIN.handle = val;
    case {'graphwinfigure'}
        % This is just the figure number, usually.
        vcSESSION.GRAPHWIN.hObject = val;
    otherwise
        error('Unknown parameter')
end
