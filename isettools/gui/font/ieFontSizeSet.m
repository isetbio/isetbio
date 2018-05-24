function fSize = ieFontSizeSet(fig, fSize)
% Set the font size of all the text in the window objects
%
% Syntax:
%   fSize = ieFontSizeSet(fig, fSize);
%
% Description:
%    Set the font size for all of the text in the window object(s).
%
%    The font size is set to all the text in the window.  The first textbox
%    in the window is the one that is assigned the fSize passed in here.
%    When the default size for the text is a little bigger or smaller in
%    the different boxes, the relative amount is preserved.
%
%    (This is part of replacing ieFontChangeSize)
%
%    There are examples contained in the code. To access, type 'edit
%    ieFontSizeSet.m' into the Command Window.
%
% Inputs:
%    fig   - Handle. A handle to the figure of, say, the scene or oi window
%    fSize - Integer. The font size you want for this window
%            If fSize is 0, we are simply refreshing the window
%            if fSize is missing, we bring up a window and ask the user
%            Otherwise, we uset the actual fSize value
%
% Outputs:
%    fSize - The font size that was set
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - Why are size preferences and min/max sizes unique (and
%      contradictory) for each function that touches on font size?]
%

% History:
%    xx/xx/15       Copyright Imageval Consulting, LLC, 2015
%    02/28/18  jnm  Formatting

% Example:
%{
   s = sceneCreate;
    ieAddObject(s);
    sceneWindow;
   fig = ieSessionGet('scene window');
   ieFontSizeSet(fig, 14);
%}

%% Set up parameters

if notDefined('fig'), error('Figure required.'); end

% Pull out the current font size preference
isetP = getpref('ISETBIO');
if checkfields(isetP, 'fontSize')
    prefSize = isetP.fontSize;
else
    prefSize = 12;  % Default preference
end

if notDefined('fSize')
    % fSize is empty or missing, so ask the user
    fSize = ieReadNumber('Enter font size (7-25): ', prefSize, ' %.0f');
    if isempty(fSize), return; end
elseif fSize == 0
    % Refresh condition. Use the ISET pref 
    fSize = prefSize;
end

% Clip to range
minSize = 7;
maxSize = 25;
fSize = ieClip(fSize, minSize, maxSize);

%% Apply the new change in the font size to the window. 

% If there is no window, we just update the preference
% Get all the children of the figure.
t = allchild(fig);

% Change the text displays
tHandles = findall(t, 'Style', 'Text');
setFontSize(tHandles, fSize);

% Change the popupmenu font sizes.
tHandles = findall(t, 'Style', 'popupmenu');
setFontSize(tHandles, fSize);

% Change the popupmenu font sizes.
tHandles = findall(t, 'Style', 'edit');
setFontSize(tHandles, fSize);

% Change the radiobutton font sizes.
tHandles = findall(t, 'Style', 'radiobutton');
setFontSize(tHandles, fSize);

% Change the pushbutton font sizes.
tHandles = findall(t, 'Style', 'pushbutton');
setFontSize(tHandles, fSize);

setpref('ISETBIO', 'fontSize', fSize);

end

function setFontSize(tHandles, fSize)
% Set the size of the font for all of this type of handle
%
% Syntax:
%   setFontSize(tHandlse, fSize)
%
% Description:
%    Set the size of the font for all handles of this type.
%
% Inputs:
%    tHandles - Handle(s). The handles to modify.
%    fSize    - Integer. The desired font size.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Current font sizes
curSize = get(tHandles, 'FontSize');

% Some of the handles are odd objects that we don't need.
if iscell(curSize)
    lst = (cell2mat(curSize) > 0);
    curSize = curSize(lst);
    tHandles = tHandles(lst);   
end

if isempty(curSize)
    % No fonts to change.
    return;
elseif length(curSize) == 1
        % Single font size case
        set(tHandles, 'FontSize', fSize);
else
    % Set as if first size is the base size and everything else is offset
    for ii = 1:length(curSize)
        if ii == 1
            offset = 0;
        else
            offset = curSize{ii} - curSize{1};
        end
        thisSize = fSize + offset;
        set(tHandles, 'FontSize', thisSize);
    end
end

end
