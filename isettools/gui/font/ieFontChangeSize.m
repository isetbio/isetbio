function dSize = ieFontChangeSize(fig, dSize)
% Change and apply the font size preference used in an ISET window
%
% Syntax:
%   dSize = ieFontChangeSize(fig, [dSize])
%
% Description:
%    The font size preference information is stored using the Matlab
%    setpref/getpref mechanism. Hence, the size is remembered across
%    sessions (or installs of ISET).
%
%    This routine is used in many ISET windows and also on startup.
%
%    The Matlab preferences are checked when a window opens. The
%    preferences are applied to the window. Hence, the preferences are
%    transferred between windows.
%
%    The change size is limited to   -6 and 6.
%    The total change is limited to -12 and 12.
%    The smallest point size is limited to 6pt on windows & 10pt on linux.
%    There is no upper limit on the max point size.
%
%    There are examples contained in the code. To access the examples, type
%    'edit ieFontChangeSize.m' into the Command Window.
%
% Inputs:
%    fig   - Handle. The figure handle
%    dSize - (Optional) Integer. The amount to change the display size for
%            text. Limits are -6 and 6. Default is empty (query user).
%
% Outputs:
%    dSize - Integer. The amount to change the display text size by.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: Assign someone to flesh out example.
%
% See Also:
%    ieFontInit, ieSessionGet/Set
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    02/28/18  jnm  Formatting

% Examples:
%{
    % Example needs to be fleshed out. Currently if called on an existing
    % window changes size of all text to system minimum.
   ieFontChangeSize(ieSessionGet('oifigure'), ...
        ieSessionGet('increasefontsize'));
%}

% This routine either receives a dSize or it gets it from the user. When
% dSize is changed here by the user, we also change the dSize in the Matlab
% preferences.
if notDefined('dSize')
    % Find the increase or decrease in font size from the user. This is
    % always relative to the current deltaFont value.
    dSize = ieReadNumber('Enter font size change (-6, 6)', 2, ' %.0f');
    if isempty(dSize), return; end

    % Keep the font size change within a reasonable range
    dSize = ieClip(dSize, -6, 6);

    % This is a new setting, so update the preferred change in font size
    % from the Matlab default. 
    oldVal = ieSessionGet('deltaFont');

    % But make sure that the change is within a reasonable range.
    newVal = ieClip(oldVal + dSize, -12, 12);
    ieSessionSet('deltaFont', newVal);
    % fprintf('Adjusting delta font to %d\n', newVal)
else
    % If we didn't read the dSize, there is no reason to change it.
end

%% Apply the new change in the font size to the text. 

% Get all the children of the figure.
t = allchild(fig);

% Change the text displays
tHandles = findall(t, 'Style', 'Text');
changeFontSize(tHandles);

% Change the popupmenu font sizes.
tHandles = findall(t, 'Style', 'popupmenu');
changeFontSize(tHandles);

% Change the popupmenu font sizes.
tHandles = findall(t, 'Style', 'edit');
changeFontSize(tHandles);

% Change the radiobutton font sizes.
tHandles = findall(t, 'Style', 'radiobutton');
changeFontSize(tHandles);

% Change the pushbutton font sizes.
tHandles = findall(t, 'Style', 'pushbutton');
changeFontSize(tHandles);

return;

function changeFontSize(tHandles)
% Algorithm for changing font size. Never let font size get smaller than 6.
%
% Syntax:
%   changeFontSize(tHandles)
%
% Description:
%    Algorithm for changing the font size. Never allow to go below 6.
%
%    We have a baseSize for each system. We have a current size for the
%    fonts, deltaFont, stored in Matlab preferences. This delta defined the
%    change from the default baseSize.
%
%    We find the difference between the currentSize and the baseSize and we
%    adjust the current size so that the difference equals deltaFont.
%
% Inputs:
%    tHandles - Handle for the figure
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: XXX - Why should we never let the font get smaller than a 6?]
%    * [Note: JNM - We say forced minimum is 6, but there are multiple
%      exceptions. Why have these?]
%
if ispc
    baseSize = 6;  % Windows
elseif isunix
    baseSize = 10;  % Linux?
else
    baseSize = 8;  % Maybe apple?
end

curSize = get(tHandles, 'FontSize');
if isempty(curSize), return;
else
    deltaFont = ieSessionGet('deltaFont');  % Total change from baseline
    if length(curSize) == 1
        desiredSize = baseSize + deltaFont;
        set(tHandles, 'FontSize', max(desiredSize, baseSize));
    else
       currentDelta = max(curSize{1} - baseSize, 0);
        for ii = 1:length(curSize) 
            % These are the base sizes of each object
            desiredSize = (curSize{ii} - currentDelta) + deltaFont;
            if isempty(desiredSize)
                set(tHandles(ii), 'FontSize', baseSize);
            else
                set(tHandles(ii), 'FontSize', max(desiredSize, baseSize));
            end
        end
    end
end

return
