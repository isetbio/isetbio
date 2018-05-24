function gMat = ieReadSmallMatrix(mSize, defMatrix, fmt, prompt, ...
    yxPosition, outVarName, cList)
% A GUI to read a small matrix in a structured way
%
% Syntax:
%   m = ieReadSmallMatrix(mSize)
%
% Description:
%    We create a figure, gMat, that contains edit boxes where the user can
%    enter matrix data.
%
%    The boxes can be slightly colored to indicate, say, the color of the
%    filter or some useful information for the user.
%
%    The code below contains examples of function usage. To access, type
%    'edit ieReadSmallMatrix.m' into the Command Window.
%
% Inputs:
%    mSize      - (Optional) Vector. 1 x 2 vector containing the matrix
%                 size. Default [3 3].
%    defMatrix  - (Optional) Matrix. Default Matrix variable. Default [].
%    fmt        - (Optional) String. Matrix display format. Default '%.2f'.
%    prompt     - (Optional) String. Usesr prompt. Default 'Enter matrix'.
%    yxPosition - (Optional) Vector. 1 x 2 vector containing the starting
%                 yx position. Default [100 500].
%    outVarName - (Optional) String. Output variable name. Default string
%                 is 'matrixData'.
%    cList      - (Optional) Array. The color list array. Default [].
%
% Outputs:
%    gMat       - The generated GUI for the matrix
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: A deficiency of this routine, compared to ieReadMatrix, is
%      that we can't just copy and paste into the window. You have to type
%      numbers (you can <TAB> through the boxes). It would be nice to add
%      the function where you can paste something into the entries when
%      there is a matrix in the clipboard.
%    * TODO: Assign someone to fix example.
%

% Examples:
%{
    % ETTBSkip - example not functioning, as sensor is not defined. Skip.
    mSize = [2, 3];
    defMatrix = [];
    fmt = [];
    prompt = 'Durations (ms)';
    yxPosition = [100, 500];
    saturation = 0.3;
    filterRGB = sensorFilterRGB(sensor, saturation);
    ieReadSmallMatrix(...
        mSize, defMatrix, fmt, prompt, yxPosition, 'myData', filterRGB);
    defMatrix = myData
%}

%%
if notDefined('mSize'), mSize = [3, 3]; end
if notDefined('defMatrix'), defMatrix = zeros(mSize); end
if notDefined('fmt'), fmt = '%.2f'; end
if notDefined('prompt'), prompt ='Enter matrix'; end
if notDefined('yxPosition'), yxPosition = [100, 500]; end
if notDefined('outVarName'), outVarName = 'matrixData'; end
if notDefined('cList'), cList = []; end

%%
% Create the figure and main window

% See if a previous session exists; if so, close it
tmp = findobj('Tag', 'cfaExposure');
if ~isempty(tmp), close(tmp); end

nRows = mSize(1);
nCols = mSize(2);
gMat.fig = figure;
gMat.mSize = mSize;

if notDefined('position')
    % figPosition = get(gMat.fig, 'Position');
    figWidth = nCols * 40 + 100; % pixels
    figHeight = nRows * 30 + 100; % pixels
    position = [yxPosition(1), yxPosition(2), figWidth, figHeight];
end

set(gMat.fig, ...
    'Tag', 'CFAexposure', ...
    'NumberTitle', 'Off', ...
    'Resize', 'on', ...
    'Position', position, ...
    'Name', prompt, ...
    'Menubar', 'None'...
    );

% Panel position within the figure.
gMat.panel = uipanel('Parent', gMat.fig, ...
    'Position', [0.01 0.21 0.98 0.78], ...
    'BackgroundColor', get(gcf, 'Color')...
    );

% Create the boxes
% Set dimensions (relative to the dimensions of the main panel) of the
% editboxes that are used to enter matrix elements
p1 = 0.9;
p2 = 0.1;
editWidth = p1 / nCols;
editHeight = p1 / nRows;
editSpacingX = p2 / (nCols + 1);
editSpacingY = p2 / (nRows + 1);

editboxCallBack = {@editboxCallback};
% Place the boxes in the window. Note the odd ordering. We want to draw the
% boxes from the upper left down to lower right so the <TAB> will work OK.
for currentRow = 1:nRows
    for currentCol = 1:nCols
        if isempty(cList)
            thisColor = [1 1 1];
        else
            thisColor = squeeze(cList(currentRow, currentCol, :));
        end

        xPos = editSpacingX + (currentCol - 1)*(editWidth + editSpacingX);
        yPos = editSpacingY + (nRows - currentRow) * ...
            (editHeight + editSpacingY);

        pointIndex = sub2ind([nRows, nCols], currentRow, currentCol);
        str = sprintf(fmt, defMatrix(currentRow, currentCol));
        gMat.exps(pointIndex).val = uicontrol(...
            'Parent', gMat.panel, ...
            'Style', 'edit', ...
            'Units', 'normalized', ...
            'Position', [0.1 0.4 0.8 0.5], ...
            'BackgroundColor', thisColor, ...
            'String', str, ...
            'Position', [xPos yPos editWidth editHeight], ...
            'Callback', editboxCallBack);

    end
end

doneCallBack = {@doneCallback, gMat, outVarName};

% Done button
gMat.buttonDone = uicontrol('Parent', gMat.fig, ...
    'Style', 'pushbutton', ...
    'String', 'Done', ...
    'Units', 'Normalized', ...
    'Position', [0.30 0.03 0.4 0.14], ...
    'BackgroundColor', get(0, 'defaultUicontrolBackgroundColor'), ...
    'Callback', doneCallBack, ...
    'KeyPressFcn', doneCallBack...
    );

% Wait for a user response
set(gMat.fig, 'pointer', 'watch');
uiwait(gMat.fig);

return;

function editboxCallback(~, ~, varargin)
% In theory shouldn't do anything....
%
% Syntax:
%    editboxCallback(~, ~, [varargin])
%
% Description:
%    Edit call back shouldn't really do anything, though it might adjust
%    the scale if we get really frisky.
%
% Inputs:
%    None required.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
return

function doneCallback(~, ~, varargin)
% Read all edit boxes when DONE button pressed.
%
% Syntax:
%   doneCallback(~, ~, [varargin])
%
% Description:
%    The DONE button from the window is the object. Reads values from all
%    the edit boxes.
%
% Inputs:
%    None required.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
gMat = varargin{1};
outVarName = varargin{2};
mSize = gMat.mSize;
nRows = mSize(1);
nCols = mSize(2);

% Loop through the edit boxes and find the variables
%
matrixData = zeros(size(mSize));
for ii = 1:nRows
    for jj = 1:nCols
        editBoxIndex = sub2ind([nRows, nCols], ii, jj);
        hEditBox = gMat.exps(editBoxIndex).val;
        matrixData(ii, jj) = str2double(get(hEditBox, 'String'));
    end
end

% Assign the matrix into the base work space
assignin('base', outVarName, matrixData);

% Resume and close the window
uiresume(gMat.fig);
set(gMat.fig, 'pointer', 'arrow');
close(gMat.fig);

return
