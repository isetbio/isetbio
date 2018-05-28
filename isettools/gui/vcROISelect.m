function [roiLocs, rect] = vcROISelect(obj, objFig)
% Select a region of interest (ROI) from an image and calculate locations
%
% Syntax:
%   [roiLocs, rect] = vcROISelect(obj, [objFig])
%
% Description:
%    The row and col locations of the region of interest (ROI) are returned
%    in the Nx2 matrix, roiLocs.
%
%    If requested, the selected rectangle (rect) determining the region of
%    interest, [colmin, rowmin, width, height], is also returned.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcROISelect.m' into the Command Window.
%
% Inputs:
%    obj     - Object. The object in question.
%    objFig  - (Optional) Handle. Handle to object figure.
%
% Outputs:
%    roiLocs - Matrix. Matrix of RoI Locations. Nx2 Matrix of locations.
%    rect    - Vector. Vector containing RoI in [x, y, width, height] where
%              x and y are the column minimum and row minimum respectively.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: See proposal for ieOpenWindow below. We should also add
%      ieRoiSelect to plan for deprecation of the vcXXX routines.
%
% See Also:
%    ieRoi2Locs
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/11/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - Skipping broken example that when fixed will require user
    % input to function as desired.
    vci = vcGetObject('VCIMAGE');
    [roiLocs, rect] = vcROISelect(vci);
    iData = vcGetROIData(vci, roiLocs, 'results');
%}

if notDefined('obj'), error('You must define an object'); end
if notDefined('objFig')
    objFig = vcGetFigure(obj);
    if isempty(objFig)
        % We should add ieAddAndSelect()
        vcAddAndSelectObject(obj);
        % Should become ieOpenWindow(obj)
        switch obj.type
            case 'scene'
                objFig = sceneWindow;
            case 'opticalimage'
                objFig = oiWindow;
            case 'sensor'
                objFig = sensorWindow;
            case 'vcimage'
                objFig = ipWindow;
            otherwise
                error('Unknown obj type %s\n', obj.type);
        end
    end
end

% Select points.
hndl = guihandles(objFig);
msg = sprintf('Drag to select a region.');
ieInWindowMessage(msg, hndl);

% Select an ROI graphically. Calculate the row and col locations.
% figure(objFig);
rect = round(getrect(objFig));
ieInWindowMessage('', hndl);

% If the user double clicks without selecting a rectangle, we treat the
% response as a single point. We do this by making the size 1, 1.
if rect(3) == 0 && rect(4) == 0
    rect(3) = 1;
    rect(4) = 1;
end

% Transform the rectangle into ROI locations
roiLocs = ieRoi2Locs(rect);

return;
