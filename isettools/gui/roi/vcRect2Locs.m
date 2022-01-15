function roiLocs = vcRect2Locs(rect)
% Obsolete. Use ieRoi2Locs. Convert rect from ISET window into RoI locs.
%
% Syntax:
%   roiLocs = vcRect2Locs(rect)
%
% Description:
%    OBSOLETE. Use this function to convert a rect from an ISET window into
%    region of interest locations.
%
%    The rect is usually selected using the Matlab grphical user interface.
%    The rect values are:
%        [colMin rowMin (width-1) (height-1)].
%    (Note that col is the x dimension and row is the y dimension).
%
%    Usually we call the routine vcROISelect directly, which then calls
%    this routine:
%      vci = vcGetObject('vcimage');
%      [roiLocs roiRect] = vcROISelect(vci);
%
%    The code below contains examples of function usage. To access, type
%    'edit vcRect2Locs.m' into the Command Window.
%
% Inputs:
%    rect    - Vector. A 1x4 vector of the format described above.
%
% Outputs:
%    roiLocs - Matrix. An Nx2 matric of (row, col) values.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    Update: Use ieRoi2Locs instead of this function.
%    vcROISelect(), vcLocs2Rect
%

% History:
%    xx/xx/04       (c) Imageval, 2004
%    05/09/18  jnm  Formatting.

% Examples:
%{
    % ETTBSkip - skipping broken example
    rect = round(getrect(ieSessionGet('vcimagewindow')));
    roiLocs = vcRect2Locs(rect);
%}

% The rect entries are (colMin, rowMin, colWidth - 1, rowWidth - 1)
% The number of data values for the dimensions follow colMax - colMin + 1
% and likewise for the rows.

cmin = rect(1);
cmax = rect(1) + rect(3);
rmin = rect(2);
rmax = rect(2) + rect(4);

[c, r] = meshgrid(cmin:cmax, rmin:rmax);
roiLocs = [r(:), c(:)];

return;