function roiLocs = ieRoi2Locs(rect)
% Should be called ieRect2Locs, converts rectangles to roi.
%
% Syntax:
%   roiLocs = ieRoi2Locs(rect)
%
% Description:
%    Convert rectangle coordinates into roi locations. The rect format is:
%
%       (column, row, width, height).
%
%    The first data point is drawn from the position [col, row].
%
%    The roiLocs are a N x 2. The roiLocs have rows in the first column and
%    columns in the 2nd.
%
%    For unfortunate historical reasons, which could be fought here, the
%    spatial size of the returned data are width + 1 and height + 1. Thus,
%    [col, row, 1, 1] returns four positions. [col, row, 0, 0] returns 1
%    position. Blame it on C and Fortran. Or me.
%
% Inputs:
%    rect    - Vector. 1 x 4 vector of the rectangular locations. Follows
%              the format [column, row, width, height].
%
% Outputs:
%    roiLocs - Vector. N x 2 vector of ROI Locations. First column contains
%              rows, and the second contains columns.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    vcROISelect, vcLocs2Rect
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/02/18  jnm  Formatting

% The rect entries are  The number of data
% values are colMax - colMin +1 and similarly for the row
cmin = rect(1);
cmax = rect(1) + rect(3);
rmin = rect(2);
rmax = rect(2) + rect(4);

[c, r] = meshgrid(cmin:cmax, rmin:rmax);
roiLocs = [r(:), c(:)];

return;
