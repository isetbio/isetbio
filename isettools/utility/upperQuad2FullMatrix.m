function fullmatrix = upperQuad2FullMatrix(upperRight, nRows, nCols)
% Duplicates the upper right quad of a matrix in the other quads
%
% Syntax:
%   fullmatrix = upperQuad2FullMatrix(upperRight, nRows, nCols)
%
% Description:
%    For a full matrix (nRows x nCols), suppose upperRight is the upper
%    right quadrant of a matrix. 
%
%    This routine duplicates the upper right quad in the other quads, 
%    mirroring the data. 
%
%    If there is an odd number of output cols, the middle data are always
%    attached to the quadrants on the right side of the data. If there is
%    an odd number of rows, the middle data are always attached to the
%    upper quadrants.
%
%    Examples in code.
%
% Inputs:
%    upperRight - Smaller matrix that will be duplicated across the 'full'
%                 matrix via replication and mirroring.
%    nRows      - Number of rows in the matrix
%    nCols      - Number of columns in the matrix
%
% Outputs:
%    fullmatrix - The expanded matrix
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - nRows and nCols appear to be 'suggestions' rather than
%      strictly enforced?]

% History:
%    xx/xx/05       (c) ImagEval Consultants, LLC, 2005.
%    11/17/17  jnm  Formatting
%    01/17/18  jnm  Formatting update to match Wiki.

% Example:
%{
    upperRight = [1, 2, 3; 4, 5, 6; 7, 8, 9];
    nRows = 5;
    nCols = 6;
    fullmatrix = upperQuad2FullMatrix(upperRight, nRows, nCols)
%}

[r, c] = size(upperRight);
if isodd(nCols)
    upperLeft = fliplr(upperRight(:, 2:c));
else
    upperLeft = fliplr(upperRight);
end

if isodd(nRows)
    lowerRight =  flipud(upperRight(1:(r - 1), :));
else
    lowerRight =  flipud(upperRight);
end

if isodd(nRows)
    lowerLeft = flipud(upperLeft(1:(r - 1), :));
else
    lowerLeft = flipud(upperLeft);
end

fullmatrix = [upperLeft, upperRight; lowerLeft, lowerRight];

end
