function [inputRow, inputCol] = inputPositions(rgcMosaic, row, col)
% Calculate the receptive field position of one RGC cell in input frame
%
% Syntax:
%   [inputRow, inputCol] = inputPositions(obj, row, col)
%
% Description:
%    The calling routine determines the spacing of the input bipolars,
%    which can vary by mosaic type. So we take that as an input argument.
%
%    The important machinery of coordinating the spatial coordinate frames
%    between layers is managed by finding the position of the cells in
%    different layers with respect to
%
%     (a) the locations (in microns) on the cone mosaic, and
%     (b) the location with respect to sample positions of cells in the
%         input layers
%
%    This is a good routine for using in visualization for demonstrating
%    spatial receptive fields in different ways.
%
%    The length of the inputRow and inputCol has to be the same as the size
%    of the receptive field row and col
%
% Inputs:
%    rgcMosaic - String. The rgc mosaic type. Select from the type list.
%    row       - Numeric. The row position of the RGC in the array.
%    col       - Numeric. The column position of the RGC in the array.
%    bipolarsPerMicron
%              - Numeric. The input layer (bipolars) has a spatial density
%                This can vary with the input cell type. Passed in from the
%                rgcSpaceDot function.
%
% Outputs:
%    inputRow  - Array. The inputRow and inputCol define the portion of the
%                (bipolar) input used to compute the inner product of the
%                spatial receptive field. These are vectors and the whole
%                (row, col) space can be built by meshgrid(inputRow,
%                inputCol) or something like that.
%    inputCol  - Array. See inputRow above.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   rgcSpaceDot
%

% History:
%    05/XX/16  JRG  (c) ISETBIO Team
%    06/13/19  JNM  Documentation pass

% The RGC center location is specified in bipolar samples.
rgcCenter = rgcMosaic.cellLocation(row, col, :);

% The sRFcenter{} is the set of weights that will be applied to the input
% layer. The positions of the weights are in the sample coordinates of the
% input. We calculate the midpoint of RF by taking half of the row size.
sRFMidPointRow = (1 / 2) * size(rgcMosaic.sRFcenter{row, col}, 1);

% The row and samples where read the input for this cell
rowStart = (rgcCenter(1) - sRFMidPointRow);
% rowEnd = (rgcCenter(1) + sRFMidPointRow);

% Repeat for the column dimension.
sRFMidPointCol = (1 / 2) * size(rgcMosaic.sRFcenter{row, col}, 2);
colStart = (rgcCenter(2) - sRFMidPointCol);
% colEnd = (rgcCenter(2) + sRFMidPointCol);

% The offset parameter is the distance in bipolar samples to the edge
% of the input.
offset = rgcMosaic.cellLocation(1, 1, :);

nRow = (0:size(rgcMosaic.sRFcenter{row, col}, 1) - 1);
nCol = (0:size(rgcMosaic.sRFcenter{row, col}, 2) - 1);
inputRow = ceil(nRow + rowStart - offset(1));
inputCol = ceil(nCol + colStart - offset(2));

% Older code.
% Check to make sure that the inputRow/Col match the size of the receptive
% field. If we never get this error, then we will delete the check.
% if (length(inputRow) ~= size(rgcMosaic.sRFcenter{row, col}, 1)) || ...
%         (length(inputCol)~=size(rgcMosaic.sRFcenter{row, col}, 2))
%     error('Dimension mismatch of input and rf size');
% end

end