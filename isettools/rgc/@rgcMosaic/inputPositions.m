function [inputRow, inputCol] = inputPositions(rgcMosaic,row,col)
% Calculate the receptive field position of one RGC cell in input frame
%
%  [inputRow, inputCol] = inputPositions(obj,row,col)
%
% Inputs:
%    rgcMosaic - any type
%    row,col - position of the RGC in the array
%    bipolarsPerMicron - The input layer (bipolars) has a spatial density
%          This can vary with the input cell type.   Passed in by
%          rgcSpaceDot.m
%
% Return
%   inputRow and inputRow define the portion of the (bipolar) input used to
%     compute the inner product of the spatial receptive field.  These are
%     vectors and the whole (row,col) space can be built by
%     meshgrid(inputRow,inputCol) or something like that
%
% The calling routine determines the spacing of the input bipolars, which
% can vary by mosaic type.  So we take that as an input argument.
%
% The important machinery of coordinating the spatial coordinate frames
% between layers is managed by finding the position of the cells in
% different layers with respect to 
%
%  (a) the locations (in microns) on the cone mosaic, and 
%  (b) the location with respect to sample positions of cells in the input
%  layers
%
% See also:  rgcSpaceDot
%
% Examples:
%    Good routine for using in visualization for demonstrating spatial
%    receptive fields in different ways.
%
% 5/2016 JRG (c) ISETBIO Team

% The length of the inputRow and inputCol has to be the same as the size of
% the receptive field row and col

% The RGC center location is specified in bipolar samples.
rgcCenter = rgcMosaic.cellLocation(row,col,:);

% The sRFcenter{} is the set of weights that will be applied to the input
% layer.  The positions of the weights are in the sample coordinates of the
% input. We calculate the midpoint of RF by taking half of the row size.
sRFMidPointRow = (1/2)*size(rgcMosaic.sRFcenter{row,col},1);

% The row and samples where read the input for this cell
rowStart = (rgcCenter(1) - sRFMidPointRow);
rowEnd   = (rgcCenter(1) + sRFMidPointRow);

% Repeat for the column dimension.
sRFMidPointCol = (1/2)*size(rgcMosaic.sRFcenter{row,col},2);
colStart = (rgcCenter(2) - sRFMidPointCol);
colEnd   = (rgcCenter(2) + sRFMidPointCol);

% The offset parameter is the distance in bipolar samples to the edge
% of the input. 
offset =rgcMosaic.cellLocation(1,1,:);

% Add the eps0 offset to each position and apply ceil and floor to ensure
% inputRow and inputCol are equal to size(rgcMosaic.sRFcenter{row,col})
eps0 = .0001;
inputRow =  ceil(rowStart - offset(1) + eps0):floor((rowEnd) - offset(1));
inputCol =  ceil(colStart - offset(2) + eps0):floor((colEnd) - offset(2));

% Check to make sure that the inputRow/Col match the size of the receptive
% field. If we never get this error, then we will delete the check.
if (length(inputRow) ~= size(rgcMosaic.sRFcenter{row,col},1)) || ...
        (length(inputCol)~=size(rgcMosaic.sRFcenter{row,col},1))
    error('Dimension mismatch of input and rf size');
end

end