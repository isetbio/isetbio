function [inputRow, inputCol] = inputPositions(rgcMosaic,row,col,bipolarsPerMicron)
% Calculate the receptive field position of one RGC cell in input frame
%
%  [inputRow, inputCol] = inputPositions(obj,row,col,bipolarsPerMicron)
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

% The length of the inputRow and inputCol has to end up being the same as
% the size of the receptive field row and col

% The RGC center location with respect to the input (bipolar) sampling grid
rgcCenter = bipolarsPerMicron .* rgcMosaic.cellLocation{row,col};

% Get midpoint of RF by taking half of the col size
% We might want to start using continuous functions that doing everything
% with gridded sample.
sRFMidPointRow = (1/2)*size(rgcMosaic.sRFcenter{1,1},1);

% The region where take the input for this cell
%
% rgcCenter indicates position of RGC on the bipolar samples. The row coord
% of the input is the RGC center minus the midpoint size of the RF
% midpoint.
rowStart = (rgcCenter(1) - sRFMidPointRow);
rowEnd   = (rgcCenter(1) + sRFMidPointRow);

% Get midpoint of RF by taking half of the row size
sRFMidPointCol = (1/2)*size(rgcMosaic.sRFcenter{1,1},2);

% stimCenterCoords indicates position of RGC on stimulus image
% The first x coord of the stimulus of interest is the RGC center minus
% the midpoint size of the RF.
colStart = (rgcCenter(2) - sRFMidPointCol);
colEnd   = (rgcCenter(2) + sRFMidPointCol);

% The offset is the distance in bipolar samples to the edge position.  In
% Matlab notation the edge positions are 1 to max of the size.  The offset
% puts the return row,col values in this range.
offset = bipolarsPerMicron .* rgcMosaic.cellLocation{1,1};
offset = floor(offset);

% If there is rounding, keep it in range
inputRow =  (ceil(rowStart):floor(rowEnd)) - offset(1);
inputCol =  (ceil(colStart):floor(colEnd)) - offset(2);

% Now, how do we check and what do we do to make sure that the inputRow/Col
% match the size of the receptive field.

%% Checking stuff
% if length(inputRow)>length(inputCol)
%     inputRow = inputRow(1:length(inputCol));
% end
% if length(inputCol)>length(inputRow)
%     inputCol = inputCol(1:length(inputRow)); 
% end

if length(inputRow)>size(rgcMosaic.sRFcenter{row,col},1) || length(inputCol)>size(rgcMosaic.sRFcenter{row,col},2) 
    inputRow = inputRow(1:size(rgcMosaic.sRFcenter{row,col},1));
    inputCol = inputCol(1:size(rgcMosaic.sRFcenter{row,col},2)); 
end
% if length(stimY)>size(rgcMosaic.sRFcenter{xcell,ycell},2); stimY = stimY(1:size(rgcMosaic.sRFcenter{xcell,ycell},2)); end;

%% Calculate the offset parameter

% if nargout == 3
%     % An offset is sometimes needed because RGC mosaics may be defined with
%     % their center coordinates not at (0,0).
%     offset(1) = ceil(rgcMosaic.cellLocation{1,1}(1));
%     offset(2) = ceil(rgcMosaic.cellLocation{1,1}(2));
% end

end