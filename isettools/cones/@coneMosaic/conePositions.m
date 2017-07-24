function val = conePositions(obj)
%CONEPOSITIONS  Return the cone sample positions in meters
%   val = CONEPOSITIONS(obj)
%
%   For other objects, such as the bipolars or RGCs to have this grid, we
%   only need to store the patternSampleSize and number of rows and cols
%   Then we can always regenerate the spatial pattern support of the cone
%   mosaic when we are in those other objects.
%
%   So when we run bipolar.compute, we attach these four numbers
%     bp.spatialFrame.row, col, sampleSizeRow, sampleSizeCol
%
%   And then we use this type of routine
%     positions = samplePositions([row,col], sampleSize);
%
% [DHB NOTE: THIS HELP TEXT NEEDS SOME HELP. I DONT' UNDERSTAND IT.  WHAT IS THE FORMAT
% OF VAL?  WHAT IS THE ROUTINE SAMPLEPOSITIONS MENTIONED ABOVE?  IT DOES NOT SEEM TO EXIST
% IN ISETBIO AS OF 01/17/17.]

% JRG/BW ISETBIO Team, 2016

x = (1:obj.cols) * obj.patternSampleSize(1); x = x - mean(x);
y = (1:obj.rows) * obj.patternSampleSize(2); y = y - mean(y);
[xx, yy] = meshgrid(x,y);
val(:,:,1) = xx;
val(:,:,2) = yy;

end