function config = sensorConfig(ISA)
% Compute the spatial configuration of color filters in an entire ISA
%
%    config = sensorConfig(ISA)
%
%  The routine computes an Nx3 matrix containing (row,col) positions and color
%  information about the pixels in a basic sensor color block.  Usually the
%  block is a 2x2 region that defines a Bayer sampling grid.  But, in
%  principle, the block could have a different shape, size, and include
%  various sampling configurations.  
%  
%  The config field contains (row,col,color) values  where(row,col) are 
%  position (meters).  The third column contains character values that
%  summarize the color property of each pixel.  These characters are
%  permissible color filter name values, such as r,g,b.  The full list is
%  defined in the routine sensorColorOrder().  (I am not sure why that
%  routine has that name.)
%
%  The non-Bayer possibilities for the data have not been tested or
%  developed yet.  This is planned for future experiments and
%  implementations. 
%
% Copyright ImagEval Consultants, LLC, 2005

disp('sensorConfig is obsolete')
evalin('caller','mfilename')

end