function cmStruct = geometryStruct(obj)
% % Return a struct with the mosaic geometry
%
% Syntax:
%    cmStruct = geometryStruct(obj)
%
% Description:
%    Return a struct with the mosaic geometry (cone positions, Delaunay 
%    triangles, and cone aperture sizes)
%       
%
% Inputs:
%    obj - The cone mosaic hex object
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% History:
%    3/3/19  NPC  ISETBIO TEAM, 2019

    % cone positions in degrees
    cmStruct.coneLocs = obj.coneLocsHexGrid() * 1e6 / obj.micronsPerDegree;
    
    % cone aperture in degrees
    cmStruct.coneApertures = obj.computeApertureDiameters()/obj.micronsPerDegree;
    
    % Delaunayn triangularization
    cmStruct.triangles = delaunayn(conePositions);
end
