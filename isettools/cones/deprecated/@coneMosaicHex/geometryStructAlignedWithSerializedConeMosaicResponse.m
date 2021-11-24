function cmStructSerialized = geometryStructAlignedWithSerializedConeMosaicResponse(obj)
% Return a struct with the mosaic geometry where cones are in same order as
% the serialized response
%
% Syntax:
%    cmStruct = geometryStructAlignedWithSerializedConeMosaicResponse(obj)
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

    
    cmStruct = obj.geometryStruct();
    coneLocsHexGrid = obj.coneLocsHexGridAlignedWithSerializedConeMosaicResponse() * 1e6 / obj.micronsPerDegree;
    
    
    % Find correspondence
    [~, idx] = pdist2(coneLocsHexGrid, cmStruct.coneLocs, 'euclidean', 'Smallest', 1);
    
    cmStructSerialized.coneLocs = cmStruct.coneLocs(idx,:);
    cmStructSerialized.coneLocsMicrons =  cmStructSerialized.coneLocs * obj.micronsPerDegree;
    cmStructSerialized.coneApertures = cmStruct.coneApertures(idx);
    
    if (isempty(cmStruct.coneTypes))
        cmStructSerialized.coneTypes = obj.pattern(find(obj.pattern > 1));
    else
        cmStructSerialized.coneTypes = cmStruct.coneTypes(idx);
    end
    
    
    % Delaunayn triangularization
    xHex = cmStructSerialized.coneLocs(:,1);
    yHex = cmStructSerialized.coneLocs(:,2);
    cmStructSerialized.triangles = delaunayn([xHex(:), yHex(:)]);
    
end