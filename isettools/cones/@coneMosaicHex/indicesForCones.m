function [lConeIndices, mConeIndices,sConeIndices]  = indicesForCones(obj)
% Return indices for the L-, M-, and S-cones
%
% Syntax:
%   [lConeIndices, mConeIndices,sConeIndices]  = indicesForCones(obj)
%
% Description:
%    Return the indices for the L-, M-, and S-cones
%
% Inputs:
%    obj        - The cone mosaic hex object
%
% Outputs:
%    lConeIndices - indices of all L-cones
%    mConeIndices - indices of all M-cones
%    sConeIndices - indices of all S-cones
%
% Optional key/value pairs:
%    None.
%

% History:
%    6/11/18  NPC  ISETBIO TEAM, 2018

    nonNullConeIndices = find(obj.pattern > 1);
    lConeIndices = find(obj.pattern(nonNullConeIndices) == 2);
    mConeIndices = find(obj.pattern(nonNullConeIndices) == 3);
    sConeIndices = find(obj.pattern(nonNullConeIndices) == 4);

end

