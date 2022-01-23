function [lConeIndices, mConeIndices,sConeIndices, nonNullConeIndices] = indicesForCones(obj)
% Return indices for the L-, M-, and S-cones
%
% Syntax:
%   [lConeIndices, mConeIndices,sConeIndices, nonNullConeIndices] = ...
%     indicesForCones(obj)
%
% Description:
%    Return the indices for the L-, M-, and S-cones within the nonNullCone
%    indices and the nonNullConeIndices
%
% Inputs:
%    obj        - The cone mosaic hex object
%
% Outputs:
%    lConeIndices - indices of all L-cones within the nonNullConeIndices
%    mConeIndices - indices of all M-cones within the nonNullConeIndices
%    sConeIndices - indices of all S-cones within the nonNullConeIndices
%    nonNullConeIndices - indices of all non-null cones
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

