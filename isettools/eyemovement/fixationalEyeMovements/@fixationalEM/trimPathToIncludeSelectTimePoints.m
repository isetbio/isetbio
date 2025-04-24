function trimPathToIncludeSelectTimePoints(obj, selectTimeIndices)
%  Method to trim the emPath to include only selectTimeIndices
%
% History:
%   5/24/2025  npc      Wrote it

    obj.emPosArcMin = obj.emPosArcMin(:,selectTimeIndices,:);
    obj.emPosMicrons = obj.emPosMicrons(:,selectTimeIndices,:);
    obj.timeAxis = obj.timeAxis(selectTimeIndices);

end
