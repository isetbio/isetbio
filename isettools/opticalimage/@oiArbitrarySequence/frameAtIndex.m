function oiFrame = frameAtIndex(obj, frameIndex)
% Compute the oi at desired index
%
% Syntax:
%   oiFrame = frameAtIndex(obj, frameIndex)
%
% Description:
%    Compute the oi at the desired frame index.
%
% Inputs:
%    obj          - Object. The oi arbitrary sequence object.
%    frameIndex   - Numeric. The index of the frame
%
% Outputs:
%    oiFrame - The optical image at frame index
%
% Optional key/value pairs:
%    None.
%

% History:
%    06/26/2020 NPC   ISETBIO Team, 2020

    oiFrame = obj.oiList{frameIndex};
end
