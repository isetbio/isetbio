function removeRFsWithinOpticNerveHead(obj)
% Remove midget RGCs located within the optic disk
%
% Syntax:
%   obj.removeRFsWithinOpticNerveHead()
%
% Description:
%    Remove midget RGCs located within the optic disk
%
% Inputs:
%    obj                 - A @mRGCmosaic object
%
% Outputs:                 None

    % Generate OD struct in microns
    [odStructMicrons, ~] = obj.inputConeMosaic.odStruct();
    
    % Find indices of cones lying inside the optic disk ellipsoid
    idxInside = obj.indicesOfRFsWithinROI(odStructMicrons);
    
    % Find indices of cones outside the optic disk
    idxOutside = setdiff(1:size(obj.rgcRFpositionsDegs,1), idxInside);

    obj.rgcRFpositionsDegs = obj.rgcRFpositionsDegs(idxOutside,:);
    obj.rgcRFpositionsMicrons = obj.rgcRFpositionsMicrons(idxOutside,:);
    obj.rgcRFspacingsDegs = obj.rgcRFspacingsDegs(idxOutside);
    obj.rgcRFspacingsMicrons = obj.rgcRFspacingsMicrons(idxOutside);
end