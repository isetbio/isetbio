function removeConesWithinOpticNerveHead(obj)
% Remove cones located within the optic disk
%
% Syntax:
%   obj.removeConesWithinOpticNerveHead()
%
% Description:
%    Remove cones located within the optic disk
%
% Inputs:
%    obj                 - A @cMosaic object
%
% Outputs:                 None

    % Generate OD struct in microns
    [odStructMicrons, ~] = obj.odStruct();


    % Find indices of cones lying inside the optic disk ellipsoid
    idxInside = obj.indicesOfConesWithinROI(odStructMicrons);
    
    % Find indices of cones outside the optic disk
    idxOutsideOpticDisk = setdiff(1:size(obj.coneRFpositionsDegs,1), idxInside);

    % Update state
    obj.updateStateGivenKeptConeIndices(idxOutsideOpticDisk);
end
                  