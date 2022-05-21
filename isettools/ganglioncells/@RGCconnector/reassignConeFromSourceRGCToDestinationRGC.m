function reassignConeFromSourceRGCToDestinationRGC(obj, ...
    indexOfConeToBeReassigned, sourceRGCIndex, destinationRGCindex)
%
% Update the connectivityMatrix, by disconnecting
%   indexOfConeToBeReassigned  FROM  rgcIndex
% and connecting 
%   indexOfConeToBeReassigned  to theTargetRGCindex:
%
% USAGE:
%           obj.reassignConeFromSourceRGCToDestinationRGC( ...
%              indexOfConeToBeReassigned, sourceRGCIndex, destinationRGCindex);
%

    % DISCONNECT cone from its RGC
    if (obj.coneConnectivityMatrix(indexOfConeToBeReassigned, sourceRGCIndex) == 1)
        obj.coneConnectivityMatrix(indexOfConeToBeReassigned, sourceRGCIndex) = 0; % disconnect
    else
        error('Cone %d was not connected to RGC %d\n', indexOfConeToBeReassigned, sourceRGCIndex);
    end
    
    % And CONNECT it to the new RGC
    obj.coneConnectivityMatrix(indexOfConeToBeReassigned, destinationRGCindex) = 1;
    
    % Update the position of the source RGC to be the centroid of its new cone inputs
    indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, sourceRGCIndex)) == 1);
    obj.rgcRFpositionsMicrons(sourceRGCIndex,:) = mean(obj.inputConeMosaic.coneRFpositionsMicrons(indicesOfConeInputs,:),1);
    obj.rgcRFpositionsDegs(sourceRGCIndex,:) = mean(obj.inputConeMosaic.coneRFpositionsDegs(indicesOfConeInputs,:),1);
            
    % Update the position of the destination RGC to be the centroid of its new cone inputs
    indicesOfConeInputs = find(squeeze(obj.coneConnectivityMatrix(:, destinationRGCindex)) == 1);
    obj.rgcRFpositionsMicrons(destinationRGCindex,:) = mean(obj.inputConeMosaic.coneRFpositionsMicrons(indicesOfConeInputs,:),1);
    obj.rgcRFpositionsDegs(destinationRGCindex,:) = mean(obj.inputConeMosaic.coneRFpositionsDegs(indicesOfConeInputs,:),1);
end