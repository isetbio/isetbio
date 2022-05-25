function transferConeFromSourceRGCToDestinationRGC(obj, ...
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
    
    % Update the centroids of the 2 RGCs
    affectedRGCindices = [sourceRGCIndex destinationRGCindex];
    obj.updateCentroidsFromInputs(affectedRGCindices);
end