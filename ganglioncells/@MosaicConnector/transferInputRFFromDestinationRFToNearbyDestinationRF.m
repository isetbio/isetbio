function transferInputRFFromDestinationRFToNearbyDestinationRF(obj, ...
    indexOfInputRFToBeReassigned, destinationRFindex, nearbyDestinationRFindex)
%
% Update the connectivityMatrix, by disconnecting
%   indexOfInputRFToBeReassigned  FROM  destinationRFindex
% and connecting 
%   indexOfInputRFToBeReassigned  to the nearbyDestinationRFindex
%
% USAGE:
%           obj.transferInputRFFromDestinationRFToNearbyDestinationRF( ...
%              indexOfInputRFToBeReassigned, destinationRFindex, nearbyDestinationRFindex);
%

    % DISCONNECT input RF from its current destination RF
    if (obj.connectivityMatrix(indexOfInputRFToBeReassigned, destinationRFindex) == 1)
        obj.connectivityMatrix(indexOfInputRFToBeReassigned, destinationRFindex) = 0; % disconnect
    else
        error('Source RF %d was not connected to destination RF %d\n', indexOfInputRFToBeReassigned, destinationRFindex);
    end
    
    % And CONNECT it to the nearby destination RF
    obj.connectivityMatrix(indexOfInputRFToBeReassigned, nearbyDestinationRFindex) = 1;
    
    % Update the centroids of the 2 destination RFs
    destinationRFList = [destinationRFindex nearbyDestinationRFindex];
    obj.updateDestinationCentroidsFromInputs(destinationRFList);
end