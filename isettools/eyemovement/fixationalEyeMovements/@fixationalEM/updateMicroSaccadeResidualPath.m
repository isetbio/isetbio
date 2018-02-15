function updateMicroSaccadeResidualPath(obj)
    % Update the microSaccadeResidualPath
    obj.microSaccadeResidualPath = circshift(obj.microSaccadeResidualPath,-1,2);
    obj.microSaccadeResidualPath(:,end) = 0;
end