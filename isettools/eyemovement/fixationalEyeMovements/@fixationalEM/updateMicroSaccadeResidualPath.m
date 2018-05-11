function updateMicroSaccadeResidualPath(obj)
% Update the microSaccadeResidualPath
%
% Syntax:
%   updateMicroSaccadeResidualPath(obj)
%
% Description:
%    Update the microSaccadeResidualPath
%
% Inputs:
%    obj - The fixational eye movement object.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
    obj.microSaccadeResidualPath = circshift(...
        obj.microSaccadeResidualPath, -1, 2);
    obj.microSaccadeResidualPath(:,end) = 0;
end