function updateMicroSaccadeResidualPath(obj)
% Update the microSaccadeResidualPath
%
% Syntax:
%   updateMicroSaccadeResidualPath(obj)
%
% Description:
%    Update the microSaccadeResidualPath by removing the first sample and
%    shiting the rest by 1.
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
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments


    obj.microSaccadeResidualPath = circshift(...
        obj.microSaccadeResidualPath, -1, 2);
    obj.microSaccadeResidualPath(:,end) = 0;
end