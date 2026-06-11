function initOutputs(obj)
%  Initialize the outputs for a fixationalEM
%
% Syntax:
%   obj.initOutputs()
%   initOutputs(obj)
%
% Description:
%    Initialize the outputs for a fixational eye movement object.
%
% Inputs:
%    obj - Object. A fixationalEM (fixational eye movement) object.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

obj.velocityArcMin = [];
obj.emPosArcMin = [];
obj.emPosMicrons  = [];
obj.emPos = [];

end