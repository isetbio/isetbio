function objType = vcEquivalentObjtype(objType)
% Translate aliases into the name used in vcSESSION variable
%
% Syntax:
%   objType = vcEquivalentObjtype(objType);
%
% Description:
%    This call translates equivalent variable names using this call, we can
%    have aliases such as OI and SENSOR instead of OPTICALIMAGE and ISA. Or
%    IMGPROC instead of VCIMAGE.
%
%    The official structure names are SCENE, OPTICALIMAGE, ISA, and VCIMAGE
%
% Inputs:
%    objType - String. A string describing the object type.
%
% Outputs:
%    objType - String. A string describing the object type.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    05/02/18  jnm  Formatting

% Get rid of case-dependence
objType = lower(objType);

% These are the aliases we use sometimes
if strcmp(objType, 'oi'), objType = 'opticalimage';
elseif strcmp(objType, 'sensor'), objType = 'isa';
elseif strcmp(objType, 'imgproc') || ...
       strcmp(objType, 'vci') || ...      % Virtual camera image
       strcmp(objType, 'ip')              % Image processor
    objType = 'vcimage';
    % Other translations belong here
end

return
