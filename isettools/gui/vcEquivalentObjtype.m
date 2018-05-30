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

% Get rid of case-dependence.  Return must be upper.
objType = upper(objType);

% These are the aliases we use sometimes
if     strcmp(objType,'OI'), objType = 'OPTICALIMAGE';
elseif strcmp(objType,'SENSOR'), objType = 'ISA';
elseif strcmp(objType,'IMGPROC') || ...
       strcmp(objType,'VCI') || ...     % Virtual camera image
       strcmp(objType,'IP')             % Image processor
    objType = 'VCIMAGE';
end

end
