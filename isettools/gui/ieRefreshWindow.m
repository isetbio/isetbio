function ieRefreshWindow(objType)
% Refresh one of the four types of windows
%
% Syntax:
%   ieRefreshWindow(objType)
%
% Description:
%    Issue a refresh to one of the ISET windows.  This routine is useful
%    when you have changed the data and would like the window to update
%    according to the new data.
%
%    The function contains examples of usage below. To access these
%    examples, type 'edit ieRefreshWindow.m' into the Command Window.
%
% Inputs:
%    objType - String. The type of the object. Acceptable values include:
%              'scene', 'opticalImage', 'isa', and 'vcImage'.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    05/02/18  jnm  Formatting

% Examples:
%{
    scene = sceneCreate;
    vcReplaceObject(scene);
    ieRefreshWindow('scene');

    oi = oiCreate;
    vcReplaceObject(oi);
    ieRefreshWindow('oi');
%}
%{
    % ETTBSkip - isa is not initialized
    vcReplaceObject(isa);
    ieRefreshWindow('isa');

%}
%{
    % ETTBSkip - vcimage is not initialized
    vcReplaceObject(vcimage);
    ieRefreshWindow('vcimage');
%}


if notDefined('objType'), error('You must specify an object type.'); end

objType = vcEquivalentObjtype(objType);

switch lower(objType)
    case {'scene'}, sceneWindow;
    case {'opticalimage'}, oiWindow;
    case {'isa'}, sensorImageWindow;
    case {'vcimage'}, vcimageWindow;
    otherwise, error('Unknown object type');
end

end