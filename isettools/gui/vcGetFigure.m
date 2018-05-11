function figNum = vcGetFigure(obj)
% Return the figure number associated with an object.
%
% Syntax:
%   figNum = vcGetFigure(obj)
%
% Description:
%    Possible ISET objects are scene, opticalimage, isa, vcimage.  This
%    routine allows us to get the figure number when the object type is not
%    identified in the main routine, such as getting an ROI for an oi,
%    scene, sensor or other type of window.
%
%    There is a separate routine for GraphWins.  But I am not sure why.
%
%    The function contains examples of usage below. To access these
%    examples, type 'edit vcGetFigure.m' into the Command Window.
%
% Inputs:
%    obj    - Object. The object to return the figure number of.
%
% Outputs:
%    figNum - Numeric. The figure number associated with the object.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/03/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - Example broken
    figNum = vcGetFigure(obj)
    figure(figNum);
    handles = guihandle(figNum);
%}

% global vcSESSION;

objType = vcGetObjectType(obj);
objType = vcEquivalentObjtype(objType);

switch lower(objType)
    case 'scene'
        figNum = ieSessionGet('scenewindow');
    case {'opticalimage'}
        figNum = ieSessionGet('oiwindow');
    case {'isa'}
        figNum = ieSessionGet('sensorwindow');
    case {'vcimage'}
        figNum = ieSessionGet('vcimagewindow');

    otherwise
        error('Unknown object type.');
end

return;