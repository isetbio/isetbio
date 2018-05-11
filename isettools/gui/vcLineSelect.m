function xy = vcLineSelect(obj, objFig)
% Select (x, y) coordinate that determines a line.
%
% Syntax:
%   xy = vcLineSelect(obj, [objFig])
%
% Description:
%    This routine uses getpts() on the ISET window corresponding to the
%    ISET object. Options for legitimate objects are SCENE, OI, SENSOR, and
%    VCIMAGE. A message is displayed in the window asking the user to
%    select points.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcLineSelect.m' into the Command Window.
%
% Inputs:
%    obj    - Object. The object/struct to select the line from.
%    objFig - Figure. The figure associated with the object..
%
% Outputs:
%    xy     - Vector. The selected point in [x, y] coordinate fashion. Each
%             individual coordinate is returned as an integer via use of
%             the round() function. If multiple points are selected, the
%             last selected point is used in the calculations.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: getpts() allows you to specify the axis, not just the figure.
%      We should probably use that addressing. We should also trap the case
%      no points returned ... is that possible?  Or out of range points?
%
% See Also:
%    vcPointSelect, vcLineSelect, vcROISelect, ieGetXYCoords
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/09/18  jnm  Formatting

% Example:
%{
    % ETTBSkip - Broken example, skip it.
    xy = vcLineSelect(vcGetObject('isa'));
%}
%{
    % ETTBSkip - Requires user input, skipping.
    s = sceneCreate;
    sc = vcAddObject(s);
    sceneWindow;
    xy = vcLineSelect(s);
%}

if notDefined('obj')
    error('You must define an object (isa, oi, scene ...)');
end
if notDefined('objFig'), objFig = vcGetFigure(obj); end

% Select points.
hndl = guihandles(objFig);
msg = sprintf('Right click to select one point.');
ieInWindowMessage(msg, hndl);

[x, y] = getpts(objFig);
nPoints = length(x);
if nPoints > 1
    warning('ISET:vcLineSelect1', '%.0f points selected', nPoints);
    xy = [round(x(end - 1)), round(y(end - 1))];
else
    xy = [round(x), round(y)];
end

ieInWindowMessage('', hndl);

end