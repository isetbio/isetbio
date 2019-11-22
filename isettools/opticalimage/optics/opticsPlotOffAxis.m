function oi = opticsPlotOffAxis(oi, val)
% Plot relative illuminantion (off-axis) fall off function
%
% Syntax:
%   oi = opticsPlotOffAxis(oi, val)
%
% Description:
%    Plot the relative illumination of the off-axis fall-off for the
%    current optical image and optics assumptions.
%
%    The data are stored in the oi{val} if that parameter is sent in.
%
% Inputs:
%    oi  - Struct. An optical image structure.
%    val - (Optional) Struct. A struct to replace the OI with?
%
% Outputs:
%    oi  - Struct. The modified optical image structure.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: Determine if need to remove the replacement line? (L52)
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/09/18  jnm  Formatting. Add check for OI.

if notDefined('oi'), error('OI is required.'); end
if notDefined('val'), val = []; end

optics = oiGet(oi, 'optics');
data = opticsGet(optics, 'cos4th');

if isempty(data)
    % The fall-off data have not been calculated for the current (shift-
    % invariant) optical image. We calculate the off-axis here and we store
    % it in optics/oi in case the calling return wants to keep the result.
    method = opticsGet(optics, 'cos4thfunction');
    if isempty(method), method = 'cos4th'; end

    % Calculating the cos4th scaling factors
    % We might check whether it exists already and only do this if the
    % cos4th slot is empty.
    optics = feval(method, optics, oi);
    oi = oiSet(oi, 'optics', optics);
    data = opticsGet(optics, 'cos4th');

    % Should probably be eliminated.
    if ~isempty(val), vcReplaceObject(oi, val); end
end

figNum = vcSelectFigure('GRAPHWIN');
plotSetUpWindow(figNum);

mesh(data);

xlabel('Row')
ylabel('Col')
zlabel('Relative intensity');
title('Off-axis intensity falloff');
grid on;
udata.x = data;
set(figNum, 'Userdata', udata);

end