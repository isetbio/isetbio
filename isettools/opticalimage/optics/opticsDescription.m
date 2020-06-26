function txt = opticsDescription(optics)
% Text description of the optics properties.
%
% Syntax:
%   txt = opticsDescription(optics)
%
% Description:
%    The text description of the properties of an optics structure.
%
% Inputs:
%    optics - (Optional) Struct. The optics structure. Default is uses
%             vcGetSelectedObject to find an optics struct.
%
% Outputs:
%    txt    - String. The description string.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    03/09/18  jnm  Formatting
%    04/07/18  dhb  Make example run.
%    06/27/19  JNM  Formatting update

% Examples:
%{
   optics = opticsCreate;
   descriptionText = opticsDescription(optics)
%}

if notDefined('optics'), [~, optics] = vcGetSelectedObject('OPTICS'); end
if ~isfield(optics, 'name'), optics.name = 'No name'; end

txt = sprintf('Optics: %s\n', opticsGet(optics, 'name'));
txt = [txt, sprintf('  Num Aper:  %0.2e  \n', opticsGet(optics, 'na'))];
txt = [txt, sprintf('  Aper Area: %0.2e m^2\n', ...
    opticsGet(optics, 'apertureArea'))];
txt = [txt, sprintf('  Aper Diam: %0.2e m\n', ...
    opticsGet(optics, 'apertureDiameter'))];

end
