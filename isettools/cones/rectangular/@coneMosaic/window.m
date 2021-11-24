function window(obj, varargin)
% Opens a cone mosaic window GUI.
%
% Syntax:
%   window(obj, [varargin])
%
% Description:
%    This function opens a cone mosaic window GUI via calling the
%    coneMosaicWindow function.
%
% Inputs:
%    obj - The cone mosaic object.
%
%   Outputs:
%   None.
%
%   Optional key/value pairs:
%   'show' - 'mean absorptions', 'cone mosaic', 'mean photocurrent'
%
%   See Also:
%    coneMosaicWindow

% History:
%    xx/xx/16  HJ   ISETBIO Team 2016
%    02/19/18  jnm  Formatting

%% 
p = inputParser;

valid = {'conemosaic','meanabsorptions','meanphotocurrent'};
p.addRequired('obj',@(x)(isa(x,'coneMosaic')));
p.addParameter('show','meanabsorptions',@(x)(ismember(ieParamFormat(x),valid)));

p.parse(obj,varargin{:});
coneMosaicWindow(obj,p.Results.show);

end