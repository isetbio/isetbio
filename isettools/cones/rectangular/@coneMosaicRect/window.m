function window(obj, varargin)
% Open a cone mosaic window GUI (older GUIDE style for now).
%
% TODO:  Update to appdesigner
%
% Syntax:
%   window(obj, [varargin])
%
% Description:
%    Oopens a cone mosaic rect window GUI via the coneMosaicWindow
%    function. 
%
% Inputs:
%    obj - The rectangular cone mosaic object.
%
% Outputs:
%   None.
%
% Optional key/value pairs:
%   'show' - 'mean absorptions', 'cone mosaic', 'mean photocurrent'
%
%  See Also:
%    coneMosaicWindow

% History:
%    xx/xx/16  HJ   ISETBIO Team 2016
%    02/19/18  jnm  Formatting

%% 
p = inputParser;

valid = {'conemosaic','meanabsorptions','meanphotocurrent'};
p.addRequired('obj',@(x)(isa(x,'coneMosaicRect')));
p.addParameter('show','meanabsorptions',@(x)(ismember(ieParamFormat(x),valid)));

p.parse(obj,varargin{:});
coneMosaicWindow(obj,p.Results.show);

end
