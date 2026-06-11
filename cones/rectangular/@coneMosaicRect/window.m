function window(obj, varargin)
% Open the App Designer rectangular cone mosaic window.
%
% Syntax:
%   window(obj, [varargin])
%
% Description:
%    Opens a rectangular cone mosaic GUI via coneRectWindow.
%
% Inputs:
%    obj - The rectangular cone mosaic object.
%
% Outputs:
%   None.
%
% Optional key/value pairs:
%   'plot type' - 'mean absorptions', 'cone mosaic', 'mean photocurrent'
%
%  See Also:
%    coneRectWindow, coneRectWindow_App


%% 
varargin = ieParamFormat(varargin);
p = inputParser;

valid = {'conemosaic','meanabsorptions','meanphotocurrent'};
p.addRequired('obj',@(x)(isa(x,'coneMosaicRect')));
p.addParameter('plottype','meanabsorptions',@(x)(ismember(ieParamFormat(x),valid)));

p.parse(obj,varargin{:});

coneRectWindow(obj,'plot type',p.Results.plottype);

end
