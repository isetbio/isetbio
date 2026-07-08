function coneRectW = coneRectWindow(cones,varargin)
% Wrapper that replaces the GUIDE coneMosaicWindow
%
% Synopsis
%   coneRectW = coneRectWindow(cones,varargin)
%
% Brief description
%   Opens a coneRectWindow interface based on the coneRectWindow_App.
%
% Inputs
%   cones:  The coneMosaicRect you want in the window.
%
% Outputs
%   coneRectW:  A coneRectWindow_App object.
%
% Optional key/value pairs
%   'plot type':  Initial display. Default is 'mean absorptions'.
%
% See also
%    sceneWindow_App, oiWindow_App
%

% Examples
%{
   cm = coneMosaicRect;
   coneRectWindow(cm);   
%}

%%
varargin = ieParamFormat(varargin);

p = inputParser;
p.addRequired('cones',@(x)(isa(x,'coneMosaicRect')));
p.addParameter('plottype','meanabsorptions',@ischar);

p.parse(cones,varargin{:});

%% See if there is a window.

% Empty, so create one and put it in the vcSESSION
coneRectW = coneRectWindow_App(cones,p.Results.plottype);

drawnow;

end
