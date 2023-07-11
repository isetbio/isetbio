function coneRectW = coneRectWindow(cones,varargin)
% Wrapper that replaces the GUIDE coneMosaicWindow
%
% Synopsis
%   coneRectW = coneRectWindow(sensor,varargin)
%
% Brief description
%   Opens a coneRectWindow interface based on the coneRectWindow_App.
%
% Inputs
%   cones:  The coneMosaicRect you want in the window.
%
% Outputs
%   coneRectW:  An coneMosaicWindow_App object.
%
% Description
%
% See also
%    sceneWindow_App, oiWindow_App
%

% Examples
%{
   coneRectWindow;
%}
%{
   cm = coneMosaicRect;
   coneRectWindow(cm);   
%}

%%
varargin = ieParamFormat(varargin);
p = inputParser;
p.addRequired('cones',@(x)(isa(x,'coneMosaicRect')));
p.addParameter('show',true,@islogical);

p.parse(cones,varargin{:});

%% See if there is a window.

% Empty, so create one and put it in the vcSESSION
coneRectW = coneRectWindow_App(cones);
% coneRectW.refresh;

if p.Results.show
    drawnow;
end

end
