function obj = piWRS(SE,varargin)
% Calls the sceneEye render and then shows the returned object
%
% Synopsis
%   obj = sceneEye.pWRS(varargin);
%
% Brief
%   Writes, Renders, and Shows the sceneEye recipe
%
% Inputs
%   obj - sceneEye object
%
% Key/val pairs
%    render type      - Usual recipe render type cell array 
%    scaleIlluminance - Scale the returned illuminance with respect to
%                       the pupil area. Default true. 
%    write - Call piWrite first. Default: true 
%            For debugging we sometimes suppress piWrite. 
%    docker wrapper - Specify a docker wrapper.  If not specified, we use
%                         dockerWrapper.humanEyeDocker
%   'name'  - Set the Scene or OI name
%   'gamma' - Set the display gamma for the window
%   'show'  - Default: true
%
% Outputs
%   obj - Returns the ISETBio/Cam scene or oi struct
%
% See also
%   sceneEye.render, sceneEye.write
%

%% Parse
varargin = ieParamFormat(varargin);
p = inputParser;

p.addRequired('SE',@(x)(isa(x,'sceneEye')));

p.addParameter('dockerwrapper',dockerWrapper.humanEyeDocker,@(x)(isa(x,'dockerWrapper')));
p.addParameter('name','',@ischar);
p.addParameter('show',true,@islogical);
p.addParameter('gamma',[],@isnumeric);
p.addParameter('renderflag','',@ischar);

p.parse(SE,varargin{:});

thisDocker  = p.Results.dockerwrapper;
g           = p.Results.gamma;
name       = p.Results.name;
show       = p.Results.show;
renderFlag = p.Results.renderflag;

%% Render

obj = SE.render('dockerwrapper',thisDocker);

%% Edit parameters

switch lower(obj.type)
    case 'scene'

        if ~isempty(name), obj = sceneSet(obj,'name',name); end
        if show
            sceneWindow(obj);
            if ~isempty(g), sceneSet(obj,'gamma',g); end
            if ~isempty(renderFlag), sceneSet(obj,'render flag',renderFlag); end
        end

    case 'opticalimage'

        if ~isempty(name), obj = oiSet(obj,'name',name); end
        if show
            oiWindow(obj);
            if ~isempty(g), oiSet(obj,'gamma',g); end
            if ~isempty(renderFlag), oiSet(obj,'render flag',renderFlag); end
        end

    otherwise
        error('Unknown ISET Object type %s',obj.type)
end

end
