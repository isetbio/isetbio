function val = osGet(obj, param)
% Get base class outer segment parameters
% 
% Parameters:
%  'noise flag' - noise-free ('0') or include noise ('1')
%  'timestep'   -  delta t
%  'cone current signal' - cone current as a function of time
%  'patch size'          - diameter of the cone mosaic patch
%  'array size'          - number of cone samples (row,col)
%
% 7/2015 JRG NC DHB

if ~exist('param','var'), error('Parameter required'); end

switch ieParamFormat(param);  % Lower case and remove spaces

    case {'noiseflag'}
        % The turn on or off the cone noise
        val = obj.noiseFlag;
        
    case{'patchsize'}
        % Diameter of the cone mosaic patch
        val = obj.patchSize;
        
    case{'timestep'}
        % Temporal step size
        val = obj.timeStep;
        
    case{'conecurrentsignal'}
        % Time varying current
        val = obj.coneCurrentSignal;
       
    case {'arraysize'}
        % Spatial samples
        sz = size(obj.coneCurrentSignal);
        val = sz(1:2);
    otherwise
        error('Unknown parameter %s\n',param);
       
end

