function val = osGet(obj, param)
% Get base class outer segment parameters
% 
% Syntax:
%    val = osGet(obj, param)
%
% Description:
%    Retrieve the base class outer segment parameters
%
% Inputs:
%	 obj   - The outer segment object
%    param - The parameter you wish to retrieve the value of. Options
%            include the following:
%         'noise flag'          - Cone noise settings. Options are:
%                                 'random', 'frozen', 'none'
%         'timestep'            - Temporal step size. The delta t
%         'cone current signal' - Cone current as a function of time
%         'patch size'          - Diameter of the cone mosaic patch
%         'array size'          - Number of cone samples (row, col)
%
% Outputs:
%    val   - The value of the requested parameter
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - Cone current signal does not appear to be supported any
%      longer? Should we remove this option?]

% History:
%    07/xx/15  JRG NC DHB  Created
%    02/12/18  jnm         Formatting

if ~exist('param', 'var'), error('Parameter required'); end

switch ieParamFormat(param)  % Lower case and remove spaces
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
        error('Unknown parameter %s\n', param);
end