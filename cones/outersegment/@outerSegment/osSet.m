function obj = osSet(obj, param, value)
% Sets isetbio outersegment object properties for base class
% 
% Syntax:
%   obj = osSet(obj, param, value)
%
% Description:
%    Set the isetbio outer segment object properties for the base class.
%
%    Examples are contained in the code. To access, type 'edit osSet.m'
%    into the Command Window.
%
% Inputs:
%    obj   - The outer segment object you wish to assign/change a property
%            value of.
%    param - The name of the parameter you wish assign/change the value of.
%            The options include:
%         'noise flag'          -  'random', 'frozen', 'none'.
%         'time step'           -  Temporal size step.
%         'patch size'          -  Diameter of cone mosaic patch.
%         'cone current signal' -  Cone current as a function of time.
%    value - The value to assign to the aforementioned parameter.
%
% Outputs:
%    obj   - The modified outer segment object.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - Cone current signal does not appear to be supported any
%      longer? Should we remove this option?]

% History:
%    08/xx/15  JRG NC DHB  Created
%    02/12/18  jnm         Formatting

% Examples:
%{
    adaptedOS = osCreate;
	adaptedOS = osSet(adaptedOS, 'noise flag', 'none');
%}

%% Check for the number of arguments and create parser object.
if ~exist('param', 'var'), error('Parameter required'); end
if ~exist('value', 'var'), error('Value required'); end

%%
switch ieParamFormat(param)
    case{'noiseflag'}
        obj.noiseFlag = value;
    case{'timestep'}
        obj.timeStep = value;
    case{'patchsize'}
        % Spatial sample spacing
        obj.patchSize = value;
    case{'conecurrentsignal'}
        obj.coneCurrentSignal = value; 
end
