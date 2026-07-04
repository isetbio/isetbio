function obj = osCompute(obj, sceneRGB, varargin)
% Pass the cone isomerizations (R*) without any temporal filtering. 
%
% Syntax:
%   obj = osCompute(obj, sceneRGB, [varargin])
%
% Description:
%    This method of @osDisplayRGB passes on the cone isomerizations (R*)
%    without any temporal filtering. This subclass is intended to be used
%    for stimulus-referred retinal ganglion cell models.
%
% Inputs:
%    obj      - The osIdentity object
%    sceneRGB - The sceneRGB data (x, y, t, 3), where the last index is the
%               RGB channel.
% 
% Outputs:
%    obj      - The modified osIdentity object, with the rgbData property
%               set to the rgb image of the scene. This is for use with
%               stimulus-referred models of retinal ganglion cell
%               processing, like the linear, LNP and GLM  subclasses of the
%               @rgc object.
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/xx/15  JRG  Created
%    02/14/18  jnm  Formatting

obj = osSet(obj, 'rgbData', sceneRGB);

end

