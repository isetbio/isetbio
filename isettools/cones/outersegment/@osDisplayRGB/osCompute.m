function obj = osCompute(obj, sceneRGB, varargin)
% osCompute: this method of @osDisplayRGB passes on the cone isomerizations
% (R*) without any temporal filtering. This subclass is intended to be used
% for stimulus-referred retinal ganglion cell models.
%
% Inputs: the osIdentity object and the sceneRGB data (x,y,t,3), where the 
% last index is the RGB channel.
% 
% Outputs: the osIdentity object, with the rgbData property set 
% to the rgb image of the scene. This is for use with stimulus-referred
% models of retinal ganglion cell processing, like the linear, LNP and GLM 
% subclasses of the @rgc object.
% 
% 8/2015 JRG

obj = osSet(obj, 'rgbData', sceneRGB);

end

