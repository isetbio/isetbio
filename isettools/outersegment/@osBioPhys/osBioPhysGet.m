function val = osBioPhysGet(obj, param, varargin)
%Get isetbio outersegment object parameters
% 
% 
% 
% 
% 6/22/15 James Golden

if ~exist('param','var') || isempty(param)
    error('Parameter field required.');
end

val = []; % initialize output variable

param = ieParamFormat(param);  % Lower case and remove spaces
switch lower(param)

%     case {'noiseflag'}
%         
%         val = obj.noiseflag;
        
        
end

