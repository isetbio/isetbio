function val = get(obj, param, varargin)
% rgcGet: a method of @rgc that gets rgc object 
% parameters using the input parser structure.
% 
%       val = @rgcLayer.get(rgc, param, varargin)
% 
% Inputs: 
%   @rgcLayer.get(property,varargin);
%
% Outputs: 
%   val - of property
% 
% Proprrties:
%   name: type of rgc object, e.g., 'macaque RGC'
%   input: 'cone current' or 'scene RGB', depends on type of outer
%           segment object created.
%   temporalEquivEcc: the temporal equivalent eccentricity, used to 
%             determine the size of spatial receptive fields.   
%   mosaic: contains rgcMosaic objects for the five most common types
%           of RGCs: onParasol, offParasol, onMidget, offMidget,
%           smallBistratified.
% 
% Example:
%   val = rgcGet(rgc1, 'name')
%   val = rgcGet(rgc1, 'input')
% 
% JRG/BW ISETBIO Team, 2015

%% Parse
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...         
        'name',...
        'input',...
        'temporalEquivEcc',...       
        'mosaic',...
        'timing',...
        'spacing',...
        'numbertrials'
    };
p.addRequired('param',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));

% Parse and put results into structure p.
p.parse(param,varargin{:});
param = ieParamFormat(p.Results.param);

%% Set key-value pairs.
switch param   
    case{'name'}
        val = obj.name;
    case{'input'}
        val = obj.input;
    case{'temporalequivecc'}        
        val = obj.temporalEquivEcc;
    case{'mosaic'}        
        val = obj.mosaic;    
    case{'timing'}
        val = obj.timing;
    case{'spacing'}
        val = obj.spacing;
    case {'numbertrials'}
        val = obj.nTrials;
end

end
